"""Run the no-reset memory experiment on IQM Emerald for several round counts.

All round counts go into ONE batch job (same calibration for every point of
the P_L(n) curve). Raw per-shot measurements are saved, so every decoding /
noise-model choice can be revisited later without spending QPU time.

    .venv/Scripts/python experiments/run_hardware.py --rounds 1 2 3 4 5 6 --shots 4000            # dry run
    .venv/Scripts/python experiments/run_hardware.py --rounds 1 2 3 4 5 6 --shots 4000 --confirm  # spends credits
    ... --dd both      # DD off and on, back to back, same layout
    ... --fake         # full dry run on IQM's fake 54-qubit backend (noisy Aer sim)
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import build_memory_circuit, describe  # noqa: E402
from qec.hardware import (DATA_DIR, coupling_edges, get_backend, job_to_arrays,  # noqa: E402
                          payload_dict, save_run, stim_to_qiskit, submit, swap_report,
                          timestamp, transpile_pinned)

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, nargs="+", default=[1, 2, 3, 4, 5, 6])
    ap.add_argument("--basis", default="Z", choices=["Z", "X"])
    ap.add_argument("--shots", type=int, default=4000)
    ap.add_argument("--dd", default="off", choices=["off", "on", "both"])
    ap.add_argument("--dd-strategy", default="standard", choices=["standard", "dense"],
                    help="DD sequence used when DD is on (see qec.hardware.dd_strategy)")
    ap.add_argument("--quantum-computer", default="emerald")
    ap.add_argument("--layout", help="default: results/layout_<quantum-computer>.json")
    ap.add_argument("--confirm", action="store_true", help="actually submit (spends credits)")
    ap.add_argument("--fake", action="store_true", help="run on IQMFakeAphrodite instead")
    ap.add_argument("--wait", type=float, default=6 * 3600,
                    help="seconds to wait for results before leaving the job pending")
    ap.add_argument("--retrieve", help="collect a previously submitted run from its pending_*.json")
    args = ap.parse_args()

    if args.retrieve:
        collect(json.loads(Path(args.retrieve).read_text()), None, args.wait)
        return

    layout_path = Path(args.layout or RESULTS / f"layout_{args.quantum_computer}.json")
    lay_info = json.loads(layout_path.read_text())
    if lay_info.get("quantum_computer", args.quantum_computer) != args.quantum_computer and not args.fake:
        raise SystemExit(f"{layout_path} was made for {lay_info['quantum_computer']}, "
                         f"not {args.quantum_computer}: run experiments/layout.py first")
    layout = {int(k): v for k, v in lay_info["layout"].items()}
    print("layout:", lay_info["physical"])

    if args.fake:
        from iqm.qiskit_iqm.fake_backends import IQMFakeAphrodite
        backend = IQMFakeAphrodite()
        edges = coupling_edges(None)
    else:
        backend = get_backend(args.quantum_computer)
        edges = coupling_edges(backend)

    tqcs = []
    for r in args.rounds:
        c = build_memory_circuit(3, r, args.basis, True, None)
        tqc = transpile_pinned(stim_to_qiskit(c), backend, layout)
        rep = swap_report(tqc, edges)
        assert rep["swaps"] == 0 and rep["non_native_pairs"] == 0, rep
        print(f"  rounds={r}: {describe(c)['gate_counts']}  -> transpiled depth {rep['depth']}, "
              f"{rep['two_qubit_gates']} CZ, 0 SWAP")
        tqcs.append(tqc)

    total = len(tqcs) * args.shots * (2 if args.dd == "both" else 1)
    print(f"{len(tqcs)} circuits x {args.shots} shots"
          f"{' x 2 (DD off/on)' if args.dd == 'both' else ''} = {total} shots total")
    if not (args.confirm or args.fake):
        print("dry run only - add --confirm to submit")
        return

    # Freeze the calibration this run was planned with (the *_latest file moves on).
    stamp = timestamp()
    if args.fake:
        from qec.hardware import calibration_from_fake
        cal = calibration_from_fake(backend)
    else:
        cal = json.loads((DATA_DIR / f"calibration_{args.quantum_computer}_latest.json").read_text())["calibration"]
    cal_file = DATA_DIR / f"calibration_{'fake' if args.fake else args.quantum_computer}_run{stamp}.json"
    cal_file.write_text(json.dumps({"fetched_at": stamp, "calibration": cal}))

    # Keep the exact circuits sent to IQM (checked by experiments/verify_circuits.py).
    pay_dir = RESULTS / "hw" / "payloads"
    pay_dir.mkdir(parents=True, exist_ok=True)
    for r, tqc in zip(args.rounds, tqcs):
        try:
            (pay_dir / f"{stamp}_{args.quantum_computer}_{args.basis}_submitted_r{r}.json").write_text(
                json.dumps(payload_dict(backend, tqc)))
        except Exception as exc:  # fake backends may not serialize
            print(f"payload not saved for rounds={r}: {exc}")

    dd_settings = {"off": [False], "on": [True], "both": [False, True]}[args.dd]
    jobs, live = [], {}
    for dd in dd_settings:
        job = backend.run(tqcs, shots=args.shots) if args.fake else submit(backend, tqcs, args.shots, dd=dd,
                                                                  dd_strategy_name=args.dd_strategy)
        print(f"submitted DD={'on' if dd else 'off'}: job {job.job_id()}", flush=True)
        jobs.append({"dd": dd, "job_id": str(job.job_id())})
        live[str(job.job_id())] = job

    pending = dict(stamp=stamp, rounds=args.rounds, basis=args.basis, shots=args.shots,
                   dd_strategy=args.dd_strategy,
                   quantum_computer=args.quantum_computer, fake=args.fake, jobs=jobs,
                   layout=lay_info["layout"], physical=lay_info["physical"],
                   calibration_file=str(cal_file), depths=[t.depth() for t in tqcs])
    pending_path = RESULTS / "hw" / f"pending_{stamp}_{args.quantum_computer}_{args.basis}.json"
    pending_path.parent.mkdir(parents=True, exist_ok=True)
    pending_path.write_text(json.dumps(pending, indent=1))
    print(f"job info saved to {pending_path}", flush=True)
    collect(pending, live, args.wait, pending_path)


def collect(pending: dict, live_jobs: dict | None, wait: float, pending_path: Path | None = None):
    """Wait for (or re-attach to) the jobs of a run and save their raw shots."""
    backend = None
    rounds = pending["rounds"]
    num_meas = [build_memory_circuit(3, r, pending["basis"], True, None).num_measurements
                for r in rounds]
    for entry in pending["jobs"]:
        jid, dd = entry["job_id"], entry["dd"]
        strat = pending.get("dd_strategy", "standard")
        tag = (f"{pending['stamp']}_{'fake' if pending['fake'] else pending['quantum_computer']}"
               f"_{pending['basis']}_dd{'on' if dd else 'off'}"
               f"{'-' + strat if dd and strat != 'standard' else ''}")
        path = RESULTS / "hw" / f"{tag}.npz"
        if path.exists():
            print(f"already saved: {path}")
            continue
        job = (live_jobs or {}).get(jid)
        if job is None:
            backend = backend or get_backend(pending["quantum_computer"])
            job = backend.retrieve_job(jid)
        try:
            arrays = job_to_arrays(job, num_meas, timeout=wait)
        except Exception as exc:
            print(f"job {jid} not finished ({type(exc).__name__}: {exc}). "
                  f"collect later with: experiments/run_hardware.py --retrieve "
                  f"{pending_path or '<pending json>'}")
            continue
        save_run(path, dd=dd, job_id=jid,
                 **{k: v for k, v in pending.items() if k not in ("jobs", "stamp")},
                 **{f"raw_r{r}": a for r, a in zip(rounds, arrays)})
        print(f"saved {path}")


if __name__ == "__main__":
    main()
