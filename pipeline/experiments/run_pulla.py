"""Pulse-level execution of the memory experiment through IQM PulLA.

The same transpiled circuits as run_hardware.py are converted to IQM pulse
circuits, compiled with PulLA's standard compiler (using the current
calibration set), and submitted as pulse playlists. Pulse-level access lets us
set compiler stages that the circuit API does not expose, here a custom
dynamical-decoupling strategy:

  --dd standard : IQM's STANDARD_DD_STRATEGY (same as the server-side DD)
  --dd dense    : up to 3 back-to-back EDD sequences (24 pulses), centre-
                  aligned, in the ~1 us mid-circuit readout window in which data
                  qubits idle (standard DD places a single 8-pulse EDD there)
  --dd off      : no DD

Results are saved in the run_hardware.py format, so analyze.py and
decoder_comparison.py work unchanged.

    .venv/Scripts/python experiments/run_pulla.py --dd dense --rounds 1 2 3 4 5 --shots 1000            # compile only
    .venv/Scripts/python experiments/run_pulla.py --dd dense --rounds 1 2 3 4 5 --shots 1000 --confirm
"""

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import build_memory_circuit  # noqa: E402
from qec.hardware import (DATA_DIR, SERVER_URL, coupling_edges, get_backend,  # noqa: E402
                          memory_to_array, save_run, stim_to_qiskit, swap_report,
                          timestamp, transpile_pinned)

RESULTS = Path(__file__).resolve().parent.parent / "results"


def dd_strategy(kind: str):
    from iqm.station_control.interface.models.circuit import DDStrategy
    if kind == "standard":
        from iqm.cpc.compiler.dd import STANDARD_DD_STRATEGY
        return STANDARD_DD_STRATEGY
    if kind == "dense":
        edd = "XYXYYXYX"
        return DDStrategy(gate_sequences=[(25, edd * 3, "center"),
                                          (17, edd * 2, "center"),
                                          (9, edd, "center"),
                                          (5, "YXYX", "center"),
                                          (2, "XX", "center")])
    return None


def compile_for_pulla(pulla, backend, tqcs, shots: int, dd: str):
    from iqm.pulla.utils_qiskit import qiskit_to_pulla

    circuits, compiler = qiskit_to_pulla(pulla, backend, tqcs)
    settings = compiler.get_settings(circuits=circuits)
    settings.set_shots(shots)
    settings.stages.dynamical_decoupling.apply_dd_strategy.dd_is_disabled = dd == "off"
    context = {"DDStrategy": dd_strategy(dd)} if dd != "off" else {}
    run_def, context = compiler.compile(circuits, settings=settings, context=context)
    return run_def, context


def drive_pulse_count(run_def, qubit: str, segment: int = -1) -> int:
    """Number of physical drive pulses on one qubit in one circuit (shows DD)."""
    pl = run_def.sweep_definition.playlist
    ch = f"{qubit}__drive.awg"
    cd = pl.channel_descriptions[ch]
    idx = pl.segments[segment].instructions[ch]
    return sum(type(cd.instruction_table[i].operation).__name__ == "IQPulse" for i in idx)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--basis", default="Z", choices=["Z", "X"])
    ap.add_argument("--shots", type=int, default=1000)
    ap.add_argument("--dd", default="standard", choices=["off", "standard", "dense"])
    ap.add_argument("--quantum-computer", default="emerald")
    ap.add_argument("--layout")
    ap.add_argument("--confirm", action="store_true")
    ap.add_argument("--wait", type=float, default=6 * 3600)
    args = ap.parse_args()

    from iqm.pulla.pulla import Pulla
    from iqm.pulla.utils_qiskit import sweep_job_to_qiskit

    layout_path = Path(args.layout or RESULTS / f"layout_{args.quantum_computer}.json")
    lay_info = json.loads(layout_path.read_text())
    layout = {int(k): v for k, v in lay_info["layout"].items()}
    backend = get_backend(args.quantum_computer)
    pulla = Pulla(SERVER_URL, quantum_computer=args.quantum_computer)

    circuits, tqcs = [], []
    for r in args.rounds:
        c = build_memory_circuit(3, r, args.basis, True, None)
        tqc = transpile_pinned(stim_to_qiskit(c), backend, layout)
        rep = swap_report(tqc, coupling_edges(backend))
        assert rep["swaps"] == 0 and rep["non_native_pairs"] == 0
        circuits.append(c)
        tqcs.append(tqc)

    run_def, context = compile_for_pulla(pulla, backend, tqcs, args.shots, args.dd)
    data_q = lay_info["physical"]["4"]          # centre data qubit
    print(f"compiled {len(tqcs)} circuits with DD={args.dd}: centre data qubit {data_q} "
          f"gets {drive_pulse_count(run_def, data_q)} drive pulses in the {args.rounds[-1]}-round circuit")
    if not args.confirm:
        print("compile only - add --confirm to submit")
        return

    job = pulla.submit_playlist(run_def, context=context)
    stamp = timestamp()
    print(f"submitted pulse job {job.job_id}", flush=True)
    pending = RESULTS / "hw" / f"pending_pulla_{stamp}_{args.quantum_computer}_{args.basis}_{args.dd}.json"
    pending.write_text(json.dumps({"job_id": str(job.job_id), "rounds": args.rounds,
                                   "dd": args.dd, "shots": args.shots}))

    t0 = time.time()
    while True:
        job.update()
        if str(job.status).lower().endswith("completed"):
            break
        if any(s in str(job.status).lower() for s in ("failed", "cancel", "abort")):
            raise SystemExit(f"job ended with status {job.status}")
        if time.time() - t0 > args.wait:
            raise SystemExit(f"still {job.status} after {args.wait}s; job id in {pending}")
        time.sleep(15)

    result = sweep_job_to_qiskit(job, shots=args.shots)
    arrays = [memory_to_array(result.get_memory(i), c.num_measurements)
              for i, c in enumerate(circuits)]
    cal = json.loads((DATA_DIR / f"calibration_{args.quantum_computer}_latest.json").read_text())
    cal_file = DATA_DIR / f"calibration_{args.quantum_computer}_run{stamp}.json"
    cal_file.write_text(json.dumps(cal))
    path = RESULTS / "hw" / f"{stamp}_{args.quantum_computer}_{args.basis}_pulla-{args.dd}.npz"
    save_run(path, rounds=args.rounds, basis=args.basis, shots=args.shots,
             dd=args.dd != "off", pulse_level=True, dd_strategy=args.dd,
             job_id=str(job.job_id), quantum_computer=args.quantum_computer, fake=False,
             layout=lay_info["layout"], physical=lay_info["physical"],
             calibration_file=str(cal_file), depths=[t.depth() for t in tqcs],
             **{f"raw_r{r}": a for r, a in zip(args.rounds, arrays)})
    print(f"saved {path}")


if __name__ == "__main__":
    main()
