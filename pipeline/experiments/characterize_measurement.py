"""Repeated mid-circuit measurement on every qubit (noise characterization).

Why: our learned time-like rates show some ancillas err 10-30x more than the
calibration says. In the no-reset surface code each ancilla is measured once
per round and never reset, so what matters is not the single-shot readout
fidelity Resonance reports but how a qubit behaves under *repeated*
measurement: does a measurement flip the state it reports (non-QND,
measurement-induced transitions), and how much does |1> decay during the
measurement window?

Circuits (all qubits in parallel, as in the QEC experiment):
  prep0: measure K times                     (qubit should stay in |0>)
  prep1: X, then measure K times             (qubit should stay in |1>)

Per qubit we report
  e0 / e1     first-measurement readout error for |0> / |1>
  up          P(m_k = 1 | m_{k-1} = 0) in prep0   - excitation under measurement
  down        P(m_k = 0 | m_{k-1} = 1) in prep1   - decay / relaxation under measurement
  flip        mean probability that two consecutive outcomes differ, averaged
              over both preparations - the rate a no-reset ancilla would show
              as a spurious time-like detection event even with perfect gates
  state_flip  P(pattern A A B B in four consecutive outcomes): a change that the
              previous AND the next measurement confirm. A single readout error
              cannot produce it (it needs two), so this isolates real state
              changes caused by (or during) measurement. Needs k >= 4.

    .venv/Scripts/python experiments/characterize_measurement.py --quantum-computer emerald            # dry run
    .venv/Scripts/python experiments/characterize_measurement.py --quantum-computer emerald --confirm
    .venv/Scripts/python experiments/characterize_measurement.py --analyze results/hw/meas_<...>.npz
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec.hardware import (DATA_DIR, get_backend, job_to_arrays, load_run,  # noqa: E402
                          qubit_names, save_run, submit, timestamp)

RESULTS = Path(__file__).resolve().parent.parent / "results"


def circuits(n: int, k: int) -> list[QuantumCircuit]:
    out = []
    for prep in (0, 1):
        qc = QuantumCircuit(n, n * k, name=f"prep{prep}")
        if prep:
            qc.x(range(n))
        for rep in range(k):
            qc.barrier()
            for q in range(n):
                qc.measure(q, rep * n + q)      # clbit = rep * n + qubit
        out.append(qc)
    return out


def analyze(run: dict) -> list[dict]:
    n, k = int(run["n_qubits"]), int(run["k"])
    names = run["names"]
    rows = []
    m0 = run["raw_prep0"].reshape(-1, k, n)
    m1 = run["raw_prep1"].reshape(-1, k, n)
    for q in range(n):
        a, b = m0[:, :, q], m1[:, :, q]
        e0 = a[:, 0].mean()
        e1 = 1 - b[:, 0].mean()
        prev0, nxt0 = a[:, :-1].ravel(), a[:, 1:].ravel()
        prev1, nxt1 = b[:, :-1].ravel(), b[:, 1:].ravel()
        up = nxt0[~prev0].mean() if (~prev0).any() else np.nan
        down = 1 - nxt1[prev1].mean() if prev1.any() else np.nan
        flip = 0.5 * ((prev0 != nxt0).mean() + (prev1 != nxt1).mean())
        sf = []
        for m in (a, b):
            before = m[:, 1:-2] == m[:, :-3]
            changed = m[:, 2:-1] != m[:, 1:-2]
            after = m[:, 3:] == m[:, 2:-1]
            sf.append((before & changed & after).mean())
        rows.append(dict(qubit=names[str(q)], e0=float(e0), e1=float(e1),
                         up=float(up), down=float(down), flip=float(flip),
                         state_flip=float(np.mean(sf))))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quantum-computer", default="emerald")
    ap.add_argument("--k", type=int, default=5, help="measurements per qubit")
    ap.add_argument("--shots", type=int, default=4000)
    ap.add_argument("--confirm", action="store_true")
    ap.add_argument("--analyze")
    args = ap.parse_args()

    if args.analyze:
        run = load_run(Path(args.analyze))
        rows = analyze(run)
        out = Path(args.analyze).with_suffix(".json")
        out.write_text(json.dumps(rows, indent=1))
        print(f"{'qubit':>6} {'e0':>6} {'e1':>6} {'up':>6} {'down':>6} {'flip':>6} {'state':>6}")
        for r in sorted(rows, key=lambda r: -r["flip"]):
            print(f"{r['qubit']:>6} {100*r['e0']:5.1f}% {100*r['e1']:5.1f}% {100*r['up']:5.1f}% "
                  f"{100*r['down']:5.1f}% {100*r['flip']:5.1f}% {100*r['state_flip']:5.1f}%")
        print(f"saved {out}")
        return

    backend = get_backend(args.quantum_computer)
    n = backend.num_qubits
    names = qubit_names(backend)
    qcs = circuits(n, args.k)
    tqcs = [transpile(qc, backend=backend, initial_layout=list(range(n)),
                      layout_method="trivial", routing_method="none", optimization_level=0)
            for qc in qcs]
    print(f"{args.quantum_computer}: {n} qubits x {args.k} measurements, 2 circuits x "
          f"{args.shots} shots = {2 * args.shots} shots")
    if not args.confirm:
        print("dry run only - add --confirm to submit")
        return

    stamp = timestamp()
    job = submit(backend, tqcs, args.shots, dd=False)
    print(f"submitted job {job.job_id()}", flush=True)
    pending = RESULTS / "hw" / f"pending_meas_{stamp}_{args.quantum_computer}.json"
    pending.parent.mkdir(parents=True, exist_ok=True)
    pending.write_text(json.dumps({"job_id": str(job.job_id()), "n": n, "k": args.k}))
    raw0, raw1 = job_to_arrays(job, [n * args.k, n * args.k], timeout=6 * 3600)
    path = RESULTS / "hw" / f"meas_{stamp}_{args.quantum_computer}.npz"
    cal = DATA_DIR / f"calibration_{args.quantum_computer}_latest.json"
    save_run(path, quantum_computer=args.quantum_computer, n_qubits=n, k=args.k,
             shots=args.shots, job_id=str(job.job_id()),
             names={str(k): v for k, v in names.items()},
             calibration_file=str(cal) if cal.exists() else None,
             raw_prep0=raw0, raw_prep1=raw1)
    print(f"saved {path}")
    for r in sorted(analyze(load_run(path)), key=lambda r: -r["flip"])[:10]:
        print(r)


if __name__ == "__main__":
    main()
