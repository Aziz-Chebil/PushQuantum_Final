"""Choose the physical qubits: SWAP-free and calibration-optimal.

1. Enumerate every SWAP-free embedding of the d=3 rotated code into the
   coupler graph (VF2 subgraph monomorphism).
2. Score each with a Stim simulation under a NoiseModel built from the live
   calibration of exactly those qubits (role-aware: data vs X/Z ancilla).
3. Transpile the hardware circuit with the winner pinned and routing disabled,
   and report SWAP count / non-native pairs (both must be 0).

    .venv/Scripts/python experiments/layout.py            # live (needs token)
    .venv/Scripts/python experiments/layout.py --offline  # snapshot + fake backend
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import RotatedSurfaceCode, build_memory_circuit  # noqa: E402
from qec.hardware import (calibration_from_fake, coupling_edges,  # noqa: E402
                          fetch_calibration, get_backend, load_calibration, qubit_names, stim_to_qiskit,
                          swap_report, transpile_pinned)
from qec.layout import draw_layout, rank_layouts  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true",
                    help="use the coupling snapshot, cached calibration and a fake backend")
    ap.add_argument("--quantum-computer", default="emerald")
    ap.add_argument("--exclude", nargs="*", default=[], help="physical qubits to avoid, e.g. QB25")
    ap.add_argument("--rounds", type=int, default=3, help="rounds used for scoring")
    ap.add_argument("--shots", type=int, default=20_000)
    ap.add_argument("--top-k", type=int, default=64,
                    help="how many placements (best by proxy cost) to simulate")
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)

    if args.offline:
        from iqm.qiskit_iqm.fake_backends import IQMFakeAphrodite
        backend = IQMFakeAphrodite()          # same 54-qubit lattice family
        edges = coupling_edges(None, args.quantum_computer)
        names = qubit_names()
        cal = load_calibration(args.quantum_computer)
        cal_source = "cached"
        if cal is None:
            # No cached Emerald calibration: use the fake backend's own error
            # profile, so the ranking is exercised against the noise that
            # run_hardware.py --fake will actually simulate.
            cal = calibration_from_fake(backend)
            cal_source = "fake-backend error profile"
    else:
        backend = get_backend(args.quantum_computer)
        edges = coupling_edges(backend)
        names = qubit_names(backend)
        cal = fetch_calibration(backend)
        cal_source = "live"

    name_to_idx = {v: k for k, v in names.items()}
    exclude = [name_to_idx[n] for n in args.exclude]
    code = RotatedSurfaceCode(3)
    ranked = rank_layouts(code, edges, cal, names, rounds=args.rounds,
                          shots=args.shots, exclude=exclude, top_k=args.top_k)
    print(f"{len(ranked)} SWAP-free placements scored (calibration: {cal_source})")
    for i, s in enumerate(ranked[:5]):
        print(f"  #{i}: predicted P_L({args.rounds} rounds) = {s.predicted_p_l:.4f}  "
              f"proxy = {s.proxy_cost:.3f}  qubits = {sorted(s.physical_qubits(names))}")
    if len(ranked) > 1 and ranked[-1].predicted_p_l == ranked[-1].predicted_p_l:
        print(f"  worst: predicted P_L = {ranked[-1].predicted_p_l:.4f}")

    best = ranked[0]
    print("\nChosen placement (stim coordinates -> physical qubit):")
    print(draw_layout(code, best.layout, names))

    circ = build_memory_circuit(3, args.rounds, "Z", True, None)
    tqc = transpile_pinned(stim_to_qiskit(circ), backend, best.layout)
    rep = swap_report(tqc, edges)
    print(f"\nTranspiled (routing disabled): {rep}")
    assert rep["swaps"] == 0 and rep["non_native_pairs"] == 0

    out = {
        "quantum_computer": args.quantum_computer,
        "layout": {str(k): v for k, v in best.layout.items()},
        "physical": {str(k): names[v] for k, v in best.layout.items()},
        "predicted_p_l": best.predicted_p_l,
        "scoring_rounds": args.rounds,
        "calibration_source": cal_source,
        "swap_report": {k: v for k, v in rep.items() if k != "ops"},
        "ranking": [{"qubits": sorted(s.physical_qubits(names)),
                     "predicted_p_l": s.predicted_p_l} for s in ranked],
    }
    out_path = RESULTS / f"layout_{args.quantum_computer}.json"
    out_path.write_text(json.dumps(out, indent=1))
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
