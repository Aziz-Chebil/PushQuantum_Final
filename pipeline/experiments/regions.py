"""Pick distinct regions of the chip for a spatial-variation study.

All SWAP-free placements are scored with the live calibration (as in
layout.py), grouped by the set of physical qubits they use, and the best
orientation of each set is kept. Starting from the region already measured
(results/layout_<qc>.json), `--n` further regions are chosen greedily to
minimise qubit overlap with the regions picked so far, so the comparison
samples genuinely different parts of the chip.

Writes results/layout_<qc>_region<k>.json (same format as layout.py), ready for
    run_hardware.py --layout results/layout_<qc>_region<k>.json
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import RotatedSurfaceCode  # noqa: E402
from qec.hardware import coupling_edges, fetch_calibration, get_backend, qubit_names  # noqa: E402
from qec.layout import draw_layout, rank_layouts  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quantum-computer", default="emerald")
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--shots", type=int, default=20_000)
    args = ap.parse_args()

    backend = get_backend(args.quantum_computer)
    edges, names = coupling_edges(backend), qubit_names(backend)
    cal = fetch_calibration(backend)
    code = RotatedSurfaceCode(3)
    ranked = rank_layouts(code, edges, cal, names, rounds=3, shots=args.shots, top_k=1000)

    best_per_set = {}
    for s in ranked:                       # ranked best-first
        key = frozenset(s.layout.values())
        best_per_set.setdefault(key, s)
    print(f"{len(ranked)} placements, {len(best_per_set)} distinct qubit sets")

    current = json.loads((RESULTS / f"layout_{args.quantum_computer}.json").read_text())
    chosen = [frozenset(current["layout"].values())]
    picks = []
    for _ in range(args.n):
        cand = [k for k in best_per_set if k not in chosen]
        # least overlap with everything chosen so far; ties -> better predicted P_L
        k = min(cand, key=lambda k: (max(len(k & c) for c in chosen),
                                     best_per_set[k].predicted_p_l))
        chosen.append(k)
        picks.append(best_per_set[k])

    for i, s in enumerate(picks, 1):
        overlap = len(frozenset(s.layout.values()) & chosen[0])
        print(f"\nregion {i}: predicted P_L(3) = {s.predicted_p_l:.4f}, "
              f"{overlap}/17 qubits shared with the measured patch")
        print(draw_layout(code, s.layout, names))
        out = {"quantum_computer": args.quantum_computer,
               "layout": {str(k): v for k, v in s.layout.items()},
               "physical": {str(k): names[v] for k, v in s.layout.items()},
               "predicted_p_l": s.predicted_p_l, "scoring_rounds": 3,
               "calibration_source": "live", "region": i}
        (RESULTS / f"layout_{args.quantum_computer}_region{i}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
