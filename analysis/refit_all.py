"""Refit every hardware memory experiment of the team with ONE fitter.

The three projects used different fit models (amplitude fixed at 1, or free up
to 2), so their published per-round rates are not directly comparable. This
script takes the logical-failure COUNTS per round count from each project and
fits them all with pipeline/qec/analysis.fit_ler:

    P_L(n) = [1 - (1 - 2 eps)^n] / 2           ("all rounds")
    same, rounds >= 2 only                      ("from r2", removes the start effect)

All rows use BP+OSD with an SI1000 (p = 1 %) decoder prior. The prior is each
project's own SI1000 implementation (pipeline: qec.NoiseModel.si1000;
circuit_optimization: SI1000 strengths adapted to its circuit; qrisp_track:
Stim generated-circuit approximation).

Our own runs are re-decoded from the raw shots; the others are read from the
counts saved in their result folders.

    pipeline/.venv/Scripts/python analysis/refit_all.py      (from the repo root)
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "pipeline"))
sys.path.insert(0, str(ROOT / "pipeline" / "experiments"))

from qec.analysis import fit_ler  # noqa: E402

CO = ROOT / "circuit_optimization" / "results"
QT = ROOT / "qrisp_track" / "results_qrisp"
ROUNDS = [1, 2, 3, 4, 5]


def pipeline_counts(fname, prior):
    """Re-decode one of our runs (BP+OSD) and return failures per round count."""
    from characterized_decoding import HW, priors
    from qec import build_memory_circuit
    from qec.decode import decode_dem, measurements_to_detections
    from qec.hardware import load_run

    run = load_run(HW / fname)
    basis, pri = run["basis"], priors(run)
    k, n = [], []
    for r in ROUNDS:
        det, obs = measurements_to_detections(run[f"raw_r{r}"], build_memory_circuit(3, r, basis, True, None))
        dem = build_memory_circuit(3, r, basis, True, pri[prior]).detector_error_model(decompose_errors=True)
        k.append(decode_dem(det, obs, dem, "bposd").errors)
        n.append(len(det))
    return k, n


def twirl_counts(folder, pipeline, arm):
    rows = json.loads((folder / pipeline / "summary.json").read_text())["rows"]
    rows = sorted(rows, key=lambda x: x["rounds"])
    return [x[arm] for x in rows]


def rerun_counts(pipeline, arm):
    base = CO / "twirl_emerald_20261122_rerun"
    k = np.zeros(5, int)
    for c in sorted(base.glob("credential_*")):
        k += twirl_counts(c, pipeline, arm)
    return list(map(int, k)), [16000] * 5


def amine(path, *keys, shots):
    d = json.loads(path.read_text())
    for key in keys:
        d = d[key]
    return [round(p * shots) for p in d], [shots] * 5


def main():
    first = CO / "twirl_emerald_20261122"
    runs = [
        # (circuit, qubits, label, counts, shots)
        ("pipeline (Aziz)", "Aziz patch", "Z, DD on, 21 Nov 15:19", *pipeline_counts("20261121-151940_emerald_Z_ddon.npz", "si1000")),
        ("pipeline (Aziz)", "Aziz patch", "Z, DD on, 21 Nov 21:08 (cal. ecec1d0a)", *pipeline_counts("20261121-210819_emerald_Z_ddon.npz", "si1000")),
        ("pipeline (Aziz)", "Aziz patch", "X, DD on, 21 Nov", *pipeline_counts("20261121-160049_emerald_X_ddon.npz", "si1000")),
        ("circuit_optimization (Ayoub)", "Aziz patch", "Qiskit control, run 1 (cal. ecec1d0a)", twirl_counts(first, "qiskit", "control"), [4000] * 5),
        ("circuit_optimization (Ayoub)", "Aziz patch", "Qrisp control, run 1", twirl_counts(first, "qrisp", "control"), [4000] * 5),
        ("circuit_optimization (Ayoub)", "Aziz patch", "Qiskit control, rerun x4 (cal. ecec1d0a)", *rerun_counts("qiskit", "control")),
        ("circuit_optimization (Ayoub)", "Aziz patch", "Qrisp control, rerun x4", *rerun_counts("qrisp", "control")),
        ("circuit_optimization (Ayoub)", "Aziz patch", "Qiskit twirled, rerun x4", *rerun_counts("qiskit", "twirl")),
        ("circuit_optimization (Ayoub)", "Aziz patch", "Qrisp twirled, rerun x4", *rerun_counts("qrisp", "twirl")),
        ("qrisp_track (Amine)", "own patch (auto)", "SI1000 re-decode", *amine(QT / "cmp_mine_si1000_perround.json", "mine_si1000", "bposd", shots=4000)),
        ("qrisp_track (Amine)", "own patch (auto)", "twirl study, control", *amine(QT / "cmp_twirl_summary.json", "bposd", "notwirl_pl", shots=2000)),
        ("qrisp_track (Amine)", "own patch (auto)", "twirl study, twirled (extra gates)", *amine(QT / "cmp_twirl_summary.json", "bposd", "twirl_pl", shots=2000)),
    ]
    out = []
    print(f"| circuit | qubits | run | shots/pt | P_L(1..5) % | eps all rounds | eps from r2 |")
    print("|---|---|---|---|---|---|---|")
    for circ, qubits, label, k, n in runs:
        a = fit_ler(ROUNDS, k, n)
        b = fit_ler(ROUNDS, k, n, min_round=2)
        pl = [round(100 * x / y, 2) for x, y in zip(k, n)]
        out.append(dict(circuit=circ, qubits=qubits, run=label, failures=list(map(int, k)), shots=list(map(int, n)),
                        eps=a.eps, eps_err=a.eps_err, eps_from_r2=b.eps, eps_from_r2_err=b.eps_err))
        print(f"| {circ} | {qubits} | {label} | {n[0]} | {' / '.join(map(str, pl))} | "
              f"{100 * a.eps:.2f} ± {100 * a.eps_err:.2f} % | {100 * b.eps:.2f} ± {100 * b.eps_err:.2f} % |")
    path = ROOT / "analysis" / "refit_all.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"\nsaved {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
