"""Decoder priors built from our own noise characterization.

For every hardware run, the same shots are decoded with three priors (the
noise model the decoder's matching graph / BP weights are built from):

  si1000        generic superconducting model, p = 1 %
  calibrated    per-qubit, per-coupler Resonance calibration
  characterized the calibrated model with each qubit's readout replaced by the
                repeated-measurement characterization: a classical misread
                (M(p)) plus a persisting state change during readout
                (NoiseModel.with_characterization)
  si1000+char   SI1000 for gates and idling, the characterization for readout
                (the gate/idle part of the calibrated model was what hurt)

    .venv/Scripts/python experiments/characterized_decoding.py
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import NoiseModel, build_memory_circuit  # noqa: E402
from qec.analysis import fit_ler  # noqa: E402
from qec.decode import decode, measurements_to_detections  # noqa: E402
from qec.hardware import load_run  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
HW = ROOT / "results" / "hw"
CHAR = {"emerald": HW / "meas_20261121-170050_emerald.json",
        "garnet": HW / "meas_20261121-170050_garnet.json"}
RUNS = [
    ("Emerald Z, DD off", "20261121-151940_emerald_Z_ddoff.npz"),
    ("Emerald Z, DD on", "20261121-151940_emerald_Z_ddon.npz"),
    ("Emerald X, DD on", "20261121-160049_emerald_X_ddon.npz"),
    ("Emerald region 1", "20261121-170221_emerald_Z_ddon.npz"),
    ("Emerald region 2", "20261121-170228_emerald_Z_ddon.npz"),
    ("Emerald region 3", "20261121-170232_emerald_Z_ddon.npz"),
    ("Emerald, dense DD", "20261121-171237_emerald_Z_ddon-dense.npz"),
    ("Emerald, circuit ref.", "20261121-181854_emerald_Z_ddon.npz"),
    ("Emerald, PulLA standard", "20261121-182431_emerald_Z_pulla-standard.npz"),
    ("Emerald, PulLA dense", "20261121-182440_emerald_Z_pulla-dense.npz"),
    ("Garnet Z, DD off", "20261121-144025_garnet_Z_ddoff.npz"),
    ("Garnet Z, DD on", "20261121-144025_garnet_Z_ddon.npz"),
]
DECODERS = ["mwpm", "correlated", "bposd"]
PRIORS = ["si1000", "calibrated", "characterized", "si1000+char"]


def priors(run) -> dict[str, NoiseModel]:
    names = {int(k): v for k, v in run["physical"].items()}
    cal = json.loads(Path(run["calibration_file"]).read_text())["calibration"]
    m_cal = NoiseModel.from_calibration(cal, names)
    rows = json.loads(CHAR[run["quantum_computer"]].read_text())
    m_si = NoiseModel.si1000(0.01)
    return {"si1000": m_si, "calibrated": m_cal,
            "characterized": m_cal.with_characterization(rows, names),
            "si1000+char": m_si.with_characterization(rows, names)}


def main():
    table = []
    print(f"{'run':<26}" + "".join(f"{p + '/' + d:>26}" for p in PRIORS for d in DECODERS))
    for label, fname in RUNS:
        run = load_run(HW / fname)
        rounds, basis = [int(r) for r in run["rounds"]], run["basis"]
        data = {r: measurements_to_detections(run[f"raw_r{r}"], build_memory_circuit(3, r, basis, True, None))
                for r in rounds}
        row = {"run": label, "file": fname}
        cells = []
        for pname, prior in priors(run).items():
            for dec in DECODERS:
                k, n = [], []
                for r in rounds:
                    det, obs = data[r]
                    res = decode(det, obs, build_memory_circuit(3, r, basis, True, prior), dec)
                    k.append(res.errors)
                    n.append(res.shots)
                fit = fit_ler(rounds, k, n)
                row[f"{pname}/{dec}"] = (fit.eps, fit.eps_err)
                cells.append(f"{100 * fit.eps:6.2f} ± {100 * fit.eps_err:.2f} %")
        table.append(row)
        print(f"{label:<26}" + "".join(f"{c:>26}" for c in cells), flush=True)

    # Paired comparison against the generic SI1000 prior (same shots).
    print()
    for dec in DECODERS:
        for p in PRIORS[1:]:
            d = [100 * (r[f"{p}/{dec}"][0] - r[f"si1000/{dec}"][0]) for r in table]
            wins = sum(r[f"{p}/{dec}"][0] < r[f"si1000/{dec}"][0] for r in table)
            best = sum(r[f"{p}/{dec}"][0] == min(r[f"{q}/{dec}"][0] for q in PRIORS) for r in table)
            print(f"{dec:>10} {p:>14} - si1000 = {np.mean(d):+.2f} pp (mean), "
                  f"better in {wins}/{len(table)}, best prior in {best}/{len(table)}")
    out = ROOT / "results" / "figures" / "characterized_decoding.json"
    out.write_text(json.dumps(table, indent=1))
    print(f"saved {out}")


if __name__ == "__main__":
    main()
