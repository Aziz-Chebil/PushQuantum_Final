"""Simulation baseline: P_L vs rounds and the fitted logical error rate.

Compares, on the *same* no-reset native-CZ circuit that runs on Emerald:
  * noise models: uniform (Stim-style), SI1000, and calibrated (if a
    calibration file exists in data/ - fetch one with experiments/layout.py)
  * reset vs no-reset ancillas
  * Z vs X memory
  * MWPM vs correlated MWPM vs no decoding

    .venv/Scripts/python experiments/simulate.py --p 0.003 --shots 50000
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import NoiseModel, RotatedSurfaceCode, build_memory_circuit  # noqa: E402
from qec.analysis import fit_ler, plot_ler  # noqa: E402
from qec.decode import decode  # noqa: E402
from qec.hardware import load_calibration, qubit_names  # noqa: E402

RESULTS = Path(__file__).resolve().parent.parent / "results"


def sweep(noise, rounds, basis, no_reset, shots, decoders, seed=0, distance=3):
    out = {dec: {"errors": [], "shots": []} for dec in decoders}
    for r in rounds:
        c = build_memory_circuit(distance, r, basis, no_reset, noise)
        det, obs = c.compile_detector_sampler(seed=seed + r).sample(shots, separate_observables=True)
        for dec in decoders:
            res = decode(det, obs, c, dec)
            out[dec]["errors"].append(res.errors)
            out[dec]["shots"].append(res.shots)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=float, default=0.003)
    ap.add_argument("--shots", type=int, default=50_000)
    ap.add_argument("--rounds", type=int, nargs="+", default=[1, 2, 3, 4, 5, 6, 8, 10])
    ap.add_argument("--distance", type=int, default=3)
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)

    models = {"uniform": NoiseModel.uniform(args.p), "si1000": NoiseModel.si1000(args.p)}
    cal = load_calibration()
    layout_file = RESULTS / "layout_emerald.json"
    if cal is not None and layout_file.exists():
        lay = {int(k): v for k, v in json.loads(layout_file.read_text())["layout"].items()}
        names = qubit_names()
        models["calibrated"] = NoiseModel.from_calibration(
            cal, {q: names[p] for q, p in lay.items()})

    decoders = ["mwpm", "correlated", "none"]
    summary = []
    fig, axes = plt.subplots(1, len(models), figsize=(5 * len(models), 4), squeeze=False)
    for ax, (mname, noise) in zip(axes[0], models.items()):
        for (basis, no_reset), color in zip(
                [("Z", True), ("Z", False), ("X", True)], ["C0", "C1", "C2"]):
            res = sweep(noise, args.rounds, basis, no_reset, args.shots, decoders,
                        distance=args.distance)
            for dec in decoders:
                k, n = np.array(res[dec]["errors"]), np.array(res[dec]["shots"])
                fit = fit_ler(args.rounds, k, n)
                fit_a = fit_ler(args.rounds, k, n, free_amplitude=True)
                row = dict(model=mname, basis=basis, no_reset=no_reset, decoder=dec,
                           rounds=args.rounds, p_l=(k / n).tolist(),
                           eps=fit.eps, eps_err=fit.eps_err,
                           eps_freeA=fit_a.eps, eps_freeA_err=fit_a.eps_err, A=fit_a.amp)
                summary.append(row)
                print(f"{mname:>10} {basis} {'no-reset' if no_reset else 'reset   '} "
                      f"{dec:>10}: {fit}   | free A: {fit_a}")
                if dec == "mwpm":
                    p = k / n
                    plot_ler(ax, args.rounds, p, np.sqrt(p * (1 - p) / n), fit, color=color,
                             label=f"{basis}, {'no reset' if no_reset else 'reset'}")
        ax.set_title(f"d={args.distance} sim, {mname} noise" +
                     (f" (p={args.p})" if mname != "calibrated" else ""))
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(RESULTS / f"sim_ler_d{args.distance}.png", dpi=150)
    (RESULTS / f"sim_summary_d{args.distance}.json").write_text(json.dumps(summary, indent=1))
    print(f"saved results/sim_ler_d{args.distance}.png and sim_summary_d{args.distance}.json")


if __name__ == "__main__":
    main()
