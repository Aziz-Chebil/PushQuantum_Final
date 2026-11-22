"""Decoder ablation on identical hardware shots.

Every decoder sees the same detection events and the same prior (DEM); only
the decoding algorithm changes. Reports eps_L per decoder and prior, plus the
per-shot decoding time (Python-loop timing for the ldpc decoders, so it is
an upper bound, not the asymptotic speed).

    .venv/Scripts/python experiments/decoder_comparison.py results/hw/<run>.npz ...
"""

import argparse
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import NoiseModel, build_memory_circuit  # noqa: E402
from qec.analysis import fit_ler  # noqa: E402
from qec.decode import decode, measurements_to_detections  # noqa: E402
from qec.hardware import load_run  # noqa: E402

DECODERS = ["none", "union_find", "mwpm", "correlated", "belief_matching", "bposd"]
FIG = Path(__file__).resolve().parent.parent / "results" / "figures"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--p", type=float, default=0.01, help="SI1000 prior strength")
    args = ap.parse_args()

    table = []
    for path in args.runs:
        run = load_run(Path(path))
        rounds, basis = [int(x) for x in run["rounds"]], run["basis"]
        priors = {"si1000": NoiseModel.si1000(args.p)}
        cal_file = run.get("calibration_file")
        if cal_file and Path(cal_file).exists():
            cal = json.loads(Path(cal_file).read_text())["calibration"]
            priors["calibrated"] = NoiseModel.from_calibration(
                cal, {int(k): v for k, v in run["physical"].items()})
        data = {}
        for r in rounds:
            clean = build_memory_circuit(3, r, basis, True, None)
            data[r] = measurements_to_detections(run[f"raw_r{r}"], clean)
        name = Path(path).stem
        print(f"\n=== {name}")
        for pname, prior in priors.items():
            for dec in DECODERS:
                if dec == "none" and pname != "si1000":
                    continue
                k, n, t_tot = [], [], 0.0
                for r in rounds:
                    det, obs = data[r]
                    t0 = time.perf_counter()
                    res = decode(det, obs, build_memory_circuit(3, r, basis, True, prior), dec)
                    t_tot += time.perf_counter() - t0
                    k.append(res.errors)
                    n.append(res.shots)
                fit = fit_ler(rounds, k, n)
                us = 1e6 * t_tot / sum(n)
                row = dict(run=name, prior=pname if dec != "none" else "-", decoder=dec,
                           eps=fit.eps, eps_err=fit.eps_err,
                           p_l=(np.array(k) / np.array(n)).tolist(), us_per_shot=us)
                table.append(row)
                print(f"  {row['prior']:>10} {dec:>16}: eps_L = {100*fit.eps:6.2f} +/- "
                      f"{100*fit.eps_err:.2f} %   {us:8.1f} us/shot")

    FIG.mkdir(parents=True, exist_ok=True)
    (FIG / "decoder_comparison.json").write_text(json.dumps(table, indent=1))

    runs = list(dict.fromkeys(r["run"] for r in table))
    fig, axes = plt.subplots(1, len(runs), figsize=(4.2 * len(runs), 3.8), squeeze=False)
    for ax, run in zip(axes[0], runs):
        rows = [r for r in table if r["run"] == run]
        labels = [f"{r['decoder']}\n{r['prior']}" for r in rows]
        ax.barh(range(len(rows)), [100 * r["eps"] for r in rows],
                xerr=[100 * r["eps_err"] for r in rows],
                color=["grey" if r["decoder"] == "none" else
                       ("C0" if r["prior"] == "si1000" else "C1") for r in rows])
        ax.set_yticks(range(len(rows)), labels, fontsize=6)
        ax.invert_yaxis()
        ax.set_xlabel("$\\varepsilon_L$ per round (%)")
        ax.set_title(run.replace("20261121-", ""), fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "decoder_comparison.png", dpi=150)
    print("\nsaved results/figures/decoder_comparison.png and .json")


if __name__ == "__main__":
    main()
