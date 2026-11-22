"""Simulation baseline vs hardware: did we predict the hardware result?

For every hardware run, the exact circuit that ran is simulated in Stim under
two noise models built for the exact physical qubits used:

  1. "calibration"   - Resonance calibration only (what we could predict
                       before spending QPU time)
  2. "calibration + measured ancilla errors" - the same, but each ancilla's
                       readout error is raised to the time-like error rate we
                       learned from that run's detection events (qec/learn.py)

Hardware and simulations are decoded with the same decoder (MWPM, calibrated
prior), so differences come from the physics, not the decoder. If model 2
closes most of the gap, the few bad ancillas we identified are the dominant
reason the hardware is worse than its calibration.

    .venv/Scripts/python experiments/compare_sim_hw.py results/hw/<run>.npz ...
"""

import argparse
import json
import sys
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import NoiseModel, RotatedSurfaceCode, build_memory_circuit  # noqa: E402
from qec.analysis import fit_ler  # noqa: E402
from qec.decode import decode, measurements_to_detections  # noqa: E402
from qec.hardware import load_run  # noqa: E402
from qec.learn import cross_validated_decode, summarize_learned  # noqa: E402

FIG = Path(__file__).resolve().parent.parent / "results" / "figures"


def learned_ancilla_errors(run, code, uniform) -> dict[int, float]:
    """Time-like error per ancilla, learned from the run's longest sweep point."""
    r = max(int(x) for x in run["rounds"])
    clean = build_memory_circuit(3, r, run["basis"], True, None)
    det, obs = measurements_to_detections(run[f"raw_r{r}"], clean)
    prior = build_memory_circuit(3, r, run["basis"], True, uniform)
    _, learned = cross_validated_decode(det, obs, prior)
    rows = summarize_learned(prior, learned, code, run["physical"])
    by_name = {f"{s.basis}@{run['physical'][str(s.ancilla)]}": s.ancilla for s in code.stabilizers}
    return {by_name[x["stabilizer"]]: x["timelike_learned"] for x in rows}


def sweep(model, rounds, basis, prior, shots, seed=0):
    out = []
    for r in rounds:
        c = build_memory_circuit(3, r, basis, True, model)
        det, obs = c.compile_detector_sampler(seed=seed + r).sample(shots, separate_observables=True)
        out.append(decode(det, obs, build_memory_circuit(3, r, basis, True, prior)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--shots", type=int, default=50_000)
    args = ap.parse_args()
    code = RotatedSurfaceCode(3)
    uniform = NoiseModel.uniform(0.01)

    fig, axes = plt.subplots(1, len(args.runs), figsize=(4.6 * len(args.runs), 4), squeeze=False)
    table = []
    for ax, path in zip(axes[0], args.runs):
        run = load_run(Path(path))
        rounds = [int(x) for x in run["rounds"]]
        basis = run["basis"]
        cal = json.loads(Path(run["calibration_file"]).read_text())["calibration"]
        phys = {int(k): v for k, v in run["physical"].items()}
        m_cal = NoiseModel.from_calibration(cal, phys)

        m_meas = deepcopy(m_cal)
        for anc, p in learned_ancilla_errors(run, code, uniform).items():
            m_meas.p_meas_q[anc] = max(m_cal.meas(anc), p)

        hw = []
        for r in rounds:
            clean = build_memory_circuit(3, r, basis, True, None)
            det, obs = measurements_to_detections(run[f"raw_r{r}"], clean)
            hw.append(decode(det, obs, build_memory_circuit(3, r, basis, True, m_cal)))

        curves = {"hardware": hw,
                  "sim: calibration": sweep(m_cal, rounds, basis, m_cal, args.shots),
                  "sim: calibration + measured ancilla errors":
                      sweep(m_meas, rounds, basis, m_cal, args.shots)}
        name = Path(path).stem
        for (label, res), style in zip(curves.items(), ["o", "--", ":"]):
            p = np.array([x.p_logical for x in res])
            e = np.array([x.stderr for x in res])
            fit = fit_ler(rounds, [x.errors for x in res], [x.shots for x in res])
            table.append({"run": name, "model": label, "p_l": p.tolist(),
                          "eps": fit.eps, "eps_err": fit.eps_err})
            lab = f"{label}  ($\\varepsilon_L$={100*fit.eps:.1f}%)"
            if style == "o":
                ax.errorbar(rounds, p, yerr=e, fmt="o", color="k", capsize=3, label=lab)
            else:
                ax.plot(rounds, p, style, marker=".", label=lab)
            print(f"{name:>34}  {label:<45} eps_L = {100*fit.eps:5.2f} +/- {100*fit.eps_err:.2f} %"
                  f"   P_L = {' '.join(f'{x:.3f}' for x in p)}")
        ax.set_title(name.replace("20261121-", ""), fontsize=9)
        ax.set_xlabel("rounds")
        ax.set_ylabel("$P_L$ (MWPM, calibrated prior)")
        ax.set_ylim(0, 0.5)
        ax.legend(fontsize=6.5)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "sim_vs_hardware.png", dpi=150)
    (FIG / "sim_vs_hardware.json").write_text(json.dumps(table, indent=1))
    print("saved results/figures/sim_vs_hardware.png and .json")


if __name__ == "__main__":
    main()
