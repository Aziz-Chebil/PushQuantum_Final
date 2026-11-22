"""Decode saved hardware runs, fit the logical error rate, make figures.

For each run file it compares decoder priors (the DEM's noise model):
  uniform p, SI1000 p, and the calibration of the exact qubits used,
and decoders: MWPM, correlated MWPM, and no decoding (raw parity).

    .venv/Scripts/python experiments/analyze.py results/hw/<run>.npz [more.npz ...]
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
from qec.analysis import eps_from_single, fit_ler, plot_ler  # noqa: E402
from qec.decode import decode, measurements_to_detections, syndrome_diagnostics  # noqa: E402
from qec.hardware import load_run  # noqa: E402
from qec.learn import cross_validated_decode, graph_edges, summarize_learned  # noqa: E402


def priors(run: dict, p: float, cal_file: str | None = None) -> dict[str, NoiseModel]:
    out = {"uniform": NoiseModel.uniform(p), "si1000": NoiseModel.si1000(p)}
    cal_file = cal_file or run.get("calibration_file")
    if cal_file and Path(cal_file).exists():
        cal = json.loads(Path(cal_file).read_text())["calibration"]
        phys = {int(k): v for k, v in run["physical"].items()}
        out["calibrated"] = NoiseModel.from_calibration(cal, phys)
    return out


def analyze(path: Path, p: float, decoders=("mwpm", "correlated", "none"),
            cal_file: str | None = None):
    run = load_run(path)
    rounds = [int(r) for r in run["rounds"]]
    basis = run["basis"]
    code = RotatedSurfaceCode(3)
    models = priors(run, p, cal_file)
    rows, diag, learned_summary = [], {}, {}
    for r in rounds:
        clean = build_memory_circuit(3, r, basis, True, None)
        det, obs = measurements_to_detections(run[f"raw_r{r}"], clean)
        diag[r] = syndrome_diagnostics(det)
        for mname, nm in models.items():
            noisy = build_memory_circuit(3, r, basis, True, nm)
            for dec in decoders:
                if dec == "none" and mname != "uniform":
                    continue
                res = decode(det, obs, noisy, dec)
                rows.append(dict(rounds=r, prior=mname if dec != "none" else "-",
                                 decoder=dec, errors=res.errors, shots=res.shots,
                                 p_l=res.p_logical, p_l_err=res.stderr,
                                 eps_single=eps_from_single(res.p_logical, r)))
        # Data-driven weights (structure from the uniform prior, probabilities
        # learned from the detection events, 2-fold cross-validated).
        if r >= 2:
            res, learned = cross_validated_decode(
                det, obs, build_memory_circuit(3, r, basis, True, models["uniform"]))
            rows.append(dict(rounds=r, prior="data", decoder="learned", errors=res.errors,
                             shots=res.shots, p_l=res.p_logical, p_l_err=res.stderr,
                             eps_single=eps_from_single(res.p_logical, r)))
            learned_summary[r] = summarize_learned(
                build_memory_circuit(3, r, basis, True, models["uniform"]),
                learned, code, run["physical"])

    fits = {}
    for key in sorted({(x["prior"], x["decoder"]) for x in rows}):
        sel = [x for x in rows if (x["prior"], x["decoder"]) == key]
        k = [x["errors"] for x in sel]
        n = [x["shots"] for x in sel]
        rr = [x["rounds"] for x in sel]
        f1 = fit_ler(rr, k, n)
        fa = fit_ler(rr, k, n, free_amplitude=True) if len(rr) >= 3 else None
        fits[key] = (f1, fa, rr, k, n)

    tag = path.stem
    print(f"\n=== {tag}  ({basis}-memory, DD={'on' if run.get('dd') else 'off'}, "
          f"qubits {sorted(run['physical'].values())}) ===")
    print(" rounds  det.fraction  error-free shots")
    for r in rounds:
        print(f"  {r:>4}     {diag[r]['detection_fraction']:.3f}         {diag[r]['error_free_shots']:.3f}")
    print(f" {'prior':>10} {'decoder':>10}   eps_L (A=1)               eps_L (free A)")
    for (prior, dec), (f1, fa, *_rest) in fits.items():
        print(f" {prior:>10} {dec:>10}   {100*f1.eps:6.2f} +/- {100*f1.eps_err:.2f} %"
              f"      {'' if fa is None else f'{100*fa.eps:6.2f} +/- {100*fa.eps_err:.2f} %  A={fa.amp:.3f}'}")

    if learned_summary:
        rl = max(learned_summary)
        cal_prior = {}
        if "calibrated" in models:
            # What the calibration-based model *expected* on the same edges.
            cal_circ = build_memory_circuit(3, rl, basis, True, models["calibrated"])
            cal_edges = graph_edges(cal_circ)
            cal_prior = {x["stabilizer"]: x["timelike_learned"] for x in summarize_learned(
                cal_circ, {"edges": cal_edges, "p": {e.dets: e.p_prior for e in cal_edges}},
                code, run["physical"])}
        print(f" time-like error per stabilizer ({rl} rounds; ancilla readout + idling):")
        print(f"   {'stabilizer':>10}  learned  calibration-model")
        for x in sorted(learned_summary[rl], key=lambda x: -x["timelike_learned"]):
            x["timelike_calibrated"] = cal_prior.get(x["stabilizer"])
            c = f"{100*cal_prior[x['stabilizer']]:5.1f} %" if x["stabilizer"] in cal_prior else "-"
            print(f"   {x['stabilizer']:>10}  {100*x['timelike_learned']:5.1f} %   {c}")

    # ---- figures
    outdir = path.parent.parent / "figures"
    outdir.mkdir(exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for i, ((prior, dec), (f1, fa, rr, k, n)) in enumerate(fits.items()):
        p_obs = np.array(k) / np.array(n)
        plot_ler(ax[0], rr, p_obs, np.sqrt(p_obs * (1 - p_obs) / np.array(n)), f1,
                 label=f"{dec} / {prior}", color=f"C{i}")
    ax[0].set_title(f"{tag}")
    ax[0].legend(fontsize=7)
    # detector firing rate per round, per stabilizer
    rmax = max(rounds)
    clean = build_memory_circuit(3, rmax, basis, True, None)
    det, _ = measurements_to_detections(run[f"raw_r{rmax}"], clean)
    coords = clean.get_detector_coordinates()
    stabs = sorted({(int(c[0]), int(c[1])) for c in coords.values()})
    grid = np.full((len(stabs), rmax + 1), np.nan)
    for j, c in coords.items():
        grid[stabs.index((int(c[0]), int(c[1]))), int(c[2])] = det[:, j].mean()
    im = ax[1].imshow(grid, aspect="auto", cmap="viridis", vmin=0, vmax=0.5)
    anc_name = {code.coords[s.ancilla]: f"{s.basis}@{run['physical'][str(s.ancilla)]}"
                for s in code.stabilizers}
    ax[1].set_yticks(range(len(stabs)), [anc_name.get(s, str(s)) for s in stabs], fontsize=7)
    ax[1].set_xlabel("round (last column = final data readout)")
    ax[1].set_title(f"detector firing rate, {rmax} rounds")
    fig.colorbar(im, ax=ax[1])
    fig.tight_layout()
    fig.savefig(outdir / f"{tag}.png", dpi=150)

    summary = {
        "run": tag, "basis": basis, "dd": run.get("dd"), "physical": run["physical"],
        "diagnostics": {r: {k: v for k, v in d.items() if k != "detector_firing_rate"}
                        for r, d in diag.items()},
        "points": rows,
        "learned_timelike": learned_summary,
        "fits": {f"{p}/{d}": {"eps": f1.eps, "eps_err": f1.eps_err,
                               "eps_freeA": fa.eps if fa else None,
                               "eps_freeA_err": fa.eps_err if fa else None,
                               "A": fa.amp if fa else None}
                 for (p, d), (f1, fa, *_r) in fits.items()},
    }
    (outdir / f"{tag}.json").write_text(json.dumps(summary, indent=1, default=float))
    print(f"saved figures/{tag}.png and .json")
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--p", type=float, default=0.01, help="strength for the uniform/SI1000 priors")
    ap.add_argument("--calibration", help="calibration json to use for the calibrated prior "
                    "(default: the one recorded with the run)")
    args = ap.parse_args()
    for r in args.runs:
        analyze(Path(r), args.p, cal_file=args.calibration)


if __name__ == "__main__":
    main()
