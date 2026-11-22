"""Decoder DEM combining calibration-free models, our characterization
experiment and error rates learned from the QEC data itself.

For every hardware run the same shots are decoded with four DEMs:

  si1000          generic model, p = 1 %
  si1000+char     SI1000 gates/idling + readout from the repeated-measurement
                  characterization (experiments/characterized_decoding.py)
  char+learned    si1000+char structure; every graph edge's probability
                  learned from the detection events (p_ij method), except the
                  time-like readout edges, which keep the characterized values
  all learned     as char+learned, but the readout edges are learned too

Learned rates are pooled over every occurrence of the same physical edge
(all bulk rounds, all round counts of the run; qec.learn.pooled_edge_probabilities):
2 000 training shots per circuit are too few to estimate each edge on its own
(simulated: no gain over the prior at that size, a gain once pooled).

Learned DEMs are 2-fold cross-validated: learned on one half of the shots of
every circuit of the run, scored on the other halves, then swapped. No DEM is
scored on the shots it was fit to.

    .venv/Scripts/python experiments/combined_dem.py
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from characterized_decoding import HW, RUNS, priors  # noqa: E402

from qec import build_memory_circuit  # noqa: E402
from qec.analysis import fit_ler  # noqa: E402
from qec.decode import decode_dem, measurements_to_detections  # noqa: E402
from qec.hardware import load_run  # noqa: E402
from qec.learn import graph_edges, pooled_edge_probabilities, reweighted_dem, timelike_edges  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RUNS = RUNS + [("Emerald Z, DD on (21:08)", "20261121-210819_emerald_Z_ddon.npz")]
DECODERS = ["mwpm", "correlated", "bposd"]
DEMS = ["si1000", "si1000+char", "char+learned", "all learned"]


def errors_per_round(data, circuits, dem_name, decoder, folds=2, seed=0):
    """Logical failures for each round count (list, in the order of `circuits`)."""
    rounds = list(circuits)
    if dem_name in ("si1000", "si1000+char"):
        return [decode_dem(*data[r], circuits[r].detector_error_model(decompose_errors=True),
                           decoder).errors for r in rounds]
    split = {r: np.random.default_rng(seed + r).permutation(len(data[r][0])) for r in rounds}
    errs = {r: 0 for r in rounds}
    for f in range(folds):
        test = {r: split[r][f::folds] for r in rounds}
        train = {r: np.setdiff1d(split[r], test[r]) for r in rounds}
        pooled = pooled_edge_probabilities([(circuits[r], data[r][0][train[r]]) for r in rounds])
        for r, learned in zip(rounds, pooled):
            c = circuits[r]
            keep = timelike_edges(c, graph_edges(c)) if dem_name == "char+learned" else set()
            dem = reweighted_dem(c, keep=keep, learned=learned)
            det, obs = data[r]
            errs[r] += decode_dem(det[test[r]], obs[test[r]], dem, decoder).errors
    return [errs[r] for r in rounds]


def main():
    table = []
    print(f"{'run':<26}" + "".join(f"{d[:12] + '/' + x[:4]:>19}" for d in DEMS for x in DECODERS))
    for label, fname in RUNS:
        run = load_run(HW / fname)
        rounds, basis = [int(r) for r in run["rounds"]], run["basis"]
        pri = priors(run)
        model = {"si1000": pri["si1000"], "si1000+char": pri["si1000+char"],
                 "char+learned": pri["si1000+char"], "all learned": pri["si1000+char"]}
        data = {r: measurements_to_detections(run[f"raw_r{r}"], build_memory_circuit(3, r, basis, True, None))
                for r in rounds}
        row, cells = {"run": label, "file": fname}, []
        for d in DEMS:
            circuits = {r: build_memory_circuit(3, r, basis, True, model[d]) for r in rounds}
            for dec in DECODERS:
                k = errors_per_round(data, circuits, d, dec)
                n = [len(data[r][0]) for r in rounds]
                fit = fit_ler(rounds, k, n)
                row[f"{d}/{dec}"] = (fit.eps, fit.eps_err)
                cells.append(f"{100 * fit.eps:5.2f}±{100 * fit.eps_err:.2f}")
        table.append(row)
        print(f"{label:<26}" + "".join(f"{c:>19}" for c in cells), flush=True)

    print()
    for dec in DECODERS:
        for d in DEMS[1:]:
            for base in ("si1000", "si1000+char"):
                if d == base:
                    continue
                diff = [100 * (r[f"{d}/{dec}"][0] - r[f"{base}/{dec}"][0]) for r in table]
                wins = sum(x < 0 for x in diff)
                print(f"{dec:>10} {d:>13} - {base:<12} = {np.mean(diff):+.2f} pp (mean), "
                      f"better in {wins}/{len(table)}")
    out = ROOT / "results" / "figures" / "combined_dem.json"
    out.write_text(json.dumps(table, indent=1))
    print(f"saved {out}")


if __name__ == "__main__":
    main()
