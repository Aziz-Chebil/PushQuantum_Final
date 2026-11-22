# Surface-code memory on IQM Emerald

A distance-3 rotated surface code (17 qubits) run as a quantum memory on **IQM Emerald
hardware**, decoded, and fitted to a logical error rate per round (ε_L). PushQuantum 2026 IQM
challenge, Qiskit and Qrisp tracks ([CHALLENGE.md](CHALLENGE.md)).

Constraints of every run: no mid-circuit reset, zero SWAPs, native CZ + PRX gates, IQM
server-side dynamical decoupling (DD).

## Result

| configuration (same 17 qubits, calibration set `ecec1d0a`) | ε_L per round |
|---|---|
| no error correction (raw logical readout) | 8.59 ± 0.15 % |
| starting circuit | 2.70 ± 0.10 % |
| **optimized circuit** | **2.12 ± 0.03 %** |

Z memory, rounds 1–5, BP+OSD decoder with SI1000 prior (p = 1 %). The optimized result pools
4 hardware campaigns (16 000 shots per round count). All numbers are refit with one common
fitter, `analysis/refit_all.py`; full table and caveats in [RESULTS.md](RESULTS.md).

Error correction lowers the error per round about 4×; the optimized circuit lowers it about
20 % further on the same qubits. The circuit comparison spans separate jobs, so drift is not
excluded; a single job containing both circuits also favours the optimized one
(5.90 % vs 8.13 % failures at 3 rounds, ~3.9σ).

## Challenge requirements

### Required

| # | requirement | where | what we did |
|---|---|---|---|
| 1 | circuit generator: distance, rounds, Z/X memory, `no_reset`, optional noise model | `pipeline/surface_code.py` → `qec/circuit.py`; `qrisp_track/surface_code.py` | `make_stim_circuit(distance, rounds, basis, no_reset, noise)` returns a `stim.Circuit`, `make_qiskit_circuit` the same as Qiskit; `make_qrisp_experiment(...)` returns a Qrisp `DetectorExperiment`. Every hardware run uses `no_reset=True`: detectors are rewritten to m_r ⊕ m_(r−2). Z and X memory both run on Emerald. |
| 2 | decode and fit ε_L with uncertainty; DEM matches the executed circuit | `pipeline/qec/decode.py`, `qec/analysis.py`, `analysis/refit_all.py` | DEM built from the same circuit that is sent to IQM (checked by payload verification, below). Each round count n = 1–5 is a separate circuit; P_L(n) fitted with `2 P_L = 1 − (1 − 2ε_L)^n` by binomial maximum likelihood, 1σ errors from the Fisher information. |

### Optional

| # | item | where | result |
|---|---|---|---|
| 3 | simulation baseline of the no-reset circuit | `pipeline/experiments/simulate.py`, `compare_sim_hw.py`; `circuit_optimization/simulate_qec.py` | No-reset is not worse than reset in simulation (0.345 vs 0.353 % at p = 0.3 %). Simulating the exact circuit on the exact qubits with calibration data predicts the Emerald hardware within ~15 % (3.84 % sim vs 3.15–3.41 % hardware, MWPM). |
| 4 | custom qubit placement, zero SWAPs | `pipeline/experiments/layout.py`; `qrisp_track/run_on_hardware.py` (vf2pp) | All SWAP-free embeddings scored with live calibration; ancillas measured to be bad are excluded; every CZ verified to be a coupler. 0 SWAPs in every run. |

### Extensions

| extension | where | result |
|---|---|---|
| dynamical decoupling | `pipeline/experiments/run_hardware.py --dd` | ε_L 3.24 → 2.82 % (Emerald), 13.07 → 11.62 % (Garnet). Adopted. |
| pulse level (PulLA) | `pipeline/experiments/run_pulla.py` | Works end to end; standard and denser DD at pulse level give 2.79 / 2.56 % vs 2.61 % at circuit level (all within ~1σ). No significant gain. |
| other decoders | `pipeline/experiments/decoder_comparison.py` | MWPM, correlated MWPM, Union-Find, belief-matching, BP+OSD on identical shots. Ranking on every Emerald run: BP+OSD < belief-matching < correlated MWPM < MWPM < Union-Find. |
| decoder weights from calibration data | `pipeline/experiments/characterized_decoding.py` | Calibration-based prior is worse than generic SI1000 (+0.19 to +0.29 pp). |
| extended noise characterization | `pipeline/experiments/characterize_measurement.py`, `qec/learn.py` | Repeated-measurement characterization identifies the failing ancillas and improves the decoder prior (−0.11 to −0.16 pp); error rates learned from QEC data help on Garnet (−1.1 to −1.3 pp), not on Emerald. |
| different 5×5 regions | `pipeline/experiments/regions.py` | 4 regions of Emerald: ε_L 2.65 / 3.00 / 3.20 / 3.87 % (BP+OSD). |
| distance 5 | `circuit_optimization/simulate_qec.py` | Not possible on hardware: Emerald has no SWAP-free d = 5 embedding. Simulated only (d = 3/5/7). |
| Garnet | `pipeline/` | One SWAP-free d = 3 patch; near threshold: 10.57 ± 0.18 % vs 12.16 % without correction. |
| NVIDIA Ising predecoder | | Not attempted. |

## Pipeline

```
Stim circuit → Qiskit / Qrisp → IQM Emerald (Resonance) → detection events → decoder → P_L(rounds) → fit → ε_L
```

| stage | choice | evidence |
|---|---|---|
| qubit placement | zero-SWAP patch scored on calibration data, excluding ancillas we measured to be bad (QB33, QB26) | ε_L 7.9 → 3.4 % (MWPM; separate runs) |
| circuit | CZ layer order reversed every second round; redundant Hadamards cancelled (94 instead of 104 rotations at 3 rounds, the minimum for this construction) | 2.70 → 2.12 % |
| scheduling | barriers between layers, so IQM's compiler cannot interleave them | failures at 3 rounds 28.7 % → 10.1 % |
| idle noise | IQM server-side DD | ε_L 3.24 → 2.82 % |
| decoder | BP+OSD; correlated MWPM as the ~40× faster option (within ~0.1 pp) | lowest error of 5 decoders on every Emerald run, identical shots |
| decoder prior | SI1000 for gates and idling + readout errors from our repeated-measurement characterization | −0.11 to −0.16 pp; better in 11/13 runs (BP+OSD), 12–13/13 (MWPM); each gain within error bars |
| verification | IQM payloads parsed back into Stim; every single X/Y/Z fault checked for identical detector response | 76/76 circuits pass; one payload downloaded from IQM's server matches the rebuild |

Tested and **not adopted**:

* Pauli twirling: no consistent effect (−0.09 pp with Qrisp, +0.05 pp with Qiskit; 20 000
  shots per arm, both within 1–2σ).
* Approximate likelihood decoder: matches BP+OSD on hardware, no gain.
* Error rates learned from the QEC data: −1.1 to −1.3 pp on IQM Garnet, no gain on Emerald.

## Limitations

* d = 3 only. Emerald has no zero-SWAP d = 5 patch, so distance scaling was not measured.
* The optimized circuit was decoded with SI1000 only, not yet with the characterized readout
  prior, and has not yet gone through the per-fault payload check.
* The Qrisp-track runs (1.79–1.92 %) are not comparable with the table above: their qubits
  were chosen automatically and not recorded.

## Repository

| folder | author | content | write-up |
|---|---|---|---|
| [`pipeline/`](pipeline) | Aziz Chebil | Qiskit-track hardware pipeline: placement, Emerald/Garnet runs, readout characterization, 5 decoders, decoder priors, payload verification | [pipeline/RESULTS.md](pipeline/RESULTS.md) |
| [`circuit_optimization/`](circuit_optimization) | Ayoub Lefi | optimized circuit, compiler/barrier audit, simulation study (d = 3/5/7), likelihood decoder, twirling study on hardware | [OPTIMIZATION.md](circuit_optimization/OPTIMIZATION.md) |
| [`qrisp_track/`](qrisp_track) | Mohamed Amine Abidi | Qrisp-track pipeline (Stim → Qrisp `DetectorExperiment` → Emerald), zero-SWAP pass manager, X basis, twirling prototype | [qrisp_track/RESULTS.md](qrisp_track/RESULTS.md) |
| [`analysis/`](analysis) | | `refit_all.py`: one fitter over the hardware data of all three folders | [RESULTS.md](RESULTS.md) |

Each folder contains its code and raw hardware data and has its own `requirements.txt`. Use
one virtual environment per folder; pinned versions differ.

## Quick start

No hardware token is needed for the commands below (Windows paths; use `.venv/bin` on
Linux/macOS).

```bash
# pipeline: tests, decoder study, common refit of all runs
cd pipeline
python -m venv .venv && .venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python -m pytest -q
.venv/Scripts/python experiments/characterized_decoding.py
.venv/Scripts/python ../analysis/refit_all.py      # reproduces the table in RESULTS.md
cd ..

# optimized circuit: correctness tests
cd circuit_optimization
python -m venv .venv && .venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python -m unittest -v verify_qec verify_optimization
cd ..

# Qrisp track: local simulation of the same code path
cd qrisp_track
python -m venv .venv && .venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python benchmark_qrisp.py
```

Hardware runs need an IQM Resonance token in `IQM_TOKEN` (environment or a git-ignored
`.env`). No token is stored in this repository. Hardware commands are in each folder's
README.

## Team

Aziz Chebil, Ayoub Lefi, Mohamed Amine Abidi, Ahmed Amine Loukil, Wassim zormati.
Built on the IQM PushQuantum challenge starter kit. Referenced papers (Google Quantum AI,
Nature 2021, 2023, 2024) are cited in
[circuit_optimization/PAPER_STUDY.md](circuit_optimization/PAPER_STUDY.md).
