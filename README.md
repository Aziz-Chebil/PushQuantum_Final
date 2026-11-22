# Surface-code memory on IQM Emerald: team final

Distance-3 rotated surface code (17 qubits) running as a quantum memory on **real IQM Emerald
hardware**: no mid-circuit reset, zero SWAPs, native CZ + PRX, server-side dynamical
decoupling (DD). Repeated stabilizer rounds are decoded and the logical error per round ε_L is
fitted. Built for the PushQuantum 2026 IQM challenge ([CHALLENGE.md](CHALLENGE.md)), Qiskit and
Qrisp tracks.

## Result

> **ε_L = 2.12 ± 0.03 % per round** (Z memory, BP+OSD, SI1000 prior; 16 000 shots per round
> count, rounds 1–5, 4 independent hardware campaigns), using the optimized circuit
> (`circuit_optimization/`) on the qubit patch chosen by the hardware pipeline (`pipeline/`).
> Without error correction the same patch loses its logical bit at 8.6 % per round.

The same patch with the starting circuit gives 2.65–2.70 %; the optimized circuit lowers
the logical error about 20 % on the same qubits and calibration set. All numbers in
[RESULTS.md](RESULTS.md) are refit with one common fitter (`analysis/refit_all.py`), because
the three projects originally used different fit models.

## The approach, stage by stage (best of the three projects)

| stage | what we use | from | evidence |
|---|---|---|---|
| qubit placement | calibration-scored, zero-SWAP patch that excludes ancillas measured to be bad (QB33, QB26) | `pipeline/` | re-placing to avoid them: ε_L 7.9 → 3.4 % (MWPM; partly drift possible) |
| circuit | CZ layer order reversed every second round, redundant Hadamards cancelled (94 instead of 104 rotations at 3 rounds, proven minimal for this construction), layer barriers kept | `circuit_optimization/` | 2.12 % vs 2.65–2.70 % on the same 17 qubits and calibration set |
| scheduling | barriers between layers so IQM's compiler cannot interleave layers | `circuit_optimization/` | without them: 28.7 % → 10.1 % failures at 3 rounds on the same patch |
| idle noise | IQM server-side dynamical decoupling | `pipeline/` | idle dephasing of a bad ancilla 18.3 % → 2.7 %; ε_L 3.24 → 2.82 % |
| decoder | BP+OSD (correlated MWPM as the fast option) | all | best of 5 decoders on identical shots, every Emerald run |
| decoder prior | SI1000 gates/idle + readout from our own repeated-measurement characterization | `pipeline/` (E5) | better than SI1000 in 12–13 of 13 runs (small, consistent) |
| learned error rates | only on noisy patches | `pipeline/` (E7) | −1.1 to −1.3 pp on Garnet, no gain on Emerald |
| verification | IQM payload parsed back to Stim; every single Pauli fault gives identical detector response | `pipeline/` (E6) | 76/76 circuits, payload identical to the one on IQM's server |
| Qrisp track | Qrisp `DetectorExperiment` pipeline with zero-SWAP native layout | `qrisp_track/` | 1.8–1.9 % per round on its own patch (see caveats) |

Tested and **not adopted**: Pauli twirling (no consistent effect: −0.09 pp with Qrisp,
+0.05 pp with Qiskit, 20 000 shots per arm), the approximate likelihood decoder (no gain over
BP+OSD on hardware), data-learned rates on Emerald.

## Repository layout

| folder | author | content |
|---|---|---|
| [`pipeline/`](pipeline) | Aziz Chebil | Qiskit-track hardware pipeline (`qec/` package): layout from calibration and measured errors, Emerald/Garnet runs, readout characterization, 5 decoders, characterized and learned decoder priors, circuit verification against IQM payloads. Write-up: [pipeline/RESULTS.md](pipeline/RESULTS.md) |
| [`circuit_optimization/`](circuit_optimization) | Ayoub Lefi | Optimized no-reset circuit (alternating CZ order, Hadamard cancellation, minimal rotation count), barrier/translation audit, simulation study (668 conditions, d = 3/5/7), likelihood decoder, matched comparison and large twirling study on hardware. Write-up: [OPTIMIZATION.md](circuit_optimization/OPTIMIZATION.md), [resume.txt](circuit_optimization/resume.txt) |
| [`qrisp_track/`](qrisp_track) | Mohamed Amine Abidi | Qrisp-track pipeline (Stim → Qrisp `DetectorExperiment` → Emerald), zero-SWAP native pass manager working around the IQM router, X-basis support, SI1000 prior, twirling prototype. Write-up: [qrisp_track/RESULTS.md](qrisp_track/RESULTS.md) |
| [`analysis/`](analysis) | | `refit_all.py`: one fitter over all three projects' hardware data |

Each folder is the project as its author ran it (code + raw hardware data), with build
leftovers removed. Each has its own `requirements.txt`; use one virtual environment per folder
(their pinned versions differ, e.g. `ldpc 2.4.1` in `qrisp_track` needs `numpy <= 2.2.6`).

## Reproduce

```bash
# hardware pipeline: tests and the main analyses (no hardware token needed)
cd pipeline && python -m venv .venv && .venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python -m pytest -q
.venv/Scripts/python experiments/characterized_decoding.py
cd .. && pipeline/.venv/Scripts/python analysis/refit_all.py      # the table in RESULTS.md

# optimized circuit: tests and simulation study
cd circuit_optimization && python -m unittest -v verify_qec verify_optimization

# Qrisp track: free local simulation of the same code path
cd qrisp_track && python benchmark_qrisp.py
```

Hardware runs need an IQM Resonance token in the environment variable `IQM_TOKEN` (or `.env`,
which is git-ignored). No token is stored in this repository.

`circuit_optimization/compare_peer.py` and `twirl_memory.py` import the hardware pipeline;
pass `--peer ../pipeline`.

## Team

Aziz Chebil, Ayoub Lefi, Mohamed Amine Abidi, Mohamed Amine Loukil, Wassim Zaourmati.
Built on the IQM PushQuantum challenge starter kit. Papers referenced (Google Quantum AI,
Nature 2021, 2023, 2024) are cited by DOI in
[circuit_optimization/PAPER_STUDY.md](circuit_optimization/PAPER_STUDY.md).
