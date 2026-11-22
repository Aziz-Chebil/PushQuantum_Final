# QEC-Emerald — rotated surface-code memory on IQM Emerald

Our solution to the IQM *Surface Code on IQM Emerald* challenge (Qiskit track).

```
Stim (native H+CZ, no-reset detectors) → Qiskit → zero-SWAP pinned transpile → IQM Emerald
      → raw shots → detection events (m2d) → MWPM / correlated MWPM → P_L(n) → binomial-MLE fit → ε_L
```

## What is ours and why it matters

| Piece | What we did | Judging criterion |
|---|---|---|
| `qec/code.py` | Own rotated-surface-code geometry, stabilizers, logicals and CX schedule (hook errors oriented so the **circuit distance equals d**, checked by Stim). | Correctness |
| `qec/circuit.py` | Own circuit generator written **directly in Emerald's native gate set (H/PRX + CZ)**. A Hadamard-frame tracker cancels redundant Hadamards, so the simulated/decoded circuit is the one the chip runs. **Proper multi-round no-reset detectors** (`m_r ⊕ m_{r-2}`, final `data ⊕ m_R ⊕ m_{R-1}`), Z and X memory. | Correctness, sophistication |
| `qec/noise.py` | Per-qubit / per-coupler noise model built from **live Resonance calibration** (1Q, CZ, readout fidelities → Pauli channels; T1/T2 → idle channels during gates and the long measurement windows). Also uniform and SI1000 baselines. | Noise-model improvements |
| `qec/layout.py` | **VF2 search over every SWAP-free embedding** of the code into the coupler graph, each scored by a Stim simulation under its own calibrated noise model (role-aware: which qubit is data / X / Z ancilla). | Low SWAPs, sophistication |
| `qec/hardware.py` | Transpile with the layout pinned and `routing_method="none"` — the transpiler *cannot* insert a SWAP, and `swap_report` independently verifies every CZ sits on a physical coupler. One batched job per sweep (all round counts share one calibration), optional **IQM server-side dynamical decoupling**, raw per-shot data saved. | Low SWAPs, advanced features |
| `qec/decode.py`, `qec/analysis.py` | MWPM and correlated MWPM (PyMatching) with the DEM built from the calibrated model; **binomial maximum-likelihood fit** of `P_L(n) = ½(1 − A(1−2ε_L)^n)` with Fisher-information errors, README form (A = 1) and free-A (SPAM-separated) form. | LER statistics |
| `tests/` | 24 tests: determinism of every detector (proves the no-reset algebra), circuit distance = d, reset/no-reset give identical detection events for the same physical error, Qiskit translation reproduces Stim bit-for-bit under injected errors (Aer), layout nativeness, fit recovers a known ε. | Correctness |

## Challenge entry points

`surface_code.py` exposes the scaffold's names on top of `qec/`:
`make_stim_circuit(distance, rounds, basis, no_reset, noise)`, `make_qiskit_circuit`,
`DEFAULT_NOISE`, `build_emerald_qubit_map`, `print_qubit_map`, `run_qiskit`,
`extract_syndromes`, `decode_hardware_results(syndromes) -> (P_L, err)` and
`fit_logical_error_rate(rounds, P_L, shots) -> (eps_L, err)`. They are tested end to
end in `tests/test_pipeline.py::test_scaffold_entry_points`.

Simulation baseline vs hardware: `experiments/compare_sim_hw.py results/hw/<run>.npz ...`
simulates the exact circuit and qubits of each run under the calibration model (and
the calibration model plus the ancilla errors learned from that run) and plots them
against the hardware points (`results/figures/sim_vs_hardware.png`).

## Quick start

```bash
python -m venv .venv && .venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python -m pytest -q                       # 24 tests, ~15 s, no token

# free, no token
.venv/Scripts/python experiments/simulate.py            # simulation baseline + figure
.venv/Scripts/python experiments/layout.py --offline    # SWAP-free placement (snapshot + fake backend)
.venv/Scripts/python experiments/run_hardware.py --fake --rounds 1 2 3 4 5 --shots 2000
.venv/Scripts/python experiments/analyze.py results/hw/<file>.npz

# real hardware: put IQM_TOKEN=... in .env (see .env.example)
.venv/Scripts/python experiments/layout.py              # live calibration → best patch
.venv/Scripts/python experiments/run_hardware.py --rounds 1 2 3 4 5 6 --shots 4000           # dry run
.venv/Scripts/python experiments/run_hardware.py --rounds 1 2 3 4 5 6 --shots 4000 --dd both --confirm
.venv/Scripts/python experiments/analyze.py results/hw/*.npz
```

## Key facts established so far

* **d = 3 fits Emerald with zero SWAPs** in 64 placements (16 qubit sets × 4 orientations);
  every stabilizer CZ is a native coupler. Transpiled depth: 11 (1 round), 29 (3 rounds).
* **d = 5 has no SWAP-free embedding** on Emerald's 54-qubit coupler graph (VF2 finds none):
  it would need routing, so we keep d = 3 on hardware and study d = 5 in simulation.
* The no-reset circuit is not worse than the reset one in simulation (no reset error, same
  detector semantics) — ε_L ≈ 0.35 %/round at uniform p = 0.3 %.

**Hardware results: see [RESULTS.md](RESULTS.md).** Headline: on IQM Emerald, ε_L = 2.8 % per
round with correction vs 8.6 % without (Z- and X-memory). Demo outline: `docs/demo_outline.md`.

## Acknowledgements

IQM's challenge scaffold (`PushQuantumQEC`) defined the interfaces. We read
[ItsSypher/Quantum-ECC](https://github.com/ItsSypher/Quantum-ECC) (an earlier
solution to this challenge) for orientation; none of its code is used here. The Emerald
coupler snapshot in `data/` is used only offline; live runs re-fetch topology and calibration.
