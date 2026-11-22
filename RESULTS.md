# Results

All results are from **IQM Emerald hardware** unless marked as simulation. Distance-3
rotated surface code, Z memory unless stated, no mid-circuit reset, zero SWAPs, DD on.
Detailed write-ups per project: [pipeline/RESULTS.md](pipeline/RESULTS.md) (experiments
E1–E7), [circuit_optimization/OPTIMIZATION.md](circuit_optimization/OPTIMIZATION.md),
[qrisp_track/RESULTS.md](qrisp_track/RESULTS.md).

## 1. All runs, one fitter

The three projects originally reported ε_L with different fit models: amplitude fixed at 1,
or free up to 2 (which lets the fit imply a negative error probability at zero rounds).
`analysis/refit_all.py` refits every run's failure counts with the same binomial maximum
likelihood fit, `P_L(n) = [1 − (1 − 2ε)^n] / 2`, over rounds 1–5, and again from round 2
(removes the start-of-circuit effect). All rows: BP+OSD, SI1000 prior with p = 1 %.
Numbers in this table supersede the per-project write-ups where they differ.

| circuit | qubits | run | shots/pt | P_L at n = 1 / 2 / 3 / 4 / 5 (%) | ε_L, rounds 1–5 | ε_L, rounds 2–5 |
|---|---|---|---|---|---|---|
| starting circuit (pipeline) | pipeline patch | Z, 21 Nov 15:19 | 4 000 | 2.88 / 5.47 / 7.78 / 9.68 / 11.30 | 2.65 ± 0.07 % | 2.64 ± 0.07 % |
| starting circuit (pipeline) | pipeline patch | Z, 21 Nov 21:08, cal. `ecec1d0a` | 2 000 | 3.55 / 4.85 / 7.40 / 10.10 / 11.75 | 2.70 ± 0.10 % | 2.62 ± 0.10 % |
| starting circuit (pipeline) | pipeline patch | X, 21 Nov 16:00 | 4 000 | 2.80 / 4.90 / 7.25 / 9.72 / 11.75 | 2.60 ± 0.07 % | 2.59 ± 0.07 % |
| optimized circuit, Qiskit | pipeline patch | run 1, cal. `ecec1d0a` | 4 000 | 1.20 / 3.83 / 5.58 / 7.75 / 10.68 | 2.04 ± 0.06 % | 2.11 ± 0.06 % |
| optimized circuit, Qrisp | pipeline patch | run 1 | 4 000 | 1.07 / 3.23 / 5.70 / 7.75 / 9.95 | 1.94 ± 0.06 % | 2.01 ± 0.06 % |
| **optimized circuit, Qiskit** | pipeline patch | **4 campaigns pooled**, cal. `ecec1d0a` | 16 000 | 1.01 / 3.98 / 6.02 / 8.09 / 10.93 | **2.12 ± 0.03 %** | 2.21 ± 0.03 % |
| optimized circuit, Qrisp | pipeline patch | 4 campaigns pooled | 16 000 | 1.01 / 4.08 / 6.11 / 8.79 / 11.11 | 2.20 ± 0.03 % | 2.29 ± 0.03 % |
| Qrisp track | own patch, not recorded | SI1000 re-decode | 4 000 | 1.18 / 3.40 / 6.17 / 7.45 / 9.20 | 1.92 ± 0.06 % | 1.98 ± 0.06 % |
| Qrisp track | own patch, not recorded | twirl study, control | 2 000 | 1.25 / 2.40 / 5.35 / 7.75 / 8.80 | 1.79 ± 0.08 % | 1.83 ± 0.08 % |

Without correction (raw logical readout) the pipeline patch gives 8.59 ± 0.15 % per round
(Z) and 10.64 ± 0.18 % (X).

**Reading the table:**

* **The optimized circuit beats the starting circuit on the same 17 qubits and the same
  calibration set** (`ecec1d0a`): 2.12 ± 0.03 % vs 2.70 ± 0.10 %, also from round 2
  (2.21 vs 2.62 %). The runs are separate jobs, so drift is not excluded, but the gap is
  ~5σ while repeated runs of the starting circuit agree within 0.05 pp (2.65 / 2.70 %).
  The circuits differ in rotation count (94 vs 104 at 3 rounds), CZ layer order and
  barrier placement; the data do not say which matters most.
* **Its advantage is largest at the start** (P_L at one round 1.0 % vs 2.9–3.6 %) but
  remains per round (rounds 2–5 fit).
* **Qiskit and Qrisp compile the optimized circuit to the same native gate counts** and give
  the same ε_L within 0.1 pp: the SDK does not matter, the circuit and the qubits do.
* **The Qrisp-track pipeline has the lowest numbers (1.79–1.92 %) but cannot be ranked**:
  its qubits were chosen automatically and not saved, it ran in a different session, its
  decoder prior approximates SI1000 through Stim's generated-circuit options, and its
  compiler removes the 5 round-1 CZs that act on qubits still in |0⟩ (19 vs 24 CZs at one
  round). A rerun on the pipeline patch would be needed to compare it.

### Matched job: both circuits in one job

`circuit_optimization/results/project_comparison/20261121T200524815997Z`. Both circuits in the
same hardware job, pipeline patch, DD on, 4 000 shots each, BP+OSD with SI1000:

| rounds | optimized circuit | starting circuit |
|---|---|---|
| 3 | 236 failures (5.90 %) | 325 (8.13 %) |
| 5 | 456 (11.40 %) | 520 (13.00 %) |

The 3-round difference is ~3.9σ; at 5 rounds ~2.2σ. In this job the starting circuit did
worse at 5 rounds (13.0 %) than in its own runs (11.3 %, 11.75 %).

## 2. What changed the logical error rate

| change | effect | evidence |
|---|---|---|
| qubit placement avoiding ancillas measured to be bad | ε_L 7.9 → 3.4 % (MWPM, calibrated prior) | pipeline Finding 3; separate runs, partly drift possible |
| layer barriers (stop IQM's compiler from interleaving layers) | failures at 3 rounds 28.7 % → 10.1 %, at 5 rounds 44.4 % → 16.1 % | circuit_optimization `results/translation_check`; same patch and calibration, different times |
| dynamical decoupling | ε_L 3.24 → 2.82 % (Emerald), 13.07 → 11.62 % (Garnet) | pipeline Finding 4 |
| optimized circuit | 2.70 → 2.12 % | section 1 |
| decoder: MWPM → BP+OSD | 3.01 → 2.65 % | pipeline decoder comparison, identical shots |
| decoder prior: SI1000 → SI1000 + characterized readout | −0.11 to −0.16 pp mean, better in 11–13 of 13 runs | pipeline E5, identical shots, each gain within error bars |
| learned error rates (p_ij), on top of that | Garnet −1.1 to −1.3 pp; Emerald none (BP+OSD +0.15 pp worse) | pipeline E7, 2-fold cross-validated |
| Pauli twirling | no consistent effect | section 3 |

## 3. Pauli twirling

`circuit_optimization/results/twirl_emerald_20261122*`: optimized circuit, pipeline patch,
20 twirled variants interleaved with 20 control blocks per round count, Paulis merged into
existing rotations (gate counts unchanged). Both runs pooled, 20 000 shots per round count per
arm:

| pipeline | control failures | twirled failures | difference | ε_L control → twirled |
|---|---|---|---|---|
| Qrisp | 6 084 | 5 859 | −225 (≈ 2.1σ) | 2.15 → 2.06 % |
| Qiskit | 5 965 | 6 093 | +128 (≈ 1.2σ) | 2.10 → 2.15 % |

Twirling helped slightly with Qrisp and hurt slightly with Qiskit on circuits with the same
native gate counts; opposite signs of 1–2σ are consistent with run-to-run variation.
**Conclusion: no consistent effect; not adopted.** The Qrisp-track prototype
(`qrisp_track/twirling.py`) inserted the Paulis as extra gates, ran control and twirled
circuits one after the other, and found ε_L 1.79 → 2.39 %; with the added gates and no
interleaving, that run does not test twirling itself.

## 4. Decoder and prior (pipeline, identical shots)

* BP+OSD < belief-matching < correlated MWPM < MWPM < Union-Find on every Emerald run.
  Correlated MWPM with SI1000 + characterized readout comes within ~0.1 pp of BP+OSD at ~40×
  lower cost.
* A prior from the device calibration is worse than generic SI1000 (+0.19 to +0.29 pp).
  The repeated-measurement characterization is useful as the readout part of the prior;
  calibration gate and idle errors are not.
* The approximate likelihood decoder (`circuit_optimization/likelihood_decoder.py`) together
  with the optimized circuit reduced failures 17 % in simulation versus forward order + MWPM.
  On hardware, with its own device-based prior, it matched BP+OSD with SI1000 on the same
  shots (5.80 vs 5.90 % at 3 rounds, 11.55 vs 11.40 % at 5 rounds).

## 5. Verification

* `pipeline/experiments/verify_circuits.py` (E6): 76 submitted circuits rebuilt and serialized
  with IQM's serializer, parsed back into Stim and checked: native gates on couplers only,
  correct bit order, noiseless detectors silent, and identical detector response to every
  single X/Y/Z fault (15 198/15 198). One job's payload downloaded from IQM's server is
  identical to the rebuild.
* `circuit_optimization/verify_qec.py`, `verify_optimization.py`, `investigate_translation.py`:
  noiseless parities, graphlike distance, single-fault and logical-string response, IQM
  measurement order after serialization; twirled variants pass ideal parity audits.
* The optimized circuit has its own audits but has **not** been run through the E6 per-fault
  check against IQM payloads.

## 6. Garnet (second device, pipeline)

Only one SWAP-free d = 3 patch exists on Garnet, and it is near threshold: best ε_L 10.57 ±
0.18 % (DD on, learned error rates, correlated MWPM) vs 12.16 % without correction.

## 7. Simulation (no hardware)

* Optimized circuit + likelihood decoder at balanced p = 0.3 %: ε_L = 0.243 % per round, 17 %
  fewer failures than forward order + MWPM (800 000 validation shots); distance 5 and 7 lower
  further at p = 0.3 % (36 / 13 / 3 failures in 2 000 shots at 8 rounds).
* These are simulation assumptions; hardware is ~9× worse.

## Limitations

* d = 3 only: Emerald has no zero-SWAP d = 5 embedding, so no distance scaling (Λ) on
  hardware.
* The best circuit and the best decoder prior have not been combined yet: the optimized
  circuit was decoded with SI1000 only, not SI1000 + characterized readout.
* Comparisons across jobs share a calibration set but not time; drift is not excluded.
* The constant-rate fit is rejected by goodness-of-fit for the large runs (round-dependent
  rate, largest at round 2); compare per-round P_L when it matters.
* Qrisp-track qubits were not recorded, so its runs cannot be compared with the others.

## Next steps

1. Decode the optimized circuit's 16 000-shot data with SI1000 + characterized readout
   (offline, raw shots are saved).
2. Run the optimized circuit through the E6 payload check.
3. One interleaved job: optimized vs starting circuit at rounds 1–5, and the Qrisp-track
   circuit on the same patch.
