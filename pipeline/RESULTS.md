# Results — d = 3 rotated surface-code memory on IQM Emerald and Garnet

All hardware runs: distance-3 rotated surface code (9 data + 8 ancilla qubits), **no
mid-circuit reset** (detectors `m_r ⊕ m_{r-2}`), circuits in native H/PRX + CZ,
**0 SWAPs** (layout pinned, routing disabled, every CZ verified to be a coupler).
All round counts of a sweep were submitted in one batch, so every point of a curve
shares one calibration. Uncertainties are 1σ binomial (points) and Fisher-information
(fits). Raw per-shot data for every run is in `results/hw/`.

**Definitions.** *P_L* = fraction of shots whose final, decoded logical value differs
from the prepared one. *ε_L* = logical error per round, from a binomial-MLE fit of
`2 P_L = 1 − (1 − 2 ε_L)^n` over n = 1…5 (the README's formula). "No correction" =
the raw parity of the logical row of data qubits.

## Headline

> **IQM Emerald, DD on: ε_L = 2.65 ± 0.07 % per round (Z-memory) and 2.60 ± 0.07 %
> (X-memory) with BP+OSD, versus 8.6 % / 10.6 % without correction — error correction
> reduces the logical error rate 3–4× in both logical bases.** With PyMatching's
> correlated MWPM: 2.82 % / 2.76 %.

The decoder prior here is SI1000 (p = 1 %). Every prior/decoder combination we ran gives
2.6–3.5 % on this data set (tables below), so the conclusion does not depend on that choice.
A prior that takes readout from our own characterization (E5) improves this slightly:
2.62 ± 0.07 % (Z) and **2.46 ± 0.07 %** (X) with BP+OSD.

**Decoding pipeline:** replacing the generic readout model with our own readout
characterization (E5) gives a small but consistent improvement on both chips. Each per-run
gain is within error bars (Emerald X: 2.60 → 2.46 %, BP+OSD), but the sign is the same in
12–13 of 13 runs with the MWPM decoders (sign test p ≈ 0.0005–0.006) and 11/13 with
BP+OSD. Adding error rates learned from the QEC data (E7) helps only Garnet's noisy patch
(−1.1 to −1.3 pp; 2 runs on the same qubits). On Emerald they give no gain at 2 000 shots
per circuit.

P_L per round count (4 000 shots each):

| rounds n | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| no correction | 0.110 ± 0.005 | 0.169 ± 0.006 | 0.211 ± 0.006 | 0.240 ± 0.007 | 0.296 ± 0.007 |
| **decoded (correlated MWPM)** | **0.030 ± 0.003** | **0.052 ± 0.004** | **0.081 ± 0.004** | **0.106 ± 0.005** | **0.123 ± 0.005** |

Figure: `results/figures/20261121-151940_emerald_Z_ddon.png` (P_L vs n with fits, and the
per-stabilizer detector-firing map).

## All hardware runs (Z-memory)

ε_L per round:

| run | qubits | shots/pt | no correction | MWPM, calibrated prior | best decoder |
|---|---|---|---|---|---|
| Garnet, DD off | only SWAP-free set | 4 000 | 12.33 ± 0.21 % | 13.63 ± 0.22 % | 11.19 ± 0.19 % (learned, E7) |
| Garnet, DD on | same | 4 000 | 12.16 ± 0.20 % | 12.55 ± 0.21 % | 10.57 ± 0.18 % (learned, E7) |
| Emerald test, DD off | calibration-best patch | 1 000 (n = 1–3) | 6.82 ± 0.36 % | 7.88 ± 0.39 % | 5.98 ± 0.37 % (learned) |
| Emerald, DD off | **QB33, QB26 excluded** | 4 000 | 9.22 ± 0.16 % | 3.41 ± 0.08 % | 3.18 ± 0.08 % (MWPM) |
| **Emerald, DD on** | same | 4 000 | 8.59 ± 0.15 % | 3.15 ± 0.08 % | **2.82 ± 0.07 %** (correlated) |

Decoders compared on identical shots: MWPM and correlated MWPM (PyMatching) with
uniform, SI1000 and per-qubit calibrated priors, and our cross-validated data-learned
weights (`qec/learn.py`). Full tables: `results/figures/*.json`.

## X-memory (Emerald, DD on, same qubits)

Same 17 qubits as the headline Z run, logical |+⟩ stored, data read out in the X basis.

| rounds n | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| no correction | 0.138 ± 0.005 | 0.199 ± 0.006 | 0.256 ± 0.007 | 0.274 ± 0.007 | 0.336 ± 0.007 |
| **decoded (correlated MWPM)** | **0.026 ± 0.003** | **0.057 ± 0.004** | **0.077 ± 0.004** | **0.100 ± 0.005** | **0.124 ± 0.005** |

**ε_L = 2.76 ± 0.07 % per round** (correlated MWPM, SI1000 prior) vs 10.64 ± 0.18 % without
correction, i.e. a 3.9× reduction. This is consistent with Z-memory (2.82 ± 0.07 %), so both
logical bases are protected equally. The calibrated prior does worse here (4.01 %, MWPM) than
for Z (3.15 %). Figure: `results/figures/20261121-160049_emerald_X_ddon.png`.

## Decoder comparison (identical shots, identical priors)

`experiments/decoder_comparison.py`, figure `results/figures/decoder_comparison.png`.
ε_L per round with the SI1000 (p = 1 %) prior; µs/shot is Python-loop timing for the
`ldpc` decoders (an upper bound, not their asymptotic speed):

| decoder | Emerald Z, DD on | Emerald Z, DD off | Emerald X, DD on | Garnet Z, DD on | Garnet Z, DD off | µs/shot |
|---|---|---|---|---|---|---|
| none | 8.59 ± 0.15 | 9.22 ± 0.16 | 10.64 ± 0.18 | 12.16 ± 0.20 | 12.33 ± 0.21 | – |
| Union-Find | 3.30 ± 0.08 | 3.49 ± 0.08 | 3.22 ± 0.08 | 14.56 ± 0.24 | 15.71 ± 0.26 | ~50 |
| MWPM | 3.01 ± 0.08 | 3.40 ± 0.08 | 3.01 ± 0.08 | 12.17 ± 0.20 | 13.38 ± 0.22 | ~9 |
| correlated MWPM | 2.82 ± 0.07 | 3.24 ± 0.08 | 2.76 ± 0.07 | 11.62 ± 0.19 | 13.07 ± 0.21 | ~11 |
| belief-matching | 2.77 ± 0.07 | 3.20 ± 0.08 | 2.71 ± 0.07 | **11.39 ± 0.19** | 12.54 ± 0.21 | ~300 |
| **BP+OSD** | **2.65 ± 0.07** | **3.12 ± 0.08** | **2.60 ± 0.07** | 11.44 ± 0.19 | **12.26 ± 0.20** | ~450 |

* Ranking on Emerald is the same in every run: BP+OSD < belief-matching < correlated
  MWPM < MWPM < Union-Find. Using the full hyperedge error model (BP+OSD) gains ~12 %
  over plain MWPM, which matches the ordering in simulation.
* Union-Find is the fastest decoder asymptotically but costs ~10 % in ε_L. On Garnet
  (above threshold) it is worse than no correction.
* A calibration-weighted prior is not better than the generic SI1000 prior on these data
  (e.g. Emerald X: MWPM 4.01 % calibrated vs 3.01 % SI1000), because calibration
  under-states the ancilla errors that dominate (see Finding 2).

## Extensions

### E1. Repeated-measurement characterization predicts which ancillas fail
`experiments/characterize_measurement.py`: every qubit of Emerald (54) and Garnet (20)
prepared in |0⟩ or |1⟩ and measured 5× without reset (2 circuits × 4 000 shots per chip).
Per qubit: first-read errors, flip rates between consecutive reads, and real state changes
(an A-A-B-B pattern, which one readout error cannot fake; validated in Aer: a pure readout
error gives 0.3 %, a real decay 4.8 %). Data: `results/hw/meas_*.json`.

* **|1⟩ decays under measurement on every qubit**: 4–8 % per readout, against < 1 % for
  |0⟩ → |1⟩. An un-reset ancilla is in |1⟩ half the time, so this is a generic cost of the
  no-reset scheme that single-shot readout calibration does not show.
* **It predicts in-code ancilla errors**: across all 32 ancillas of the four Emerald regions
  below, the per-qubit flip rate correlates with the time-like error learned from QEC data
  with **Pearson r = 0.84**. Resonance's readout calibration reaches only r = 0.37. Emerald
  QB33 (excluded from our main patch after the test run) is the worst qubit on the chip
  (31 % decay per readout).
* Garnet QB11 is *not* anomalous under repeated measurement (flip < 3 %). Its ~25 % in-code
  error therefore comes from its gates or neighbours, not its readout.

### E2. Spatial variation across Emerald
`experiments/regions.py` picks the 3 qubit sets with least overlap with our main patch and
with each other (from 19 SWAP-free sets on today's coupler graph). Same circuits, DD on,
rounds 1–5, 2 000 shots per point (main patch: 4 000). Data: `results/figures/regions_emerald.json`.

| region | ε_L, BP+OSD | ε_L, correlated MWPM | hottest ancillas (learned) |
|---|---|---|---|
| main patch | **2.65 ± 0.07 %** | 2.82 ± 0.07 % | QB19 4.9 %, QB30 3.9 % |
| region 3 | 3.00 ± 0.11 % | 3.26 ± 0.11 % | QB12 7.5 %, QB19 5.0 % |
| region 2 | 3.20 ± 0.11 % | 3.49 ± 0.12 % | QB12 6.9 %, QB19 6.9 % |
| region 1 | 3.87 ± 0.12 % | 4.11 ± 0.13 % | QB33 10.8 %, QB10 8.8 % |

* The same code varies by 46 % in ε_L depending on its placement on the chip.
* The hottest ancilla of every region is one of the four worst qubits in E1 (QB33, QB10,
  QB19, QB12).
* Predicting the region ranking (`results/figures/regions_prediction.json`, 30 000
  simulated shots per point): the calibration-only model gets the order wrong (Spearman
  −0.2; it even ranks the main patch worst). Adding E1's flip rates gets it nearly right
  (Spearman +0.8). Absolute values are overestimated, and with 4 regions this is
  indicative, not conclusive.

### E3. Dynamical-decoupling strategy and reproducibility
IQM's server-side DD pass accepts a custom `DDStrategy`. Same calibration, same patch,
rounds 1–5 × 1 000 shots each (`results/figures/dd_strategy_comparison.json`):

| DD strategy | pulses on centre data qubit (5 rounds) | ε_L, BP+OSD | ε_L, correlated MWPM |
|---|---|---|---|
| off (compiled for comparison) | 14 | – | – |
| standard (one EDD per readout window) | 94 | 2.64 ± 0.14 % | 2.84 ± 0.15 % |
| dense (up to 3 EDD back to back) | 174 | 2.75 ± 0.14 % | 2.97 ± 0.15 % |

* Denser DD gives no further gain: standard DD already removes the idle dephasing. What
  remains (|1⟩ decay during readout, E1) is not something DD can remove.
* The standard-DD run reproduces the headline two hours later (2.64 % vs 2.65 %).

### E4. Pulse level (PulLA)
`experiments/run_pulla.py` converts our transpiled circuits with `qiskit_to_pulla`,
compiles them with PulLA's standard compiler (with the DD strategies of E3 set as compiler
stages) and submits the pulse playlist. Pulse submission needs an organisation token with
pulse-level access; a personal-account token is refused ("Personal account does not have
pulse-level access enabled").

Same calibration, same patch, rounds 1–5 × 1 000 shots each, submitted together
(`results/figures/pulse_vs_circuit.json`):

| path | ε_L, BP+OSD | ε_L, correlated MWPM | ε_L, no correction |
|---|---|---|---|
| circuit level (`backend.run`), standard DD | 2.61 ± 0.14 % | 2.85 ± 0.15 % | 8.27 % |
| pulse level (PulLA), standard DD | 2.79 ± 0.15 % | 3.00 ± 0.15 % | 8.74 % |
| pulse level (PulLA), dense DD | 2.56 ± 0.14 % | 2.71 ± 0.14 % | 7.91 % |

* The pulse-level path works end to end and reproduces the circuit-level result. Detector
  firing after round 1 is 0.078 vs 0.076, and all three rates agree within 1σ. This is
  expected: with the standard compiler, PulLA produces the same pulse schedule the server
  would.
* Dense DD at pulse level is the lowest point, but the difference is not significant here
  (≈ 1σ), consistent with E3.
* The compiled schedules give Emerald's real timings for our circuit: PRX 24 ns, CZ 40 ns,
  **readout 1.04 µs**, initial relaxation wait 400 µs. The mid-circuit readout is ~80 % of
  each syndrome round, so the remaining error budget is dominated by readout-time effects
  (E1). A pulse-level gain would therefore have to come from the readout itself (shorter
  integration, |1⟩-protecting readout), which needs recalibrated discrimination and was
  not attempted.

### E5. Decoder prior from our own readout characterization
`experiments/characterized_decoding.py`, data `results/figures/characterized_decoding.json`.
The E1 measurements were fed back into the decoder weights. Every model so far treated a
readout error as a single state flip before `M`. Here it is replaced, per qubit, by the two
effects E1 separates: a classical misread that corrupts only the recorded bit (Stim `M(p)`),
and a real state change during readout that persists in the no-reset circuit (`X_ERROR`
after `M`) (`NoiseModel.with_characterization`). E.g. ancilla QB19: 1.5 % misread, 4.1 %
state change. Four priors decode the same shots of all 12 hardware runs:

* **SI1000**: generic, p = 1 % (readout 5 %).
* **calibrated**: per-qubit and per-coupler Resonance calibration.
* **characterized**: calibrated, with readout from E1.
* **SI1000 + char.**: SI1000 for gates and idling, readout from E1.

ε_L per round (%), BP+OSD:

| run | SI1000 | calibrated | characterized | SI1000 + char. |
|---|---|---|---|---|
| Emerald Z, DD off | 3.12 ± 0.08 | 3.25 ± 0.08 | **3.08 ± 0.08** | 3.10 ± 0.08 |
| Emerald Z, DD on | 2.65 ± 0.07 | 2.88 ± 0.07 | 2.79 ± 0.07 | **2.62 ± 0.07** |
| Emerald X, DD on | 2.60 ± 0.07 | 3.55 ± 0.08 | 3.47 ± 0.08 | **2.46 ± 0.07** |
| Emerald region 1 | 3.87 ± 0.12 | 4.20 ± 0.13 | **3.86 ± 0.12** | **3.86 ± 0.12** |
| Emerald region 2 | 3.20 ± 0.11 | 3.41 ± 0.12 | 3.24 ± 0.11 | **3.03 ± 0.11** |
| Emerald region 3 | 3.00 ± 0.11 | 3.18 ± 0.11 | 3.18 ± 0.11 | **2.95 ± 0.11** |
| Emerald, dense DD | **2.75 ± 0.14** | 2.82 ± 0.15 | 2.80 ± 0.15 | 2.79 ± 0.15 |
| Emerald, circuit ref. | **2.61 ± 0.14** | 2.81 ± 0.15 | 2.74 ± 0.14 | **2.61 ± 0.14** |
| Emerald, PulLA standard | **2.79 ± 0.15** | 3.03 ± 0.15 | 2.87 ± 0.15 | 2.81 ± 0.15 |
| Emerald, PulLA dense | 2.56 ± 0.14 | 2.70 ± 0.14 | 2.71 ± 0.14 | **2.54 ± 0.14** |
| Garnet Z, DD off | 12.26 ± 0.20 | 12.19 ± 0.20 | 11.94 ± 0.20 | **11.81 ± 0.20** |
| Garnet Z, DD on | 11.44 ± 0.19 | 11.60 ± 0.19 | 11.27 ± 0.19 | **11.01 ± 0.18** |

Mean change vs SI1000 over the 12 runs (pp; negative = better) and runs where it is better:

| prior | MWPM | correlated MWPM | BP+OSD |
|---|---|---|---|
| calibrated | +0.19 (2/12) | +0.29 (0/12) | +0.23 (1/12) |
| characterized | +0.07 (5/12) | +0.15 (4/12) | +0.09 (4/12) |
| **SI1000 + char.** | **−0.15 (12/12)** | **−0.16 (11/12)** | **−0.11 (10/12, 1 tie)** |

* **Our readout characterization improves decoding when it replaces the generic readout
  model.** SI1000 + char. is the best of the four priors in 11/12 runs with MWPM and
  correlated MWPM and in 9/12 with BP+OSD. Each per-run gain is within the error bars, but
  the priors decode the same shots, and the gain has the same sign in 12/12 (MWPM) and
  11/12 (correlated) runs, which is very unlikely by chance (two-sided sign test p ≈ 0.0005 / 0.006).
* **New best results:** Emerald X 2.46 ± 0.07 % (BP+OSD; previously 2.60 %), Garnet
  11.01 % (DD on) and 11.81 % (DD off), previously 11.39 % and 12.26 %. Emerald Z
  (2.62 %) is unchanged within errors (2.65 %).
* **Fast decoders gain most.** Correlated MWPM with SI1000 + char. reaches 2.69 % / 2.60 %
  (Z / X) on Emerald, close to BP+OSD with SI1000 (2.65 % / 2.60 %) at ~40× lower cost
  (~11 vs ~450 µs/shot).
* **The calibrated gate/idle model is what hurts.** Swapping readout into the calibrated
  model helps only a little (calibrated → characterized: −0.12 to −0.14 pp), because the
  calibrated gate and idle terms cost more (up to ~0.9 pp on Emerald X). The
  characterization is useful; calibration as a decoder prior is not.
* The next step would be to characterize the gate and idle terms as well (or learn them
  from data, as the p_ij weights do) instead of taking them from SI1000 or calibration.

### E6. The circuits sent to IQM are the circuits we decode
`experiments/verify_circuits.py`, payloads `results/hw/payloads/`. For every hardware run
(76 circuits, 16 runs, Emerald and Garnet) the circuits are rebuilt exactly as submitted
(same code, layout, transpiler settings and seed; the rebuilt depth equals the depth saved
at submission for all 76) and serialized with IQM's own serializer
(`backend.serialize_circuit`, which `backend.run` uses to build the job payload). That
payload is parsed back into Stim independently of our Stim → Qiskit translation and
checked against the reference circuit the decoder uses:

| check | result |
|---|---|
| only `prx`, `cz`, `measure`, `barrier`; every CZ on a coupler; exactly the 17 planned qubits; every `prx` a Clifford | 76/76 |
| classical bit j measures the qubit the reference measures as record j (how the decoder reads the data) | 76/76 |
| noiseless, through the same conversion as the hardware data: no detector or observable fires | 76/76 |
| every single X, Y, Z error on every qubit at the start and after every round gives the same detection events and logical flip as the reference | 15 198/15 198 |

* **Negative test:** the checker catches a payload with one CZ removed, two classical bits
  swapped or one rotation changed. Removing a CZ still passes the noiseless test; only
  the error-response test catches it.
* **Checked against IQM's server:** the original payloads of the earlier jobs can no longer
  be downloaded with the current token ("job not found"), so those checks use rebuilt
  payloads. To close that gap, one new job was run on the same 17 qubits
  (`20261121-210819_emerald_Z_ddon`, DD on, rounds 1–5 × 2 000 shots). Its payload,
  downloaded back from IQM's server (`jobs/{id}/payload`), passes every check above and is
  **identical, instruction by instruction, to the rebuilt payload**, so the rebuild
  reproduces what the server receives. The server also reports the calibration set it
  compiled with (`ecec1d0a-…`) and `dd_mode` enabled. The job reproduces the earlier
  result on these qubits: ε_L = 2.70 ± 0.10 % (BP+OSD, SI1000; 2.60 ± 0.10 % with
  SI1000 + char.) vs 2.65 ± 0.07 % at 15:19 and 2.61 ± 0.14 % at 18:18.
* **Not covered:** IQM's server-side compilation after the payload (DD insertion, pulse
  scheduling). IQM's API returns only the calibration-set ID for a circuit job, not the
  compiled schedule. `run_hardware.py` now saves the exact payload of every submission.
* The hardware results support the same conclusion independently: a wrong circuit would
  fire detectors close to 50 % of the time, while the measured detection rates and
  logical error rates match the simulation (Finding 5).

### E7. Learning the remaining error rates from the QEC data
**Result:** on top of the characterized readout, learned error rates lower ε_L by
1.1–1.3 pp on Garnet (MWPM decoders; 2 runs, same patch) and give no gain on Emerald:
±0.05 pp with MWPM decoders, +0.15 pp (worse) with BP+OSD. Whether more training shots
would help on Emerald was not tested.

`experiments/combined_dem.py`, data `results/figures/combined_dem.json`. E5 left the gate
and idle terms at SI1000. Here every edge of the SI1000 + char. decoding graph is
re-weighted with its probability learned from the detection events themselves (p_ij
method, `qec.learn`), in two variants:

* **char + learned**: every edge learned except the time-like readout edges, which keep
  the E1 characterization.
* **all learned**: the readout edges are learned too.

Learned DEMs are 2-fold cross-validated (learned on half of the shots of every circuit,
scored on the other half, then swapped), so no shot is decoded with rates fit to it.
13 runs (the 12 of E5 plus the 21:08 job of E6).

**Pooling is required.** A first version learned each edge separately per circuit; with
~1 000 training shots per circuit it was *worse* than SI1000 + char. on Emerald by
0.6–0.8 pp. The final version pools every occurrence of the same physical edge over all
bulk rounds and all round counts of a run (`pooled_edge_probabilities`). In simulation
this turned no gain into a gain (9.65 % vs 9.80 % for the prior).

Change vs SI1000 + char. (pp, negative = better; runs where better):

| | MWPM | correlated MWPM | BP+OSD |
|---|---|---|---|
| Emerald (11 runs), char + learned | +0.05 (5/11) | +0.01 (5/11) | +0.14 (2/11) |
| Emerald (11 runs), all learned | +0.05 (5/11) | +0.00 (6/11) | +0.15 (0/11) |
| Garnet (2 runs), char + learned | −1.19 (2/2) | −1.12 (2/2) | −0.27 (2/2) |
| Garnet (2 runs), all learned | −1.18 (2/2) | **−1.29 (2/2)** | −0.37 (2/2) |

ε_L per round (%) on Garnet:

| run | SI1000, BP+OSD | SI1000 + char., BP+OSD | all learned, correlated | all learned, BP+OSD |
|---|---|---|---|---|
| Garnet Z, DD off | 12.26 ± 0.20 | 11.81 ± 0.20 | **11.19 ± 0.19** | 11.41 ± 0.19 |
| Garnet Z, DD on | 11.44 ± 0.19 | 11.01 ± 0.18 | **10.57 ± 0.18** | 10.68 ± 0.18 |

* **Garnet: learning the error rates helps substantially.** Learned weights cut ε_L by
  1.1–1.3 pp for the matching decoders and give new best results, 10.57 % (DD on) and
  11.19 % (DD off), previously 11.01 % / 11.81 %. With learned weights fast correlated
  MWPM beats BP+OSD with any fixed prior. This matches Finding 5: Garnet's errors are
  dominated by a few bad ancillas that neither SI1000 nor calibration describes.
* **Emerald: no gain.** For MWPM and correlated MWPM the learned DEMs tie SI1000 + char.
  (within ±0.05 pp, about half the runs better); BP+OSD loses ~0.15 pp. On this patch the
  SI1000 + characterized-readout prior already describes the noise about as well as
  2 000 shots per circuit can measure it. Pooling removed the earlier 0.6–0.8 pp loss
  but did not create a gain.
* **Learning the readout edges adds nothing over the characterization.** char + learned
  and all learned agree within 0.1 pp everywhere: the E1 characterization measures the
  readout terms as well as the QEC data does, without spending any QEC shots on it.
* Both learned variants still beat plain SI1000 with correlated MWPM in 13/13 runs
  (−0.32 / −0.35 pp mean), but on Emerald that gain comes from the characterized readout
  (E5), not from the learned gate/idle terms.
* BP+OSD gains least: it already infers hyperedge correlations that the learned graph
  weights only partly capture, and it is more sensitive to the noise of a learned rate.

## Findings

### 1. Error correction works on Emerald, not on Garnet
On Emerald the decoder lowers ε_L from 8.6 % to 2.8 %. On Garnet decoding with a fixed
prior barely beats raw readout (ε_L ≈ 11–12 % vs 12.2–12.3 %): Garnet's patch is close to or
above threshold for this circuit, so the decoder's corrections are almost as often wrong as
right. Error rates learned from the data (E7) give the clearest gain there (10.57 % vs
12.16 %, DD on), but still only ~1.15× vs 3× on Emerald.

### 2. Calibration underestimates specific ancillas by up to 30×
Time-like error per stabilizer (ancilla readout + idling), learned from the detection
events with the p_ij method, against what the calibration-based model expects:

| device | ancilla | learned (DD off) | learned (DD on) | calibration model |
|---|---|---|---|---|
| Garnet | Z@QB11 | 26.6 % | 24.8 % | 2.5 % |
| Garnet | X@QB7 | 21.3 % | 22.2 % | 3.8 % |
| Garnet | X@QB5 | 20.0 % | 9.8 % | 1.8 % |
| Garnet | X@QB15 | 12.2 % | 9.3 % | 2.8 % |
| Emerald (test) | X@QB33 | 14.5 % | – | 1.2 % |
| Emerald (test) | Z@QB26 | 10.7 % | – | 1.0 % |
| Emerald | X@QB39 | 18.3 % | 2.7 % | 0.6 % |
| Emerald | other 7 ancillas | 1.4–4.5 % | 2.0–4.9 % | 0.6–3.0 % |

The method was validated in simulation first: it recovers a planted 25 % ancilla error
as 25.1 % (`tests/test_pipeline.py::test_learned_weights_recover_hot_ancilla`).

### 3. Placement from measured errors halves the logical error rate
Emerald has 64 SWAP-free placements of the d = 3 code (16 qubit sets × 4 orientations;
Garnet has exactly one set). Re-placing the code to exclude the two ancillas measured as
bad in the test run (QB33, QB26) lowered ε_L, same decoder (MWPM, calibrated prior),
from 7.9 % to 3.4 % per round, and raised the fraction of syndrome-free shots after
2 rounds from 7 % to 28 %. (The test run had fewer shots and round counts, and a
different calibration day-time, so part of the gain may be drift.)

### 4. Dynamical decoupling separates three error mechanisms
Switching IQM's server-side DD on and off, on identical circuits and qubits:

* **idle dephasing — fixed by DD**: Emerald QB39 18.3 % → 2.7 %; Garnet QB5 20.0 % → 9.8 %.
* **energy relaxation — not fixed by DD**: Garnet QB7 (T1 = 9.6 µs, T2 = 5.2 µs) 21 % → 22 %.
* **neither (DD-insensitive, not in calibration)**: Garnet QB11 (calibrated readout 98.8 %,
  T1 = 51 µs) stays at ~25 %. Candidates: measurement-induced excitation of the un-reset
  ancilla, leakage, or a TLS defect.

Net effect of DD on ε_L (same decoder, correlated MWPM with SI1000 prior): Emerald
3.24 → 2.82 %, Garnet 13.07 → 11.62 %; fraction of syndrome-free shots after 2 rounds
Emerald 28 % → 38 %, Garnet 1.3 % → 3.6 %.

### 5. Simulation baseline vs hardware
`experiments/compare_sim_hw.py` simulates the exact circuit on the exact qubits of each
run (`results/figures/sim_vs_hardware.png`; all decoded with MWPM + calibrated prior):

| run | sim: calibration | sim: calibration + measured ancilla errors | hardware |
|---|---|---|---|
| Emerald DD off | 3.84 % | 5.33 % | 3.41 % |
| Emerald DD on | 3.84 % | 5.21 % | 3.15 % |
| Garnet DD off | 6.85 % | 11.44 % | 13.63 % |
| Garnet DD on | 6.85 % | 10.88 % | 12.55 % |

* On the Emerald patch the calibration model predicts the hardware to within ~15 %
  (hardware is slightly *better*). Adding the learned rates over-counts there, because on
  a good patch they mostly re-measure errors the model already contains.
* On Garnet the calibration model misses half the error; adding the measured ancilla
  errors closes ~70 % of the gap and reproduces P_L for n = 2…5 almost exactly
  (0.284 / 0.351 / 0.393 predicted vs 0.302 / 0.354 / 0.396 measured, n = 3…5).
  The bad ancillas are the dominant cause.
* Garnet's n = 1 point (0.143) is 2× above both models (0.064): an additional
  preparation/readout-boundary error not present in the calibration.

### 6. Comparison with prior work on the same chip
An earlier solution to this challenge
([ItsSypher/Quantum-ECC](https://github.com/ItsSypher/Quantum-ECC), Emerald, d = 3,
1 round) reports P_L = 0.364 (circuit level) and 0.050 (pulse level via PulLA). Ours:
**P_L = 0.030 ± 0.003 at 1 round, circuit level + DD**. Different calibration dates, so part
of the difference is the device.

## Simulation (no hardware)
`experiments/simulate.py`, the no-reset native circuit, 50 000 shots per point,
n = 1–10: ε_L = 0.35 % (MWPM) / 0.32 % (correlated) at uniform p = 0.3 %; no-reset is not
worse than reset (0.345 vs 0.353 %). On IQM's 54-qubit fake backend (noisy Aer), the full
hardware path (transpile → run → m2d → decode) gives ε_L ≈ 8.6 %.

## Limitations
* d = 5 has no SWAP-free embedding on Emerald (VF2 finds none), so no distance scaling
  (Λ) is shown; d = 3 alone cannot demonstrate operation below threshold.
* X-memory was run with DD on only.
* Data-learned error rates help only on Garnet. On Emerald, even pooled over rounds
  (E7), they tie SI1000 + characterized readout; more shots per circuit would be needed to
  learn the gate/idle terms more precisely than SI1000 already describes them.
* The headline picks the best of 11 prior/decoder combinations; the spread on that data
  set is 2.6–3.5 % (Union-Find excluded: 3.2–3.3 %).
* Pulse-level runs used PulLA's standard pulse shapes; no custom waveforms or readout
  changes were tried (E4).
* The simulation noise model still uses assumed durations (CZ 80 ns, readout 1.5 µs) rather
  than the compiled ones (40 ns, 1.04 µs) measured in E4.
* Garnet calibration does not report readout QND-ness, so the QB11 mechanism is not
  identified.
