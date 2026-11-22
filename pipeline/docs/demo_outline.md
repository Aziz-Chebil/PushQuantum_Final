# 3-minute demo video — outline

Target: ~3:00, screen recording + voice. Numbers are from `RESULTS.md`.

## 0:00–0:20  The task
- One logical qubit in a d = 3 rotated surface code (17 qubits) on IQM Emerald,
  syndrome rounds without mid-circuit reset, decoded offline, logical error rate per round.
- Headline on screen: **ε_L = 2.8 % per round with correction vs 8.6 % without (3× better), in both Z and X memory.**

## 0:20–0:55  Our pipeline (show `README.md` table / repo tree)
- Own circuit generator written directly in Emerald's native gates (H + CZ), with
  Hadamard-frame tracking so redundant gates cancel.
- No-reset detectors `m_r ⊕ m_{r-2}`, which Stim proves correct (26 tests).
- Zero-SWAP placement: VF2 finds all 64 embeddings on Emerald, each scored by a
  calibrated simulation. Show the `layout.py` output: `swaps: 0, non_native_pairs: 0`.

## 0:55–1:40  The main result (show `results/figures/20261121-151940_emerald_Z_ddon.png`)
- Left panel: P_L vs rounds. Blue = no correction, others = decoders. Fit gives ε_L.
- Read one number: after 1 round only 3 % of shots are wrong (0.030), against 11 % raw.
- Mention X-memory gives the same rate (2.76 %).
- Comparison: prior work on this chip reported P_L = 0.364 (circuit) / 0.050 (pulse) after 1 round.

## 1:40–2:30  What limits the error rate (show heatmap + learned-rate table)
- Learned each ancilla's error from the data (p_ij): calibration underestimates some by up to 30×.
- Excluding the two bad ancillas found on Emerald halved ε_L (7.9 % → 3.4 %).
- DD on/off separates mechanisms: QB39 dephasing (18 % → 3 %, DD fixes), Garnet QB7
  T1 decay (DD can't fix), Garnet QB11 unexplained (not in calibration).
- Show `sim_vs_hardware.png`: on Garnet, adding the measured bad ancillas to the model
  reproduces the hardware curve.

## 2:30–3:00  Garnet, limits, wrap-up
- Garnet: ε_L ≈ 12 %, above threshold, so decoding does not help. Same pipeline, a different chip, one command line.
- Limits: d = 5 does not fit SWAP-free on Emerald, so no distance scaling yet.
- Close: reproducible from raw shots (`results/hw/*.npz`), and everything reruns with
  `experiments/*.py`.
