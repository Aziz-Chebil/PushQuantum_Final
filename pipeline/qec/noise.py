"""Circuit-level noise models for simulation and for decoder weights.

A NoiseModel assigns an error strength to every operation *per qubit / per
pair*, so it can be driven either by one uniform number (for baselines) or by
the live Emerald calibration of the exact physical qubits the code runs on.
The decoder's detector error model is built from the same object, which means
calibration data directly sets the matching-graph edge weights.

Conversions from benchmarked fidelities to Pauli-channel strengths
------------------------------------------------------------------
* Average gate fidelity F of an n-qubit depolarizing channel with total Pauli
  probability p:  F = 1 - p * 2^n / (2^n + 1)  ->  p = (1 - F) (2^n + 1) / 2^n.
  So p1 = 1.5 (1 - F_1q) and p2 = 1.25 (1 - F_cz).
* Readout fidelity f -> symmetric flip probability 1 - f.
* Idling for time t with T1, T2 -> Pauli-twirled amplitude + phase damping:
  px = py = (1 - e^{-t/T1}) / 4,  pz = (1 - e^{-t/T2}) / 2 - px.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


def _pair(a, b):
    return frozenset((a, b))


@dataclass
class NoiseModel:
    p1: float = 0.0                    # default 1q depolarizing (after H / PRX)
    p2: float = 0.0                    # default 2q depolarizing (after CZ)
    p_meas: float = 0.0                # default readout flip
    p_reset: float = 0.0               # default init error (X flip after R)
    # Per-qubit / per-pair overrides, keyed by circuit (dense) qubit index.
    p1_q: dict[int, float] = field(default_factory=dict)
    p2_pair: dict[frozenset, float] = field(default_factory=dict)
    p_meas_q: dict[int, float] = field(default_factory=dict)
    p_reset_q: dict[int, float] = field(default_factory=dict)
    # Measured readout model (from experiments/characterize_measurement.py).
    # p_meas above is a state flip BEFORE the measurement, which in a no-reset
    # circuit persists into later rounds. Real readout has two separate parts:
    # a classical misread that corrupts only the recorded bit (M(p) in Stim),
    # and a real state change during the readout (|1> decay), which persists.
    p_meas_record_q: dict[int, float] = field(default_factory=dict)
    p_meas_state_q: dict[int, float] = field(default_factory=dict)
    # Idle noise. t1/t2 in microseconds; None disables idling errors.
    t1_q: dict[int, float] = field(default_factory=dict)
    t2_q: dict[int, float] = field(default_factory=dict)
    t1: float | None = None
    t2: float | None = None
    # Flat depolarizing charged to qubits idling through a measurement window
    # (used by SI1000; adds to any T1/T2 idle noise).
    p_idle_meas: float = 0.0
    # Operation durations in microseconds (IQM-typical defaults; override from
    # the backend target when connected).
    dur_1q: float = 0.04
    dur_2q: float = 0.08
    dur_meas: float = 1.5

    # ------------------------------------------------------------ factories
    @classmethod
    def uniform(cls, p: float) -> "NoiseModel":
        """Stim-style uniform circuit noise (every channel strength = p)."""
        return cls(p1=p, p2=p, p_meas=p, p_reset=p)

    @classmethod
    def si1000(cls, p: float) -> "NoiseModel":
        """Superconducting-inspired model (Gidney & Fowler 2021): cheap 1q
        gates, expensive readout, idle noise during measurement."""
        return cls(p1=p / 10, p2=p, p_meas=5 * p, p_reset=2 * p, p_idle_meas=2 * p)

    @classmethod
    def emerald_typical(cls) -> "NoiseModel":
        """Order-of-magnitude Emerald numbers, used for anything the
        calibration does not cover: 1q 0.05%, CZ 0.5%, readout 2%, init 1%,
        T1 = 40 us, T2 = 30 us."""
        return cls(p1=5e-4, p2=5e-3, p_meas=2e-2, p_reset=1e-2, t1=40.0, t2=30.0)

    @classmethod
    def from_calibration(
        cls,
        calibration: dict,
        qubit_names: dict[int, str],
        default: "NoiseModel | None" = None,
        reset_error_scale: float = 0.5,
    ) -> "NoiseModel":
        """Build a per-qubit model from Resonance calibration data.

        Parameters
        ----------
        calibration : dict in the iqm-qubit-selector format
            {"CZ": {"['QB1', 'QB2']": F}, "1Q": {"QB1": F}, "readout": {...},
             "t1": {"QB1": us}, "t2": {...}}
        qubit_names : circuit qubit index -> physical name ("QB17")
        default     : fallback model for missing entries
        reset_error_scale : initial-state error as a fraction of the readout
            error (thermal population is not reported separately).
        """
        base = default or cls.emerald_typical()
        m = cls(p1=base.p1, p2=base.p2, p_meas=base.p_meas, p_reset=base.p_reset,
                t1=base.t1, t2=base.t2, p_idle_meas=base.p_idle_meas,
                dur_1q=base.dur_1q, dur_2q=base.dur_2q, dur_meas=base.dur_meas)
        name_to_q = {v: k for k, v in qubit_names.items()}

        for name, f in calibration.get("1Q", {}).items():
            if name in name_to_q:
                m.p1_q[name_to_q[name]] = 1.5 * (1 - f)
        for name, f in calibration.get("readout", {}).items():
            if name in name_to_q:
                q = name_to_q[name]
                m.p_meas_q[q] = 1 - f
                m.p_reset_q[q] = reset_error_scale * (1 - f)
        for key, f in calibration.get("CZ", {}).items():
            a, b = _parse_pair(key)
            if a in name_to_q and b in name_to_q:
                m.p2_pair[_pair(name_to_q[a], name_to_q[b])] = 1.25 * (1 - f)
        for name, t in calibration.get("t1", {}).items():
            if name in name_to_q:
                m.t1_q[name_to_q[name]] = t
        for name, t in calibration.get("t2", {}).items():
            if name in name_to_q:
                m.t2_q[name_to_q[name]] = t
        return m

    # ------------------------------------------------------------- lookups
    def gate1(self, q: int) -> float:
        return self.p1_q.get(q, self.p1)

    def gate2(self, a: int, b: int) -> float:
        return self.p2_pair.get(_pair(a, b), self.p2)

    def meas(self, q: int) -> float:
        return self.p_meas_q.get(q, self.p_meas)

    def reset(self, q: int) -> float:
        return self.p_reset_q.get(q, self.p_reset)

    def meas_record(self, q: int) -> float:
        return self.p_meas_record_q.get(q, 0.0)

    def meas_state(self, q: int) -> float:
        return self.p_meas_state_q.get(q, 0.0)

    def with_characterization(self, rows: list[dict], qubit_names: dict[int, str],
                              floor: float = 1e-3) -> "NoiseModel":
        """Copy of this model whose readout comes from the repeated-measurement
        characterization (rows of characterize_measurement.analyze).

        Per qubit: state_flip (a change confirmed before and after) becomes a
        persisting state flip during readout; the remaining part of the
        consecutive-outcome flip rate, flip ~ 2 p_record + p_state, becomes a
        classical misread. The calibration-based pre-measurement flip is
        removed for those qubits, since the two parts replace it.
        """
        from copy import deepcopy

        m = deepcopy(self)
        by_name = {r["qubit"]: r for r in rows}
        for q, name in qubit_names.items():
            r = by_name.get(name)
            if r is None:
                continue
            p_state = max(r["state_flip"], floor)
            p_record = max((r["flip"] - p_state) / 2, floor)
            m.p_meas_state_q[q] = p_state
            m.p_meas_record_q[q] = p_record
            m.p_meas_q[q] = 0.0
        return m

    def idle(self, q: int, duration_us: float) -> tuple[float, float, float] | None:
        """Pauli-twirled (px, py, pz) for idling `duration_us`, or None."""
        t1 = self.t1_q.get(q, self.t1)
        t2 = self.t2_q.get(q, self.t2)
        if t1 is None or t2 is None or duration_us <= 0:
            return None
        t2 = min(t2, 2 * t1)                 # physical bound
        px = (1 - math.exp(-duration_us / t1)) / 4
        pz = max((1 - math.exp(-duration_us / t2)) / 2 - px, 0.0)
        if px + pz < 1e-12:
            return None
        return (px, px, pz)


def _parse_pair(key) -> tuple[str, str]:
    if isinstance(key, (tuple, list)):
        return key[0], key[1]
    parts = [p.strip(" '\"") for p in str(key).strip("[]()").split(",")]
    return parts[0], parts[1]
