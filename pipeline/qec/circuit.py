"""Surface-code memory-experiment circuits in Stim, written in IQM's native gate set.

Why not stim.Circuit.generated?
-------------------------------
1. Native gates. Emerald executes PRX and CZ, not CX. We compile every CX to
   H.CZ.H on the target *while* tracking a Hadamard frame per qubit, so that
   back-to-back Hadamards cancel and only the Hadamards the hardware really
   executes appear in the circuit. The noisy Stim model (and therefore the
   decoder) then sees exactly the operations that run on the chip.
2. No-reset detectors. Emerald has no mid-circuit reset. Without reset an
   ancilla's raw outcome is m_r = s_r XOR m_{r-1}, so the detector comparing
   syndromes of rounds r and r-1 is m_r XOR m_{r-2} (with m_0 = 0), and the
   final data-vs-ancilla detector uses s_R = m_R XOR m_{R-1}. Stim simulates
   the un-reset ancilla physically, so the determinism check in the tests is a
   real proof that these detector definitions are right.
3. Per-qubit noise. Every gate, idle window and readout can carry its own
   calibrated error strength (see noise.py).

Hadamard-frame bookkeeping
--------------------------
We keep a set F of qubits whose physical state equals H applied to the
logical state. A logical H just toggles membership (no gate). A logical
CX(c, t) needs t in F and c not in F, then runs as one physical CZ. A
Z-basis measurement needs the qubit outside F. Physical H gates are emitted
only to repair those preconditions.
"""

from __future__ import annotations

import stim

from .code import RotatedSurfaceCode
from .noise import NoiseModel


class _Builder:
    def __init__(self, n_qubits: int, noise: NoiseModel | None, native: str):
        self.c = stim.Circuit()
        self.n = n_qubits
        self.noise = noise
        self.native = native
        self.frame: set[int] = set()
        self.num_meas = 0

    # ---------------------------------------------------------------- noise
    def _idle(self, busy: set[int], duration: float):
        if self.noise is None:
            return
        for q in range(self.n):
            if q in busy:
                continue
            ch = self.noise.idle(q, duration)
            if ch is not None:
                self.c.append("PAULI_CHANNEL_1", [q], list(ch))
            if duration >= self.noise.dur_meas and self.noise.p_idle_meas > 0:
                self.c.append("DEPOLARIZE1", [q], self.noise.p_idle_meas)

    # ---------------------------------------------------------- physical ops
    def h(self, qubits):
        qubits = sorted(qubits)
        if not qubits:
            return
        self.c.append("H", qubits)
        if self.noise is not None:
            for q in qubits:
                if self.noise.gate1(q) > 0:
                    self.c.append("DEPOLARIZE1", [q], self.noise.gate1(q))
            self._idle(set(qubits), self.noise.dur_1q)
        self.c.append("TICK")

    def two_qubit(self, gate: str, pairs):
        flat = [q for p in pairs for q in p]
        self.c.append(gate, flat)
        if self.noise is not None:
            for a, b in pairs:
                if self.noise.gate2(a, b) > 0:
                    self.c.append("DEPOLARIZE2", [a, b], self.noise.gate2(a, b))
            self._idle(set(flat), self.noise.dur_2q)
        self.c.append("TICK")

    def reset(self, qubits):
        qubits = sorted(qubits)
        self.c.append("R", qubits)
        self.frame -= set(qubits)
        if self.noise is not None:
            for q in qubits:
                if self.noise.reset(q) > 0:
                    self.c.append("X_ERROR", [q], self.noise.reset(q))
        self.c.append("TICK")

    def measure(self, qubits, duration_idle: bool = True) -> dict[int, int]:
        """Z-basis measurement; returns {qubit: absolute record index}."""
        qubits = sorted(qubits)
        if self.native == "cz":
            self.h(self.frame & set(qubits))
        self.frame -= set(qubits)
        if self.noise is not None:
            for q in qubits:
                if self.noise.meas(q) > 0:
                    self.c.append("X_ERROR", [q], self.noise.meas(q))
        if self.noise is not None and any(self.noise.meas_record(q) > 0 for q in qubits):
            # Per-qubit classical misread; same record order as one M.
            for q in qubits:
                self.c.append("M", [q], self.noise.meas_record(q))
        else:
            self.c.append("M", qubits)
        if self.noise is not None:
            for q in qubits:
                if self.noise.meas_state(q) > 0:
                    self.c.append("X_ERROR", [q], self.noise.meas_state(q))
        if self.noise is not None and duration_idle:
            self._idle(set(qubits), self.noise.dur_meas)
        rec = {q: self.num_meas + i for i, q in enumerate(qubits)}
        self.num_meas += len(qubits)
        self.c.append("TICK")
        return rec

    # ----------------------------------------------------------- logical ops
    def logical_h(self, qubits):
        if self.native == "cz":
            self.frame ^= set(qubits)
        else:
            self.h(qubits)

    def logical_cx_layer(self, pairs):
        if self.native == "cx":
            self.two_qubit("CX", pairs)
            return
        fix = set()
        for ctrl, tgt in pairs:
            if ctrl in self.frame:
                fix.add(ctrl)
            if tgt not in self.frame:
                fix.add(tgt)
        self.h(fix)
        self.frame ^= fix
        self.two_qubit("CZ", pairs)

    # ------------------------------------------------------------ annotations
    def rec(self, absolute: int) -> stim.GateTarget:
        return stim.target_rec(absolute - self.num_meas)

    def detector(self, records, coord):
        self.c.append("DETECTOR", [self.rec(r) for r in records], coord)


def build_memory_circuit(
    distance: int = 3,
    rounds: int = 3,
    basis: str = "Z",
    no_reset: bool = True,
    noise: NoiseModel | None = None,
    native: str = "cz",
    code: RotatedSurfaceCode | None = None,
) -> stim.Circuit:
    """Rotated surface-code memory experiment.

    Parameters
    ----------
    distance : odd code distance (3 on a 17-qubit Emerald patch)
    rounds   : number of syndrome-extraction rounds R >= 1
    basis    : "Z" (store |0>_L) or "X" (store |+>_L)
    no_reset : True -> ancillas are never reset (required on Emerald)
    noise    : None for the hardware circuit, a NoiseModel for simulation/DEM
    native   : "cz" (H + CZ, what Emerald runs) or "cx" (textbook CX circuit)
    """
    basis = basis.upper()
    if basis not in ("X", "Z"):
        raise ValueError("basis must be 'X' or 'Z'")
    if rounds < 1:
        raise ValueError("rounds must be >= 1")
    if native not in ("cz", "cx"):
        raise ValueError("native must be 'cz' or 'cx'")

    code = code or RotatedSurfaceCode(distance)
    b = _Builder(code.num_qubits, noise, native)
    for q, (x, y) in code.coords.items():
        b.c.append("QUBIT_COORDS", [q], [x, y])

    x_anc = [s.ancilla for s in code.stabilizers_of("X")]
    layers = code.cx_layers()

    # Initialisation. Qubits start in |0>; X memory prepares |+> on the data.
    b.reset(range(code.num_qubits))
    if basis == "X":
        b.logical_h(code.data_qubits)

    history: dict[int, list[int]] = {s.ancilla: [] for s in code.stabilizers}

    def syndrome_records(anc: int, r: int) -> list[int]:
        """Records whose XOR is the syndrome s_r (r is 1-based)."""
        m = history[anc]
        if r < 1:
            return []
        if no_reset and r >= 2:
            return [m[r - 1], m[r - 2]]
        return [m[r - 1]]

    for r in range(1, rounds + 1):
        b.logical_h(x_anc)
        for layer in layers:
            b.logical_cx_layer(layer)
        b.logical_h(x_anc)
        rec = b.measure(code.ancillas)
        for anc in code.ancillas:
            history[anc].append(rec[anc])
        if not no_reset:
            b.reset(code.ancillas)

        for st in code.stabilizers:
            if r == 1 and st.basis != basis:
                continue           # first-round opposite-basis syndromes are random
            recs = syndrome_records(st.ancilla, r) + syndrome_records(st.ancilla, r - 1)
            recs = _xor_reduce(recs)
            x, y = st.coord
            b.detector(recs, [x, y, r - 1])

    # Final transversal data readout in the memory basis.
    if basis == "X":
        b.logical_h(code.data_qubits)
    data_rec = b.measure(code.data_qubits, duration_idle=False)

    for st in code.stabilizers_of(basis):
        recs = [data_rec[q] for q in st.data] + syndrome_records(st.ancilla, rounds)
        x, y = st.coord
        b.detector(_xor_reduce(recs), [x, y, rounds])

    b.c.append("OBSERVABLE_INCLUDE",
               [b.rec(data_rec[q]) for q in code.logical[basis]], 0)
    return b.c


def _xor_reduce(records: list[int]) -> list[int]:
    """Drop records that appear an even number of times (they cancel)."""
    out = []
    for r in records:
        if r in out:
            out.remove(r)
        else:
            out.append(r)
    return sorted(out)


def describe(circuit: stim.Circuit) -> dict:
    """Gate counts and sizes - handy for reports and sanity checks."""
    counts: dict[str, int] = {}
    for inst in circuit.flattened():
        if inst.name in ("H", "M", "R"):
            counts[inst.name] = counts.get(inst.name, 0) + len(inst.targets_copy())
        elif inst.name in ("CZ", "CX"):
            counts[inst.name] = counts.get(inst.name, 0) + len(inst.targets_copy()) // 2
    return {
        "qubits": circuit.num_qubits,
        "measurements": circuit.num_measurements,
        "detectors": circuit.num_detectors,
        "observables": circuit.num_observables,
        "gate_counts": counts,
    }
