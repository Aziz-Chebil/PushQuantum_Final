"""Rotated surface code geometry.

Coordinates follow the usual convention: data qubits sit at odd (x, y) and
measurement ancillas at even (x, y). Every ancilla touches the (up to) four
data qubits diagonally adjacent to it, so after a 45 degree rotation the code's
interaction graph is a patch of the square lattice - which is exactly the
coupler graph of IQM Emerald. That is what makes a SWAP-free embedding possible.

Qubits are indexed densely: data qubits 0..d^2-1 first, then ancillas.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

# CX schedules (offset of the data qubit relative to the ancilla, per layer).
# X and Z ancillas use mirrored "N" / "Z" orders. This makes every pair of
# overlapping stabilizers commute under interleaving, and it points the hook
# errors (an ancilla fault after two of the four CX gates spreads to two data
# qubits) perpendicular to the logical operator they would otherwise shorten.
# tests/test_circuit.py checks both properties with Stim.
Z_ORDER = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
X_ORDER = [(1, 1), (-1, 1), (1, -1), (-1, -1)]


@dataclass(frozen=True)
class Stabilizer:
    ancilla: int          # dense qubit index of the measurement ancilla
    basis: str            # "X" or "Z"
    coord: tuple[int, int]
    data: tuple[int, ...]  # dense data-qubit indices in its support


@dataclass
class RotatedSurfaceCode:
    distance: int
    coords: dict[int, tuple[int, int]] = field(init=False)
    data_qubits: list[int] = field(init=False)
    ancillas: list[int] = field(init=False)
    stabilizers: list[Stabilizer] = field(init=False)

    def __post_init__(self):
        d = self.distance
        if d < 3 or d % 2 == 0:
            raise ValueError(f"distance must be odd and >= 3, got {d}")

        self.coords = {}
        coord_to_q = {}
        for y in range(d):
            for x in range(d):
                q = len(self.coords)
                c = (2 * x + 1, 2 * y + 1)
                self.coords[q] = c
                coord_to_q[c] = q
        self.data_qubits = list(range(d * d))

        self.ancillas = []
        self.stabilizers = []
        for y in range(0, 2 * d + 1, 2):
            for x in range(0, 2 * d + 1, 2):
                basis = "X" if ((x + y) // 2) % 2 == 0 else "Z"
                on_lr = x in (0, 2 * d)
                on_tb = y in (0, 2 * d)
                if on_lr and on_tb:
                    continue                    # corners carry no stabilizer
                if on_lr and basis == "X":
                    continue                    # left/right are Z boundaries
                if on_tb and basis == "Z":
                    continue                    # top/bottom are X boundaries
                support = tuple(
                    coord_to_q[(x + dx, y + dy)]
                    for dx, dy in (X_ORDER if basis == "X" else Z_ORDER)
                    if (x + dx, y + dy) in coord_to_q
                )
                q = len(self.coords)
                self.coords[q] = (x, y)
                self.ancillas.append(q)
                self.stabilizers.append(Stabilizer(q, basis, (x, y), support))

        assert len(self.stabilizers) == d * d - 1

    # ------------------------------------------------------------------ props
    @property
    def num_qubits(self) -> int:
        return len(self.coords)

    def stabilizers_of(self, basis: str) -> list[Stabilizer]:
        return [s for s in self.stabilizers if s.basis == basis]

    @cached_property
    def logical(self) -> dict[str, tuple[int, ...]]:
        """Minimum-weight logical operator supports {"X": (...), "Z": (...)}.

        A logical Z must commute with every X stabilizer (even overlap) and
        vice versa; the candidates are the straight rows and columns of data.
        """
        d = self.distance
        rows = [tuple(y * d + x for x in range(d)) for y in range(d)]
        cols = [tuple(y * d + x for y in range(d)) for x in range(d)]

        def commutes(support, other_basis):
            s = set(support)
            return all(len(s & set(st.data)) % 2 == 0
                       for st in self.stabilizers_of(other_basis))

        out = {}
        for basis, other in (("Z", "X"), ("X", "Z")):
            out[basis] = next(c for c in rows + cols if commutes(c, other))
        return out

    def cx_layers(self) -> list[list[tuple[int, int]]]:
        """Four layers of (control, target) pairs for one syndrome round.

        X ancillas control CX onto data; Z ancillas are CX targets of data.
        """
        layers = []
        for k in range(4):
            layer = []
            for st in self.stabilizers:
                x, y = st.coord
                dx, dy = (X_ORDER if st.basis == "X" else Z_ORDER)[k]
                c = (x + dx, y + dy)
                q = next((q for q in st.data if self.coords[q] == c), None)
                if q is None:
                    continue
                layer.append((st.ancilla, q) if st.basis == "X" else (q, st.ancilla))
            layers.append(layer)
        return layers

    def interaction_edges(self) -> set[frozenset[int]]:
        """All qubit pairs that need a two-qubit gate (for hardware embedding)."""
        return {frozenset(p) for layer in self.cx_layers() for p in layer}

    def draw(self) -> str:
        """ASCII picture of the patch (D = data, X/Z = ancilla type)."""
        d = self.distance
        grid = [["  ." for _ in range(2 * d + 1)] for _ in range(2 * d + 1)]
        for q in self.data_qubits:
            x, y = self.coords[q]
            grid[y][x] = f"D{q:<2}"[:3].rjust(3)
        for st in self.stabilizers:
            x, y = st.coord
            grid[y][x] = f"{st.basis}{st.ancilla}".rjust(3)
        return "\n".join(" ".join(r) for r in grid)
