"""Zero-SWAP placement of the surface code on the hardware coupler graph.

The code's interaction graph (ancilla-data pairs that share a CZ) is embedded
into the coupler graph by subgraph monomorphism (VF2). Every embedding found
is SWAP-free by construction. Among all of them we pick the one with the
lowest *predicted logical error*: each candidate gets its own calibration-
driven NoiseModel and is simulated in Stim and decoded with the matching
decoder - so the choice accounts for which physical qubit plays which role
(data vs X/Z ancilla), not just an average fidelity of the patch.
"""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np
from networkx.algorithms.isomorphism import GraphMatcher

from .circuit import build_memory_circuit
from .code import RotatedSurfaceCode
from .decode import decode
from .noise import NoiseModel


def code_graph(code: RotatedSurfaceCode) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from(range(code.num_qubits))
    g.add_edges_from(tuple(e) for e in code.interaction_edges())
    return g


def hardware_graph(edges, exclude=()) -> nx.Graph:
    g = nx.Graph()
    g.add_edges_from(edges)
    g.remove_nodes_from(exclude)
    return g


def find_embeddings(code: RotatedSurfaceCode, edges, exclude=(), limit: int = 100_000
                    ) -> list[dict[int, int]]:
    """All SWAP-free placements {code qubit -> physical index}."""
    gm = GraphMatcher(hardware_graph(edges, exclude), code_graph(code))
    out = []
    for m in gm.subgraph_monomorphisms_iter():      # hardware node -> code node
        out.append({c: h for h, c in m.items()})
        if len(out) >= limit:
            break
    return out


def is_native(code: RotatedSurfaceCode, layout: dict[int, int], edges) -> bool:
    es = {frozenset(e) for e in edges}
    return all(frozenset((layout[a], layout[b])) in es
               for a, b in (tuple(e) for e in code.interaction_edges()))


@dataclass
class ScoredLayout:
    layout: dict[int, int]
    predicted_p_l: float
    proxy_cost: float

    def physical_qubits(self, names: dict[int, str]) -> list[str]:
        return [names[self.layout[q]] for q in sorted(self.layout)]


def noise_for_layout(layout, calibration, names, default=None) -> NoiseModel:
    phys_names = {q: names[p] for q, p in layout.items()}
    return NoiseModel.from_calibration(calibration, phys_names, default=default)


def proxy_cost(code, layout, calibration, names) -> float:
    """Fast additive cost: summed error of every operation in one round."""
    nm = noise_for_layout(layout, calibration, names)
    cost = sum(nm.gate2(a, b) for a, b in (tuple(e) for e in code.interaction_edges()))
    cost += sum(nm.meas(q) for q in code.ancillas)
    cost += sum(2 * nm.gate1(q) for q in range(code.num_qubits))
    for q in code.data_qubits:
        ch = nm.idle(q, nm.dur_meas)
        cost += sum(ch) if ch else 0
    return float(cost)


def rank_layouts(
    code: RotatedSurfaceCode,
    edges,
    calibration: dict | None,
    names: dict[int, str],
    rounds: int = 3,
    basis: str = "Z",
    shots: int = 20_000,
    top_k: int = 16,
    exclude=(),
    default_noise: NoiseModel | None = None,
    seed: int = 0,
) -> list[ScoredLayout]:
    """Rank SWAP-free placements by simulated logical error under calibration.

    Without calibration all placements are equivalent; they are returned
    unranked (predicted_p_l = nan).
    """
    cands = find_embeddings(code, edges, exclude)
    if not cands:
        raise RuntimeError("no SWAP-free embedding of this code on the given coupler graph")
    if calibration is None:
        return [ScoredLayout(l, float("nan"), float("nan")) for l in cands]

    scored = sorted(((proxy_cost(code, l, calibration, names), l) for l in cands),
                    key=lambda t: t[0])
    out = []
    for cost, lay in scored[:top_k]:
        nm = noise_for_layout(lay, calibration, names, default_noise)
        circ = build_memory_circuit(code.distance, rounds, basis, True, nm, code=code)
        det, obs = circ.compile_detector_sampler(seed=seed).sample(shots, separate_observables=True)
        out.append(ScoredLayout(lay, decode(det, obs, circ).p_logical, cost))
    out.sort(key=lambda s: s.predicted_p_l)
    return out


def draw_layout(code: RotatedSurfaceCode, layout: dict[int, int], names: dict[int, str]) -> str:
    d = code.distance
    rows = []
    for y in range(2 * d + 1):
        cells = []
        for x in range(2 * d + 1):
            q = next((q for q, c in code.coords.items() if c == (x, y)), None)
            if q is None:
                cells.append("      .")
            else:
                kind = "D" if q in code.data_qubits else next(
                    s.basis for s in code.stabilizers if s.ancilla == q)
                cells.append(f"{kind}:{names[layout[q]]:>5}")
        rows.append(" ".join(cells))
    return "\n".join(rows)
