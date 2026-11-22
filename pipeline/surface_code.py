"""Challenge-scaffold entry points, implemented on top of the `qec` package.

The IQM scaffold (surface_code.py, run_on_hardware.py, build_emerald_qubit_rotated.py)
defines these names and signatures; they are provided here so the pipeline can
be driven exactly as the challenge README describes:

    from surface_code import (make_stim_circuit, make_qiskit_circuit, DEFAULT_NOISE,
                              build_emerald_qubit_map, print_qubit_map,
                              run_qiskit, decode_hardware_results, fit_logical_error_rate)

    circ = make_stim_circuit(3, rounds=3, basis="Z", no_reset=True, noise=None)
    syndromes = run_qiskit(distance=3, rounds=3, shots=4000)           # needs IQM_TOKEN
    ler, err = decode_hardware_results(syndromes)

Track: Qiskit (Stim circuit -> Qiskit -> IQM Resonance -> m2d -> PyMatching).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import stim
from qiskit import QuantumCircuit

from qec import NoiseModel, RotatedSurfaceCode, build_memory_circuit
from qec.analysis import fit_ler
from qec.decode import decode, measurements_to_detections, syndrome_diagnostics
from qec.hardware import (coupling_edges, get_backend, job_to_arrays, qubit_names,
                          stim_to_qiskit, submit, swap_report, transpile_pinned)
from qec.layout import draw_layout, find_embeddings, is_native

RESULTS = Path(__file__).resolve().parent / "results"

#: Same keys as the scaffold's DEFAULT_NOISE (Stim's four generator knobs).
DEFAULT_NOISE = dict(
    after_clifford_depolarization=0.003,
    after_reset_flip_probability=0.003,
    before_measure_flip_probability=0.003,
    before_round_data_depolarization=0.003,
)


def _to_noise_model(noise) -> NoiseModel | None:
    """Accept None, a NoiseModel, or a scaffold-style dict of Stim knobs."""
    if noise is None or isinstance(noise, NoiseModel):
        return noise
    p_gate = noise.get("after_clifford_depolarization", 0.0)
    return NoiseModel(
        p1=p_gate, p2=p_gate,
        p_reset=noise.get("after_reset_flip_probability", 0.0),
        p_meas=noise.get("before_measure_flip_probability", 0.0),
        # depolarization of the data qubits once per round = our idle channel
        # charged to every qubit that waits through the ancilla measurement
        p_idle_meas=noise.get("before_round_data_depolarization", 0.0),
    )


def _distance_of(circuit: stim.Circuit) -> int:
    # our circuits use 2 d^2 - 1 qubits (d^2 data + d^2 - 1 ancillas)
    return int(round(math.sqrt((circuit.num_qubits + 1) / 2)))


# ---------------------------------------------------------------- 1. circuits
def make_stim_circuit(distance: int, rounds: int, basis: str = "Z",
                      no_reset: bool = True, noise=None) -> stim.Circuit:
    """Rotated surface-code memory experiment in IQM-native gates (H + CZ).

    noise : None (hardware circuit), a qec.NoiseModel (per-qubit, e.g. from
            calibration), or a dict with the scaffold's DEFAULT_NOISE keys.
    no_reset=True rewrites the detectors to m_r XOR m_(r-2) (see qec/circuit.py).
    """
    return build_memory_circuit(distance, rounds, basis, no_reset, _to_noise_model(noise))


def make_qiskit_circuit(distance: int, rounds: int, basis: str = "Z",
                        no_reset: bool = True) -> QuantumCircuit:
    """The same experiment as a Qiskit circuit (clbit j = Stim record j)."""
    return stim_to_qiskit(make_stim_circuit(distance, rounds, basis, no_reset))


def make_qrisp_experiment(*args, **kwargs):
    raise NotImplementedError("This solution follows the Qiskit track; "
                              "use make_stim_circuit / make_qiskit_circuit.")


def get_circuit_info(circuit: stim.Circuit) -> dict:
    code = RotatedSurfaceCode(_distance_of(circuit))
    return {
        "num_qubits": circuit.num_qubits,
        "num_data_qubits": len(code.data_qubits),
        "num_anc_qubits": len(code.ancillas),
        "num_measurements": circuit.num_measurements,
        "num_detectors": circuit.num_detectors,
        "num_observables": circuit.num_observables,
        "circuit_distance": len(circuit.detector_error_model(decompose_errors=True)
                                .shortest_graphlike_error())
        if circuit.num_detectors and _has_noise(circuit) else None,
    }


def _has_noise(circuit: stim.Circuit) -> bool:
    return any(i.name in ("DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR", "PAULI_CHANNEL_1")
               for i in circuit.flattened())


# ------------------------------------------------------- 4. qubit placement
def build_emerald_qubit_map(stim_circuit: stim.Circuit, quantum_computer: str = "emerald",
                            layout_file: str | Path | None = None) -> dict[int, int]:
    """Circuit qubit index -> physical (0-based Qiskit) qubit index, SWAP-free.

    Uses the calibration-ranked placement saved by experiments/layout.py
    (results/layout_<qc>.json) when present; otherwise the first SWAP-free
    embedding on the cached coupler graph.
    """
    path = Path(layout_file or RESULTS / f"layout_{quantum_computer}.json")
    if path.exists():
        info = json.loads(path.read_text())
        return {int(k): v for k, v in info["layout"].items()}
    code = RotatedSurfaceCode(_distance_of(stim_circuit))
    embs = find_embeddings(code, coupling_edges(None, quantum_computer), limit=1)
    if not embs:
        raise RuntimeError("no SWAP-free placement on this device")
    return embs[0]


def print_qubit_map(stim_circuit: stim.Circuit, qubit_map: dict[int, int] | None = None,
                    quantum_computer: str = "emerald", edges=None):
    """Table of the placement plus a check that every CZ pair is a coupler."""
    qubit_map = qubit_map or build_emerald_qubit_map(stim_circuit, quantum_computer)
    edges = edges or coupling_edges(None, quantum_computer)
    names = qubit_names()
    code = RotatedSurfaceCode(_distance_of(stim_circuit))
    role = {s.ancilla: f"{s.basis}-anc" for s in code.stabilizers}
    print(f"{'circuit q':>9} {'coord':>8} {'role':>6} {'qiskit idx':>10} {'Resonance':>9}")
    for q in sorted(qubit_map):
        print(f"{q:>9} {str(code.coords[q]):>8} {role.get(q, 'data'):>6} "
              f"{qubit_map[q]:>10} {names[qubit_map[q]]:>9}")
    pairs = [tuple(e) for e in code.interaction_edges()]
    es = {frozenset(e) for e in edges}
    bad = [(a, b) for a, b in pairs if frozenset((qubit_map[a], qubit_map[b])) not in es]
    print(draw_layout(code, qubit_map, names))
    print(f"CZ pairs: {len(pairs) - len(bad)} native, {len(bad)} non-native (need SWAP)")
    return not bad


# ------------------------------------------------------- 5. hardware execution
def run_qiskit(distance: int = 3, rounds: int = 1, shots: int = 5_000, token: str | None = None,
               quantum_computer: str = "emerald", no_reset: bool = True,
               basis: str = "Z", dd: bool = False) -> dict:
    """Stim -> Qiskit -> pinned zero-SWAP transpile -> IQM Resonance -> syndromes.

    Returns the scaffold's syndromes dict (det_events, obs_flips, raw_meas,
    num_detectors, num_shots, syndrome_weight_per_shot, detector_firing_rate)
    plus "stim_circuit", "rounds", "basis" so decode_hardware_results can build
    a matching DEM. For multi-round sweeps with persistence use
    experiments/run_hardware.py.
    """
    circ = make_stim_circuit(distance, rounds, basis, no_reset, noise=None)
    backend = get_backend(quantum_computer, token)
    layout = build_emerald_qubit_map(circ, quantum_computer)
    assert is_native(RotatedSurfaceCode(distance), layout, coupling_edges(backend))
    tqc = transpile_pinned(stim_to_qiskit(circ), backend, layout)
    rep = swap_report(tqc, coupling_edges(backend))
    print(f"transpiled depth {rep['depth']}, {rep['two_qubit_gates']} CZ, {rep['swaps']} SWAP")
    job = submit(backend, [tqc], shots, dd=dd)
    raw = job_to_arrays(job, [circ.num_measurements])[0]
    return extract_syndromes(raw, circ) | {"stim_circuit": circ, "rounds": rounds,
                                          "basis": basis, "job_id": str(job.job_id())}


def extract_syndromes(raw_meas: np.ndarray, stim_circuit: stim.Circuit) -> dict:
    det, obs = measurements_to_detections(raw_meas, stim_circuit)
    diag = syndrome_diagnostics(det)
    return {"det_events": det, "obs_flips": obs, "raw_meas": raw_meas,
            "num_detectors": det.shape[1], "num_shots": len(det),
            "syndrome_weight_per_shot": det.sum(axis=1),
            "detector_firing_rate": diag["detector_firing_rate"],
            "stim_circuit": stim_circuit}


# ------------------------------------------------------- 6. offline decoding
def decode_hardware_results(syndromes: dict, noise=None, decoder: str = "mwpm"
                            ) -> tuple[float, float]:
    """Decode with PyMatching; return (P_L, 1-sigma error) for this round count.

    noise   : decoder prior (NoiseModel or scaffold dict); default uniform p=1%.
    decoder : "mwpm", "correlated" or "none".
    The DEM is built from the circuit that ran, with the prior's noise added,
    so its detectors match the data by construction.
    """
    circ = syndromes["stim_circuit"]
    d = _distance_of(circ)
    rounds = syndromes.get("rounds") or _rounds_of(circ, d)
    basis = syndromes.get("basis", "Z")
    prior = _to_noise_model(noise) or NoiseModel.uniform(0.01)
    noisy = build_memory_circuit(d, rounds, basis, True, prior) \
        if _is_no_reset(circ) else build_memory_circuit(d, rounds, basis, False, prior)
    res = decode(syndromes["det_events"], syndromes["obs_flips"], noisy, decoder)
    return res.p_logical, res.stderr


def fit_logical_error_rate(rounds, p_logical, shots, free_amplitude: bool = False):
    """Fit 2 P_L = 1 - (1 - 2 eps_L)^n over several round counts.
    Returns (eps_L, 1-sigma error)."""
    errors = np.rint(np.asarray(p_logical) * np.asarray(shots)).astype(int)
    f = fit_ler(rounds, errors, shots, free_amplitude=free_amplitude)
    return f.eps, f.eps_err


def _rounds_of(circ: stim.Circuit, d: int) -> int:
    n_anc = d * d - 1
    return (circ.num_measurements - d * d) // n_anc


def _is_no_reset(circ: stim.Circuit) -> bool:
    return sum(1 for i in circ.flattened() if i.name == "R") == 1
