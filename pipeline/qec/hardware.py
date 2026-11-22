"""Everything that touches IQM Resonance: backend, topology, calibration,
transpilation with a pinned zero-SWAP layout, job submission and readout.

Only the functions that take a `backend` need a token. Transpilation against a
fake/offline target and all verification helpers work without one.
"""

from __future__ import annotations

import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath

import numpy as np
import stim
from qiskit import QuantumCircuit, transpile

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
SERVER_URL = "https://resonance.iqm.tech"


# ------------------------------------------------------------------ backend
def get_token(token: str | None = None) -> str:
    """Token from the argument, $IQM_TOKEN, or a git-ignored .env file."""
    if token:
        return token
    if os.environ.get("IQM_TOKEN"):
        return os.environ["IQM_TOKEN"]
    env = DATA_DIR.parent / ".env"
    if env.exists():
        for line in env.read_text().splitlines():
            if line.strip().startswith("IQM_TOKEN="):
                return line.split("=", 1)[1].strip().strip("'\"")
    raise RuntimeError("No IQM Resonance token: pass token=, set IQM_TOKEN, or "
                       "put IQM_TOKEN=... in QEC-Emerald/.env")


def get_backend(quantum_computer: str = "emerald", token: str | None = None):
    from iqm.qiskit_iqm import IQMProvider

    # The SDK rejects a token given both as an argument and as an environment
    # variable, and the calibration fetcher only reads the environment - so we
    # pass it exclusively via IQM_TOKEN.
    os.environ["IQM_TOKEN"] = get_token(token)
    provider = IQMProvider(SERVER_URL, quantum_computer=quantum_computer)
    return provider.get_backend()


# --------------------------------------------------------- topology / calib
def coupling_edges(backend=None, quantum_computer: str = "emerald") -> list[tuple[int, int]]:
    """Undirected coupler list (0-based). Live if a backend is given, else the
    cached snapshot in data/."""
    if backend is not None:
        edges = {tuple(sorted(e)) for e in backend.coupling_map.get_edges()}
        return sorted(edges)
    snap = DATA_DIR / f"{quantum_computer}_coupling_snapshot.json"
    return [tuple(e) for e in json.loads(snap.read_text())["coupling_map"]]


def qubit_names(backend=None, n: int = 54) -> dict[int, str]:
    if backend is not None:
        return {i: backend.index_to_qubit_name(i) for i in range(backend.num_qubits)}
    return {i: f"QB{i + 1}" for i in range(n)}


def fetch_calibration(backend, save: bool = True) -> dict:
    """Latest calibration fidelities via iqm-qubit-selector (dict of dicts:
    CZ, CLIFFORD, 1Q, readout, readout_qndness, t1, t2). Cached to data/."""
    from iqm.qubit_selector.qubit_selector import CalibrationDataManager

    cal = CalibrationDataManager().get_calibration_fidelities(backend)
    cal = {k: {str(kk): float(vv) for kk, vv in v.items()} for k, v in cal.items()}
    if save:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        path = DATA_DIR / f"calibration_{backend.name}_{stamp}.json"
        path.write_text(json.dumps({"fetched_at": stamp, "calibration": cal}, indent=1))
        (DATA_DIR / f"calibration_{backend.name}_latest.json").write_text(path.read_text())
    return cal


def calibration_from_fake(backend) -> dict:
    """Calibration-format dict from an IQM fake backend's error profile, so the
    calibration-aware code paths can be exercised offline against the same
    noise the fake backend simulates. Aer depolarizing parameter lam ->
    average fidelity F = 1 - lam (d-1)/d."""
    ep = backend.error_profile
    cal = {"1Q": {}, "CZ": {}, "readout": {}, "t1": {}, "t2": {}}
    for q, lam in ep.single_qubit_gate_depolarizing_error_parameters["prx"].items():
        cal["1Q"][q] = 1 - lam / 2
    for (a, b), lam in ep.two_qubit_gate_depolarizing_error_parameters["cz"].items():
        cal["CZ"][str([a, b])] = 1 - 3 * lam / 4
    for q, errs in ep.readout_errors.items():
        cal["readout"][q] = 1 - (errs["0"] + errs["1"]) / 2
    for q, t in ep.t1s.items():
        cal["t1"][q] = t / 1000.0            # ns -> us
    for q, t in ep.t2s.items():
        cal["t2"][q] = t / 1000.0
    return cal


def load_calibration(quantum_computer: str = "emerald") -> dict | None:
    path = DATA_DIR / f"calibration_{quantum_computer}_latest.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())["calibration"]


# ---------------------------------------------------------- stim -> qiskit
def stim_to_qiskit(circuit: stim.Circuit, barriers: bool = True) -> QuantumCircuit:
    """Convert our native (H, CZ, M, R) Stim circuit to Qiskit.

    * Qubit i of the Stim circuit becomes Qiskit qubit i (our circuits are
      already dense), and measurement record j is written to clbit j, so the
      bitstring order equals Stim's measurement record order.
    * Noise and annotations are dropped.
    * Leading R instructions are dropped (qubits start in |0> after the
      device's own initialisation); any later R becomes a Qiskit reset and is
      reported, because Emerald does not support mid-circuit reset.
    * TICKs that follow a measurement become barriers, so each syndrome round
      stays a separate block and the transpiler cannot slide gates across a
      mid-circuit measurement window.
    """
    qc = QuantumCircuit(circuit.num_qubits, circuit.num_measurements)
    m = 0
    started = False
    last_was_meas = False
    for inst in circuit.flattened():
        name = inst.name
        qs = [t.value for t in inst.targets_copy() if t.is_qubit_target]
        if name in ("QUBIT_COORDS", "DETECTOR", "OBSERVABLE_INCLUDE", "SHIFT_COORDS",
                    "DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR", "Y_ERROR", "Z_ERROR",
                    "PAULI_CHANNEL_1", "PAULI_CHANNEL_2"):
            continue
        if name == "TICK":
            if barriers and last_was_meas:
                qc.barrier()
            last_was_meas = False
            continue
        if name == "R":
            if started:
                for q in qs:
                    qc.reset(q)
            continue
        started = True
        last_was_meas = False
        if name == "H":
            for q in qs:
                qc.h(q)
        elif name == "CZ":
            for a, b in zip(qs[::2], qs[1::2]):
                qc.cz(a, b)
        elif name == "CX":
            for a, b in zip(qs[::2], qs[1::2]):
                qc.cx(a, b)
        elif name == "M":
            for q in qs:
                qc.measure(q, m)
                m += 1
            last_was_meas = True
        else:
            raise NotImplementedError(f"stim_to_qiskit: unsupported instruction {name}")
    assert m == circuit.num_measurements
    return qc


# ------------------------------------------------------------ transpilation
def transpile_pinned(qc: QuantumCircuit, backend, layout: dict[int, int],
                     optimization_level: int = 3, seed: int = 7) -> QuantumCircuit:
    """Transpile with the code pinned to `layout` and routing disabled.

    routing_method="none" makes Qiskit raise instead of inserting a SWAP, so a
    successful return is itself the proof that the circuit is SWAP-free.
    """
    init = [layout[q] for q in range(qc.num_qubits)]
    return transpile(qc, backend=backend, initial_layout=init,
                     layout_method="trivial", routing_method="none",
                     optimization_level=optimization_level, seed_transpiler=seed)


def swap_report(tqc: QuantumCircuit, edges) -> dict:
    """Count SWAPs and two-qubit gates on non-coupled pairs in a transpiled
    circuit (both must be zero for a native embedding)."""
    edge_set = {frozenset(e) for e in edges}
    # Physical index = position in tqc.qubits. (find_bit can return the
    # virtual index for IQM-transpiled circuits, so we do not rely on it.)
    pos = {q: i for i, q in enumerate(tqc.qubits)}
    n2 = non_native = swaps = 0
    for inst in tqc.data:
        if inst.operation.name == "swap":
            swaps += 1
        if inst.operation.num_qubits == 2 and inst.operation.name != "barrier":
            n2 += 1
            pair = frozenset(pos[q] for q in inst.qubits)
            if pair not in edge_set:
                non_native += 1
    ops = dict(tqc.count_ops())
    return {"two_qubit_gates": n2, "swaps": swaps, "non_native_pairs": non_native,
            "depth": tqc.depth(), "ops": ops}


def payload_dict(backend, tqc: QuantumCircuit) -> dict:
    """The circuit exactly as IQM receives it: serialized with the backend's
    own serializer (the one backend.run uses), as a JSON-ready dict."""
    c = backend.serialize_circuit(tqc)
    return {"name": c.name,
            "instructions": [{"name": i.name, "locus": list(i.locus), "args": i.args}
                             for i in c.instructions]}


# --------------------------------------------------------------- execution
def dd_strategy(kind: str):
    """DD strategies for IQM's server-side DD pass.

    standard : IQM's default (CPMG / XY4 / one EDD per idle window)
    dense    : up to 3 back-to-back EDD sequences (24 pulses), centre-aligned,
               in long idle windows such as the ~1 us mid-circuit readout
    """
    from iqm.station_control.interface.models.circuit import DDStrategy

    if kind == "standard":
        return None
    if kind == "dense":
        edd = "XYXYYXYX"
        return DDStrategy(gate_sequences=[(25, edd * 3, "center"), (17, edd * 2, "center"),
                                          (9, edd, "center"), (5, "YXYX", "center"),
                                          (2, "XX", "center")])
    raise ValueError(kind)


def dd_options(enabled: bool, strategy: str = "standard"):
    """CircuitCompilationOptions with IQM's server-side dynamical decoupling."""
    from iqm.iqm_client import CircuitCompilationOptions, DDMode

    return CircuitCompilationOptions(dd_mode=DDMode.ENABLED if enabled else DDMode.DISABLED,
                                     dd_strategy=dd_strategy(strategy) if enabled else None)


def submit(backend, tqcs: list[QuantumCircuit], shots: int, dd: bool = False,
           dd_strategy_name: str = "standard"):
    return backend.run(tqcs, shots=shots,
                       circuit_compilation_options=dd_options(dd, dd_strategy_name))


def memory_to_array(memory: list[str], num_measurements: int) -> np.ndarray:
    """Per-shot Qiskit bitstrings -> (shots, num_measurements) bool array in
    clbit order (Qiskit strings are little-endian: last char = clbit 0)."""
    clean = [s.replace(" ", "") for s in memory]
    arr = np.array([list(s[::-1]) for s in clean], dtype="<U1") == "1"
    if arr.shape[1] != num_measurements:
        raise ValueError(f"expected {num_measurements} clbits, got {arr.shape[1]}")
    return arr


def counts_to_array(counts: dict[str, int], num_measurements: int) -> np.ndarray:
    """Fallback when per-shot memory is unavailable (shot order is lost, which
    is fine for decoding)."""
    rows = []
    for bits, c in counts.items():
        row = memory_to_array([bits], num_measurements)[0]
        rows.append(np.repeat(row[None, :], c, axis=0))
    return np.concatenate(rows, axis=0)


def job_to_arrays(job, circuits_num_meas: list[int], timeout: float | None = None
                  ) -> list[np.ndarray]:
    result = job.result(timeout=timeout) if timeout is not None else job.result()
    out = []
    for i, nm in enumerate(circuits_num_meas):
        try:
            out.append(memory_to_array(result.get_memory(i), nm))
        except Exception:
            out.append(counts_to_array(result.get_counts(i), nm))
    return out


def save_run(path: Path, **arrays_and_meta):
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {k: v for k, v in arrays_and_meta.items() if not isinstance(v, np.ndarray)}
    arrays = {k: v for k, v in arrays_and_meta.items() if isinstance(v, np.ndarray)}
    np.savez_compressed(path, meta=json.dumps(meta, default=str), **arrays)


def load_run(path: Path) -> dict:
    z = np.load(path, allow_pickle=False)
    out = {k: z[k] for k in z.files if k != "meta"}
    out.update(json.loads(str(z["meta"])))
    cal = out.get("calibration_file")
    if cal:
        # Runs store the absolute path of the machine that took them; resolve it
        # by file name in data/ so the analysis works from any clone.
        # PureWindowsPath splits on both "\" and "/" on every OS.
        local = DATA_DIR / PureWindowsPath(cal).name
        out["calibration_file"] = str(local if local.exists() else cal)
    return out


def timestamp() -> str:
    return time.strftime("%Y%m%d-%H%M%S")
