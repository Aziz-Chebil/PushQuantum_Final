"""Check that the circuits sent to IQM are the circuits we decode.

For every hardware run, the circuits are rebuilt exactly as run_hardware.py
built them (same code, layout, transpiler settings and seed) and serialized
with IQM's own serializer (backend.serialize_circuit, which backend.run uses
to build the job payload). That payload is then parsed back into a Stim
circuit, independently of stim_to_qiskit, and checked against the reference
Stim circuit the decoder uses:

  1. structure   only prx / cz / measure / barrier; every CZ on a coupler;
                 exactly the planned physical qubits; every prx a Clifford
                 (angles multiples of pi/2), so Stim can simulate it exactly
  2. readout     clbit j is a measurement of the qubit the reference circuit
                 measures as record j (the decoder reads clbit j as record j)
  3. noiseless   sampled without noise and converted to detection events
                 through the same path as the hardware data
                 (measurements_to_detections): no detector or observable
                 ever fires
  4. errors      every single X, Y and Z error on every code qubit at the
                 start and after every syndrome round produces exactly the
                 same detection events and logical flip in both circuits

The rebuilt transpiled depth is compared with the depth saved when the job
was submitted. The serialized payloads are written to results/hw/payloads/.

Where this token can still see the job, the original payload is downloaded
from IQM's server instead (jobs/{id}/payload) and checked, and it must also
be identical to the rebuilt one.

Needs a token (read-only: fetches the backend description and job payloads,
submits nothing).

    .venv/Scripts/python experiments/verify_circuits.py
"""

import glob
import json
import sys
from pathlib import Path

import numpy as np
import stim

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from qec import build_memory_circuit  # noqa: E402
from qec.decode import measurements_to_detections  # noqa: E402
from qec.hardware import (coupling_edges, get_backend, load_run,  # noqa: E402
                          payload_dict, stim_to_qiskit, transpile_pinned)

ROOT = Path(__file__).resolve().parent.parent
HW = ROOT / "results" / "hw"
OUT = HW / "payloads"


# ------------------------------------------------- IQM payload -> Stim
_CLIFFORDS = {}


def _clifford_gates():
    """Map each single-qubit Clifford (as a canonical tableau string) to a
    short Stim gate sequence."""
    if not _CLIFFORDS:
        names = ["I", "X", "Y", "Z", "H", "S", "S_DAG", "SQRT_X", "SQRT_X_DAG",
                 "SQRT_Y", "SQRT_Y_DAG", "H_XY", "H_YZ", "C_XYZ", "C_ZYX"]
        for a in names:
            for b in names:
                c = stim.Circuit()
                for g in (a, b):
                    if g != "I":
                        c.append(g, [0])
                key = str(stim.Tableau.from_circuit(c))
                _CLIFFORDS.setdefault(key, [g for g in (a, b) if g != "I"])
    return _CLIFFORDS


def prx_to_stim(angle: float, phase: float) -> list[str]:
    """PRX(theta, phi) = exp(-i theta/2 (cos phi X + sin phi Y)) as Stim gates.
    Raises if it is not a Clifford."""
    for v in (angle, phase):
        k = v / (np.pi / 2)
        if abs(k - round(k)) > 1e-6:
            raise ValueError(f"non-Clifford prx angle={angle} phase={phase}")
    n = np.cos(phase) * np.array([[0, 1], [1, 0]]) + np.sin(phase) * np.array([[0, -1j], [1j, 0]])
    u = np.cos(angle / 2) * np.eye(2) - 1j * np.sin(angle / 2) * n
    u = np.round(u, 12)
    tab = stim.Tableau.from_unitary_matrix(u, endian="little")
    return _clifford_gates()[str(tab)]


def payload_to_stim(instructions, name_to_index: dict[str, int]):
    """Return (stim circuit, clbit of each measurement in order, list of
    measurement counts at each barrier)."""
    c = stim.Circuit()
    clbits, barrier_marks = [], []
    for ins in instructions:
        qs = [name_to_index[q] for q in ins.locus]
        if ins.name == "prx":
            for g in prx_to_stim(ins.args["angle"], ins.args["phase"]):
                c.append(g, qs)
        elif ins.name == "cz":
            c.append("CZ", qs)
        elif ins.name == "measure":
            c.append("M", qs)
            clbits.append(int(ins.args["key"].rsplit("_", 1)[1]))
            c.append("TICK")
        elif ins.name == "barrier":
            barrier_marks.append(len(clbits))
            c.append("TICK")
        else:
            raise ValueError(f"unexpected instruction {ins.name}")
    return c, clbits, barrier_marks


def insert_after_measurements(c: stim.Circuit, n_meas: int, gate: str, q: int) -> stim.Circuit:
    """Copy of `c` with `gate` on qubit q right after the n_meas-th
    measurement and before the next non-measurement operation (n_meas = 0:
    at the very start, after any initial reset)."""
    out, seen, done = stim.Circuit(), 0, False
    for inst in c.flattened():
        if not done and seen == n_meas and inst.name not in ("M", "TICK", "R", "QUBIT_COORDS"):
            out.append(gate, [q])
            done = True
        out.append(inst)
        if inst.name == "M":
            seen += len(inst.targets_copy())
    if not done:
        raise ValueError(f"no insertion point after {n_meas} measurements")
    return out


def detections(c_unannotated: stim.Circuit, clbits, ref: stim.Circuit, shots: int, seed: int):
    """Sample `c_unannotated` (measurements in payload order) noiselessly,
    put them in clbit order and convert with the reference circuit's
    detectors, exactly like the hardware data."""
    m = c_unannotated.compile_sampler(seed=seed).sample(shots)
    arr = np.zeros_like(m)
    arr[:, clbits] = m
    return measurements_to_detections(arr, ref)


def reference_detections(c: stim.Circuit, ref: stim.Circuit, shots: int, seed: int):
    """Detection events of `c` (the reference with an injected error) through
    the same path, relative to the error-free reference `ref`."""
    m = c.compile_sampler(seed=seed).sample(shots)
    return measurements_to_detections(m, ref)


# ------------------------------------------------------------- per run
def verify_circuit(ref, tqc, backend, layout, edges, saved_depth, payload=None):
    """Check one circuit. `payload` is the IQM circuit to check (downloaded
    from the server); if None, `tqc` is serialized with IQM's serializer."""
    report = {"depth": tqc.depth(), "saved_depth": saved_depth}
    if payload is None:
        payload = backend.serialize_circuit(tqc)
    name_to_index = {backend.index_to_qubit_name(i): i for i in range(backend.num_qubits)}
    ops = [ins.name for ins in payload.instructions]

    # 1. structure
    edge_set = {frozenset(e) for e in edges}
    used = {name_to_index[q] for ins in payload.instructions if ins.name != "barrier" for q in ins.locus}
    cz = [frozenset(name_to_index[q] for q in ins.locus) for ins in payload.instructions if ins.name == "cz"]
    report["ops"] = {o: ops.count(o) for o in sorted(set(ops))}
    report["structure"] = (set(ops) <= {"prx", "cz", "measure", "barrier"}
                           and all(e in edge_set for e in cz)
                           and used == set(layout.values()))
    pc, clbits, marks = payload_to_stim(payload.instructions, name_to_index)  # raises if non-Clifford

    # 2. readout mapping: clbit j measures layout[qubit of reference record j]
    ref_meas_qubits = [t.value for inst in ref.flattened() if inst.name == "M"
                       for t in inst.targets_copy()]
    pay_meas_qubits = [name_to_index[ins.locus[0]] for ins in payload.instructions if ins.name == "measure"]
    by_clbit = dict(zip(clbits, pay_meas_qubits))
    report["readout"] = (sorted(clbits) == list(range(ref.num_measurements))
                         and all(by_clbit[j] == layout[q] for j, q in enumerate(ref_meas_qubits)))

    # 3. noiseless
    det, obs = detections(pc, clbits, ref, 256, seed=1)
    report["noiseless"] = not det.any() and not obs.any()

    # 4. single Pauli errors at the start and after every syndrome round
    #    (in the reference: after the same number of measurements)
    boundaries = [0] + [m for m in marks if m < ref.num_measurements]
    n_checked = n_bad = n_fire = 0
    for nm in boundaries:
        for q, phys in layout.items():
            for g in ("X", "Y", "Z"):
                d_r, o_r = reference_detections(insert_after_measurements(ref, nm, g, q), ref, 4, seed=2)
                d_p, o_p = detections(insert_after_measurements(pc, nm, g, phys), clbits, ref, 4, seed=2)
                n_checked += 1
                # Detection events of a single Pauli are deterministic (random
                # first-round outcomes cancel in the detectors).
                n_fire += int(d_r[0].sum() + o_r[0].sum() > 0)
                if not ((d_r == d_r[0]).all() and (d_p == d_r[0]).all() and (o_p == o_r[0]).all()
                        and (o_r == o_r[0]).all()):
                    n_bad += 1
    report["errors_checked"] = n_checked
    report["errors_mismatch"] = n_bad
    report["errors_detected"] = n_fire
    report["ok"] = (report["structure"] and report["readout"] and report["noiseless"]
                    and n_bad == 0 and (saved_depth is None or report["depth"] == saved_depth))
    return report, payload


def server_payload(qc: str, job_id: str):
    """The job's original payload from IQM's server, or None if this token
    cannot see the job."""
    from iqm.iqm_client import IQMClient

    from qec.hardware import SERVER_URL
    try:
        cl = IQMClient(SERVER_URL, quantum_computer=qc)
        req = cl._get_submit_circuits_payload(job_id)
        job = cl._iqm_server_client.get_job(job_id)
    except Exception:
        return None
    return {"circuits": list(req.circuits), "dd_mode": str(req.dd_mode),
            "calibration_set_id": str(job.compilation.calibration_set_id) if job.compilation else None}


def _same(a, b) -> bool:
    return [(i.name, tuple(i.locus), i.args) for i in a.instructions] ==            [(i.name, tuple(i.locus), i.args) for i in b.instructions]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    backends, edges = {}, {}
    all_ok = True
    for f in sorted(glob.glob(str(HW / "2026*.npz"))):
        run = load_run(Path(f))
        if run.get("fake"):
            continue
        qc = run["quantum_computer"]
        if qc not in backends:
            backends[qc] = get_backend(qc)
            edges[qc] = coupling_edges(backends[qc])
        be = backends[qc]
        layout = {int(k): v for k, v in run["layout"].items()}
        depths = run.get("depths") or [None] * len(run["rounds"])
        srv = server_payload(qc, run["job_id"])
        print(Path(f).name, "| payload: " + (
            f"downloaded from IQM (dd_mode {srv['dd_mode']}, calibration set {srv['calibration_set_id']})"
            if srv else "rebuilt (job not visible to this token)"), flush=True)
        for k, (r, saved) in enumerate(zip(run["rounds"], depths)):
            ref = build_memory_circuit(3, int(r), run["basis"], True, None)
            tqc = transpile_pinned(stim_to_qiskit(ref), be, layout)
            downloaded = srv["circuits"][k] if srv else None
            rep, payload = verify_circuit(ref, tqc, be, layout, edges[qc], saved, downloaded)
            if srv:
                rep["same_as_rebuilt"] = _same(downloaded, be.serialize_circuit(tqc))
                rep["ok"] &= rep["same_as_rebuilt"]
            all_ok &= rep["ok"]
            (OUT / f"{Path(f).stem}_r{r}.json").write_text(json.dumps(payload_dict(be, tqc)))
            print(f"  rounds={r}: {'OK ' if rep['ok'] else 'FAIL'} depth {rep['depth']} (saved {rep['saved_depth']}), "
                  f"structure={rep['structure']} readout={rep['readout']} noiseless={rep['noiseless']} "
                  f"errors {rep['errors_checked'] - rep['errors_mismatch']}/{rep['errors_checked']} match "
                  f"({rep['errors_detected']} visible)"
                  f"{'  identical to rebuilt=' + str(rep['same_as_rebuilt']) if srv else ''}  "
                  f"{rep['ops']}", flush=True)
    print("\nALL CIRCUITS VERIFIED" if all_ok else "\nSOME CHECKS FAILED")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
