"""Decoding detection events into logical error probabilities.

All decoders consume the same inputs: a (shots, num_detectors) bool array of
detection events, the (shots, num_observables) observed logical flips, and a
detector error model (DEM) built from a *noisy* version of exactly the circuit
that produced the data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pymatching
import stim


@dataclass
class DecodeResult:
    decoder: str
    shots: int
    errors: int

    @property
    def p_logical(self) -> float:
        return self.errors / self.shots

    @property
    def stderr(self) -> float:
        p = self.p_logical
        return float(np.sqrt(max(p * (1 - p), 1e-12) / self.shots))

    def wilson_interval(self, z: float = 1.0) -> tuple[float, float]:
        """Wilson score interval (z=1 -> ~68%), well-behaved near p=0."""
        n, p = self.shots, self.p_logical
        denom = 1 + z * z / n
        centre = (p + z * z / (2 * n)) / denom
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
        return float(centre - half), float(centre + half)

    def __repr__(self):
        return (f"{self.decoder}: P_L = {self.p_logical:.4f} +/- {self.stderr:.4f} "
                f"({self.errors}/{self.shots})")


def make_matcher(dem: stim.DetectorErrorModel, correlated: bool = False) -> pymatching.Matching:
    if correlated:
        return pymatching.Matching.from_detector_error_model(dem, enable_correlations=True)
    return pymatching.Matching.from_detector_error_model(dem)


def decode(
    det_events: np.ndarray,
    obs_flips: np.ndarray,
    noisy_circuit: stim.Circuit,
    decoder: str = "mwpm",
) -> DecodeResult:
    """Decode and count logical failures.

    decoder : "mwpm"       - PyMatching minimum-weight perfect matching
              "correlated" - PyMatching with correlated (X/Y-aware) matching
              "union_find" - weighted Union-Find (peeling) on the matching graph
              "bposd"      - belief propagation + ordered-statistics decoding on
                             the full (hyperedge) error model
              "belief_matching" - BP to reweight the graph, then MWPM
              "none"       - no correction: raw logical flip rate (a baseline
                             that shows how much the decoder actually helps)
    """
    if np.asarray(det_events).shape[1] != noisy_circuit.num_detectors:
        raise ValueError(
            f"data has {np.asarray(det_events).shape[1]} detectors but the DEM circuit has "
            f"{noisy_circuit.num_detectors}: the DEM must come from the circuit that ran")
    return decode_dem(det_events, obs_flips,
                      noisy_circuit.detector_error_model(decompose_errors=True), decoder)


def decode_dem(det_events: np.ndarray, obs_flips: np.ndarray,
               dem: stim.DetectorErrorModel, decoder: str = "mwpm") -> DecodeResult:
    """Like `decode`, from a (decomposed) detector error model directly, e.g.
    one whose probabilities were learned from data (qec.learn.reweighted_dem)."""
    det_events = np.asarray(det_events, dtype=bool)
    obs_flips = np.asarray(obs_flips, dtype=bool).reshape(len(det_events), -1)
    if det_events.shape[1] != dem.num_detectors:
        raise ValueError(f"data has {det_events.shape[1]} detectors, DEM has {dem.num_detectors}")

    if decoder == "none":
        pred = np.zeros_like(obs_flips)
    elif decoder in ("union_find", "bposd", "belief_matching"):
        pred = _decode_ldpc(det_events, dem, decoder)
    else:
        correlated = decoder == "correlated"
        matcher = make_matcher(dem, correlated=correlated)
        # Correlated matching must be requested both when building the graph
        # and when decoding (two-pass reweighting happens at decode time).
        pred = matcher.decode_batch(det_events, enable_correlations=correlated)
    errors = int(np.any(pred != obs_flips, axis=1).sum())
    return DecodeResult(decoder, len(det_events), errors)


def measurements_to_detections(raw_meas: np.ndarray, circuit: stim.Circuit):
    """Raw measurement records -> (det_events, obs_flips) using the circuit's
    own DETECTOR / OBSERVABLE_INCLUDE definitions (no-reset aware)."""
    conv = circuit.compile_m2d_converter()
    return conv.convert(measurements=np.asarray(raw_meas, dtype=bool),
                        separate_observables=True)


def syndrome_diagnostics(det_events: np.ndarray) -> dict:
    det_events = np.asarray(det_events, dtype=bool)
    w = det_events.sum(axis=1)
    return {
        "mean_syndrome_weight": float(w.mean()),
        "detection_fraction": float(det_events.mean()),
        "error_free_shots": float((w == 0).mean()),
        "detector_firing_rate": det_events.mean(axis=0),
    }


def _decode_ldpc(det_events: np.ndarray, dem: stim.DetectorErrorModel, decoder: str) -> np.ndarray:
    """Decoders from the `ldpc` / `beliefmatching` packages (pip install ldpc beliefmatching)."""
    from beliefmatching import BeliefMatching, detector_error_model_to_check_matrices

    if decoder == "belief_matching":
        return BeliefMatching.from_detector_error_model(dem, max_bp_iters=20).decode_batch(det_events)

    import ldpc

    mats = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)
    if decoder == "bposd":
        H, L, priors = mats.check_matrix, mats.observables_matrix, mats.priors
        dec = ldpc.BpOsdDecoder(H, error_channel=list(priors), max_iter=30,
                                bp_method="product_sum", osd_method="osd_cs", osd_order=7)
        decode_one = dec.decode
    else:
        # Union-Find runs on the graphlike decomposition: each edge's
        # probability is the XOR-combination of the hyperedges that contain it.
        H, L = mats.edge_check_matrix, mats.edge_observables_matrix
        h2e = mats.hyperedge_to_edge_matrix.tocsr()
        q = np.ones(H.shape[1])
        for e in range(H.shape[1]):
            for h in h2e.getrow(e).indices:
                q[e] *= 1 - 2 * mats.priors[h]
        p_edge = np.clip(0.5 * (1 - q), 1e-9, 0.5 - 1e-9)
        llrs = np.log((1 - p_edge) / p_edge)
        dec = ldpc.UnionFindDecoder(H, uf_method="peeling")
        decode_one = lambda s: dec.decode(s, llrs=llrs)  # noqa: E731

    L = L.tocsr()
    out = np.zeros((len(det_events), L.shape[0]), dtype=bool)
    for k, syndrome in enumerate(det_events.astype(np.uint8)):
        err = decode_one(syndrome)
        out[k] = (L @ err) % 2
    return out
