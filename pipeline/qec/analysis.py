"""Logical error rate per round from a sweep over round counts.

Model (Google Quantum AI, Nature 2021 / 2023): the logical fidelity
F(n) = 1 - 2 P_L(n) decays exponentially with the number of rounds n,

    F(n) = A (1 - 2 eps_L)^n,

where eps_L is the logical error rate per round. The challenge README uses
A = 1; letting A float absorbs state-preparation and final-readout error so
that eps_L reflects only the per-round cost. We report both.

Fitting is a binomial maximum likelihood over all (n, errors, shots) points,
with the 1-sigma uncertainty from the observed Fisher information (inverse
Hessian of the negative log-likelihood).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize


def p_logical_model(n, eps, amp=1.0):
    n = np.asarray(n, dtype=float)
    return 0.5 * (1 - amp * (1 - 2 * eps) ** n)


def eps_from_single(p_l: float, n: int) -> float:
    """Per-round rate from one round count, A = 1 (README closed form)."""
    f = 1 - 2 * p_l
    if f <= 0:
        return 0.5
    return 0.5 * (1 - f ** (1 / n))


@dataclass
class LERFit:
    eps: float
    eps_err: float
    amp: float
    amp_err: float
    free_amplitude: bool
    rounds: np.ndarray
    p_obs: np.ndarray

    def __repr__(self):
        s = f"eps_L = {100*self.eps:.3f} +/- {100*self.eps_err:.3f} % per round"
        if self.free_amplitude:
            s += f"  (A = {self.amp:.3f} +/- {self.amp_err:.3f})"
        return s


def fit_ler(rounds, errors, shots, free_amplitude: bool = False,
            min_round: int = 1) -> LERFit:
    """Binomial-MLE fit of P_L(n) = (1 - A (1-2 eps)^n) / 2.

    min_round : drop round counts below this (Google excludes the first
                cycle(s), which carry extra boundary/state-prep error).
    """
    rounds = np.asarray(rounds, dtype=float)
    errors = np.asarray(errors, dtype=float)
    shots = np.asarray(shots, dtype=float)
    keep = rounds >= min_round
    r, k, N = rounds[keep], errors[keep], shots[keep]

    def unpack(theta):
        eps = 0.5 / (1 + np.exp(-theta[0]))          # eps in (0, 0.5)
        amp = 1 / (1 + np.exp(-theta[1])) if free_amplitude else 1.0
        return eps, amp

    def nll(theta):
        eps, amp = unpack(theta)
        p = np.clip(p_logical_model(r, eps, amp), 1e-12, 1 - 1e-12)
        return -np.sum(k * np.log(p) + (N - k) * np.log(1 - p))

    # Initial guess from the pooled single-round estimates.
    eps0 = np.clip(np.median([eps_from_single(ki / Ni, int(ri))
                              for ri, ki, Ni in zip(r, k, N)]), 1e-5, 0.45)
    th0 = [np.log(2 * eps0 / (1 - 2 * eps0))]
    if free_amplitude:
        th0.append(2.0)
    res = minimize(nll, th0, method="Nelder-Mead",
                   options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 20000})
    eps, amp = unpack(res.x)

    # Uncertainty: numerical Hessian in the natural parameters (eps, amp).
    def nll_nat(v):
        e = v[0]
        a = v[1] if free_amplitude else 1.0
        p = np.clip(p_logical_model(r, e, a), 1e-12, 1 - 1e-12)
        return -np.sum(k * np.log(p) + (N - k) * np.log(1 - p))

    x0 = np.array([eps, amp]) if free_amplitude else np.array([eps])
    H = _hessian(nll_nat, x0)
    try:
        cov = np.linalg.inv(H)
        errs = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        errs = np.full(len(x0), np.nan)
    return LERFit(eps=float(eps), eps_err=float(errs[0]), amp=float(amp),
                  amp_err=float(errs[1]) if free_amplitude else 0.0,
                  free_amplitude=free_amplitude, rounds=rounds,
                  p_obs=errors / shots)


def _hessian(f, x, rel=1e-4):
    x = np.asarray(x, dtype=float)
    n = len(x)
    h = np.maximum(np.abs(x) * rel, 1e-7)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei = np.zeros(n); ei[i] = h[i]
            ej = np.zeros(n); ej[j] = h[j]
            H[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej)
                       + f(x - ei - ej)) / (4 * h[i] * h[j])
    return H


def plot_ler(ax, rounds, p_obs, p_err, fit: LERFit | None = None, label="", color=None):
    """P_L vs rounds with optional fitted curve (Fig. 3 style, 2023 paper)."""
    ax.errorbar(rounds, p_obs, yerr=p_err, fmt="o", color=color, label=label,
                capsize=3, ms=5)
    if fit is not None:
        n = np.linspace(min(rounds), max(rounds), 200)
        ax.plot(n, p_logical_model(n, fit.eps, fit.amp), "-", color=color, alpha=0.8,
                label=f"fit: $\\varepsilon_L$ = {100*fit.eps:.2f}%")
    ax.set_xlabel("rounds $n$")
    ax.set_ylabel("logical error probability $P_L$")
    ax.axhline(0.5, color="grey", lw=0.8, ls=":")
