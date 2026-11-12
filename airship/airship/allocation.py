"""Best-effort bounded allocation: wrench command -> per-propeller thrusts.

    min_u  || W (M u - tau_cmd) ||^2 + eps ||u||^2     s.t.  u_min <= u <= u_max

If tau_cmd is achievable the result reproduces it (up to the tiny eps
regularisation, which picks the smallest thrusts among exact solutions). If
not, it returns the closest achievable wrench, weighting errors by W (default
priorities: Fx, Fy and yaw first). It never fails: a failed propeller is
simply u_min = u_max = 0.

Solved exactly with a primal active-set method for the box-constrained
quadratic program (numba-compiled, warm-started from the previous solution),
because scipy's lsq_linear costs ~0.6 ms per call, too slow for training.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

DEFAULT_W = (2.0, 2.0, 1.0, 0.5, 0.5, 2.0)      # Fx Fy Fz Mx My Mz (brief)


@njit(cache=True)
def box_qp(H, b, lo, hi, u0, max_iter=200):
    """min 0.5 u'Hu - b'u  s.t.  lo <= u <= hi, H positive definite.

    Primal active-set method started from u0 (clipped into the box). Returns u.
    """
    n = b.shape[0]
    u = np.minimum(np.maximum(u0, lo), hi)
    # status: 0 free, -1 at lower bound, +1 at upper bound
    st = np.zeros(n, dtype=np.int64)
    for i in range(n):
        if u[i] <= lo[i]:
            st[i] = -1
            u[i] = lo[i]
        elif u[i] >= hi[i]:
            st[i] = 1
            u[i] = hi[i]
    for _ in range(max_iter):
        free = np.where(st == 0)[0]
        # minimiser over the free variables with bounded ones fixed
        target = u.copy()
        if free.shape[0] > 0:
            rhs = b[free] - H[free][:, st != 0] @ u[st != 0]
            target[free] = np.linalg.solve(H[free][:, free], rhs)
        # step toward the target until a free variable hits a bound
        alpha, block, block_st = 1.0, -1, 0
        for i in free:
            d = target[i] - u[i]
            if d < 0 and target[i] < lo[i]:
                a = (lo[i] - u[i]) / d
                if a < alpha:
                    alpha, block, block_st = a, i, -1
            elif d > 0 and target[i] > hi[i]:
                a = (hi[i] - u[i]) / d
                if a < alpha:
                    alpha, block, block_st = a, i, 1
        for i in free:
            u[i] += alpha * (target[i] - u[i])
        if block >= 0:
            st[block] = block_st
            u[block] = lo[block] if block_st < 0 else hi[block]
            continue
        # at the subspace minimum: release the bound with the worst multiplier
        g = H @ u - b
        worst, wi = 1e-12, -1
        for i in range(n):
            if st[i] == -1 and -g[i] > worst:       # pushing up from the lower bound
                worst, wi = -g[i], i
            elif st[i] == 1 and g[i] > worst:       # pushing down from the upper bound
                worst, wi = g[i], i
        if wi < 0:
            return u
        st[wi] = 0
    return u


@dataclass
class AllocationResult:
    u: np.ndarray           # per-propeller thrust [N]
    achieved: np.ndarray    # M u
    gap: np.ndarray         # tau_cmd - M u
    saturated: np.ndarray   # bool per propeller: at a bound (and the bound is not 0 for push-only idle)


class Allocator:
    def __init__(self, M, u_min, u_max, W=DEFAULT_W, eps=1e-8):
        self.M = np.ascontiguousarray(M, dtype=float)
        self.W2 = np.asarray(W, dtype=float) ** 2
        self.eps = float(eps)
        n = self.M.shape[1]
        self.set_bounds(u_min, u_max)
        self.H = self.M.T @ (self.W2[:, None] * self.M) + self.eps * np.eye(n)
        self.u_prev = np.zeros(n)

    def set_bounds(self, u_min, u_max):
        n = self.M.shape[1]
        self.lo = np.broadcast_to(np.asarray(u_min, dtype=float), (n,)).copy()
        self.hi = np.broadcast_to(np.asarray(u_max, dtype=float), (n,)).copy()

    def reset(self):
        self.u_prev[:] = 0.0

    def solve(self, tau_cmd):
        tau = np.asarray(tau_cmd, dtype=float)
        b = self.M.T @ (self.W2 * tau)
        u = box_qp(self.H, b, self.lo, self.hi, self.u_prev)
        self.u_prev = u
        achieved = self.M @ u
        span = self.hi - self.lo
        sat = (span > 0) & ((u >= self.hi - 1e-9) | ((u <= self.lo + 1e-9) & (self.lo < 0)))
        return AllocationResult(u, achieved, tau - achieved, sat)
