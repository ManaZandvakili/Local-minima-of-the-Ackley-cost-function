"""Wind: mean wind + Ornstein-Uhlenbeck turbulence + discrete 1-cosine gusts.

World frame (z up). Everything is driven by one seeded generator, so every
controller and layout flown with the same seed sees exactly the same wind.
The OU process is a first-order approximation of Dryden turbulence: each axis
has standard deviation sigma and correlation time tau.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class WindConfig:
    mean_speed_m_s: float = 0.0
    mean_dir_deg: float = 0.0                 # direction the wind blows toward, from +x, CCW
    turbulence_sigma_m_s: tuple = (0.5, 0.5, 0.25)
    turbulence_tau_s: float = 2.0
    gust_rate_per_min: float = 0.0
    gust_amplitude_m_s: tuple = (1.0, 3.0)
    gust_duration_s: tuple = (2.0, 6.0)

    @classmethod
    def from_config(cls, d):
        d = dict(d or {})
        for k in ("turbulence_sigma_m_s", "gust_amplitude_m_s", "gust_duration_s"):
            if k in d:
                d[k] = tuple(d[k])
        return cls(**d)


class WindModel:
    def __init__(self, cfg, seed, duration_s=600.0):
        self.cfg = cfg
        self.rng = np.random.default_rng(seed)
        a = np.deg2rad(cfg.mean_dir_deg)
        self.mean = cfg.mean_speed_m_s * np.array([np.cos(a), np.sin(a), 0.0])
        self.sigma = np.asarray(cfg.turbulence_sigma_m_s, dtype=float)
        self.turb = self.rng.normal(size=3) * self.sigma          # start in the stationary distribution
        self.t = 0.0
        self.gusts = self._schedule(duration_s)

    def _schedule(self, duration_s):
        """Gust list (t0, duration, peak vector), drawn up front from the seed."""
        c, out, t = self.cfg, [], 0.0
        if c.gust_rate_per_min <= 0:
            return out
        while True:
            t += self.rng.exponential(60.0 / c.gust_rate_per_min)
            if t >= duration_s:
                return out
            dur = self.rng.uniform(*c.gust_duration_s)
            amp = self.rng.uniform(*c.gust_amplitude_m_s)
            ang = self.rng.uniform(0, 2 * np.pi)
            out.append((t, dur, amp * np.array([np.cos(ang), np.sin(ang), 0.0])))

    def gust(self, t):
        g = np.zeros(3)
        for t0, dur, peak in self.gusts:
            if t0 <= t < t0 + dur:
                g += 0.5 * peak * (1.0 - np.cos(2 * np.pi * (t - t0) / dur))
        return g

    def sample(self, n, dt):
        """Wind vectors at the start of the next n steps of length dt, shape (n, 3)."""
        out = np.empty((n, 3))
        tau = self.cfg.turbulence_tau_s
        a = np.exp(-dt / tau)
        b = self.sigma * np.sqrt(1.0 - a * a)
        for k in range(n):
            out[k] = self.mean + self.turb + self.gust(self.t)
            self.turb = a * self.turb + b * self.rng.normal(size=3)   # exact OU update
            self.t += dt
        return out
