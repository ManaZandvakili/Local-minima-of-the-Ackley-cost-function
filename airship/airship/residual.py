"""Residual RL support: the PID gains each training layout flies with.

In residual mode the RL action is a correction added to the PID wrench:
    wrench = PID(x) + residual_scale * action
The PID part uses that layout's tuned gains (Phase 3). Perturbed and
one-motor-failure variants use their parent layout's gains; layouts that were
never tuned (most Pareto layouts) use the median of all tuned gains.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from .controllers.pid import GAIN_NAMES, PIDGains


class GainsBook:
    def __init__(self, pid_dir):
        self.gains = {}
        for f in sorted((Path(pid_dir) / "gains").glob("*.json")):
            g = json.load(open(f))
            self.gains[g["layout"]] = PIDGains(**g["gains"])
        vals = np.array([[getattr(g, k) for k in GAIN_NAMES] for g in self.gains.values()])
        self.median = PIDGains.from_vector(np.median(vals, axis=0))

    @staticmethod
    def label(layout_name):
        """Pool name "<run>:<layout>[_pertK|_minusI]" -> PID-run label."""
        run, _, name = layout_name.rpartition(":")
        name = re.sub(r"(_pert\d+|_minus\d+)+$", "", name)
        m = re.search(r"_(maxN\d+)_", run)
        return f"{m.group(1)}_{name}" if m else name

    def for_layout(self, layout_name):
        return self.gains.get(self.label(layout_name), self.median)
