"""Load YAML configs and build the objects they describe."""
from __future__ import annotations

from pathlib import Path

import yaml

from .geometry import Hull

CONFIG_DIR = Path(__file__).resolve().parent.parent / "configs"


def load_config(name_or_path="sphere"):
    p = Path(name_or_path)
    if not p.suffix:
        p = CONFIG_DIR / f"{name_or_path}.yaml"
    return yaml.safe_load(p.read_text())


def hull_from_config(cfg):
    h = cfg["hull"]
    if h["shape"] == "sphere":
        return Hull.sphere(h["diameter_m"] / 2.0)
    return Hull(h["shape"], tuple(h["semi_axes_m"]))
