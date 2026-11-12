"""Mass, centre of gravity and inertia of the airship from its mass budget.

All quantities in the body frame: origin at the hull centre (= centre of
buoyancy), z up. Components are idealised as thin spherical shells, solid
spheres, points on the z axis, or horizontal rings, so the products of inertia
vanish and the CG lies on the z axis.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _component_inertia(kg, kind, R, z=0.0, ring=0.0):
    """Mass, z of its CG, and diagonal inertia about the body origin."""
    if kind == "shell":
        i = 2.0 / 3.0 * kg * R**2
        return kg, 0.0, np.array([i, i, i])
    if kind == "solid_sphere":
        i = 2.0 / 5.0 * kg * R**2
        return kg, 0.0, np.array([i, i, i])
    if kind == "point":
        return kg, z, np.array([kg * z**2, kg * z**2, 0.0])
    if kind == "ring":
        ixx = kg * (z**2 + ring**2 / 2.0)
        return kg, z, np.array([ixx, ixx, kg * ring**2])
    raise ValueError(f"unknown component kind {kind!r}")


@dataclass(frozen=True)
class MassProperties:
    mass_kg: float              # total mass, including helium and trim
    buoyancy_kg: float          # rho * V (displaced-air mass)
    cg_m: np.ndarray            # centre of gravity [m], body frame
    inertia_kgm2: np.ndarray    # 3x3 inertia about the body origin (hull centre)
    added_mass_kg: float        # translational added mass, each axis
    added_inertia_kgm2: np.ndarray
    trim_kg: float
    breakdown: dict             # name -> (kg, z_m)

    @property
    def h_cg_m(self):
        """How far the CG sits below the hull centre [m]."""
        return -float(self.cg_m[2])


def from_config(cfg):
    """Build mass properties from configs/sphere.yaml."""
    R = cfg["hull"]["diameter_m"] / 2.0
    rho = cfg["air"]["density_kg_m3"]
    V = cfg["hull"]["volume_m3"]
    mcfg = cfg["mass"]
    buoy = rho * V
    total, mz, inertia, parts = 0.0, 0.0, np.zeros(3), {}
    for name, c in mcfg["components"].items():
        kind = c["kind"]
        if name == "helium" and not mcfg.get("helium_in_inertia", True):
            kg, z, i = c["kg"], 0.0, np.zeros(3)
        else:
            kg, z, i = _component_inertia(c["kg"], kind, R, c.get("z_m", 0.0), c.get("ring_m", 0.0))
        total += kg; mz += kg * z; inertia += i; parts[name] = (kg, z)
    trim = 0.0
    t = mcfg.get("trim", {})
    if t.get("mode", "neutral") == "neutral":
        trim = buoy - total
        if trim < 0:
            raise ValueError(f"components ({total:.2f} kg) are heavier than the buoyancy ({buoy:.2f} kg)")
        kg, z, i = _component_inertia(trim, t["kind"], R, t.get("z_m", 0.0), t.get("ring_m", 0.0))
        total += kg; mz += kg * z; inertia += i; parts["trim_ballast"] = (kg, z)
    return MassProperties(
        mass_kg=total, buoyancy_kg=buoy, cg_m=np.array([0.0, 0.0, mz / total]),
        inertia_kgm2=np.diag(inertia), added_mass_kg=mcfg.get("added_mass_coeff", 0.5) * rho * V,
        added_inertia_kgm2=np.diag(np.asarray(mcfg.get("added_inertia_kgm2", [0, 0, 0]), dtype=float)),
        trim_kg=trim, breakdown=parts,
    )
