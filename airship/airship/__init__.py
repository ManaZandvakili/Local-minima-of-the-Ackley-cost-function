"""Optimized propeller layouts and closed-loop control for a spherical airship.

Everything in this package works in physical units:
    positions  [m]   body frame, origin at the hull centre (centre of buoyancy),
                     z up, gondola at -z
    directions [-]   unit thrust vectors
    thrust     [N]
    torque     [N*m]

The conference-paper code (reference/Final.py) is kept unchanged so its
results stay reproducible; functions reused from it say so in their docstrings.
"""
