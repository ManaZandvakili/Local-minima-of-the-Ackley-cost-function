"""Controllers that output a 6-D wrench command [Fx, Fy, Fz, Mx, My, Mz] (body frame).

The wrench goes through the same allocation (airship.allocation) for every
controller, so PID and RL are compared fairly.
"""
