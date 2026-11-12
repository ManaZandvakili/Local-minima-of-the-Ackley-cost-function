"""3-D animation of an episode (brief, section 10) -> MP4 (ffmpeg) or GIF.

Shows the translucent hull, the gondola and its keep-out cone, every propeller
with a thrust arrow (length proportional to thrust; grey -> orange -> red as it
nears its limit; failed motors as a black X), the trajectory trail, the target,
a wind arrow, the time, and position/yaw error and sum of squared thrusts.
Side-by-side mode puts two runs (same wind seed) next to each other.

Input: the `log` dict returned by closed_loop.run_episode(..., log=True).
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

from ..dynamics import euler_from_quat, quat_to_rot

GREY, ORANGE, RED = "#9aa3ab", "#e07b00", "#d62728"


def _sphere(R, n=18):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, n), np.linspace(0, np.pi, n // 2 + 1))
    return R * np.cos(u) * np.sin(v), R * np.sin(u) * np.sin(v), R * np.cos(v)


class _Panel:
    def __init__(self, ax, log, layout, title, extent, gondola_drop=0.25, trail_s=20.0, dt=0.1):
        self.ax, self.log, self.L, self.dt = ax, log, layout, dt
        self.R = layout.semi_axes_m[0]
        self.extent, self.trail = extent, int(trail_s / dt)
        self.drop = gondola_drop
        self.umax = np.maximum(np.abs(layout.u_max), np.abs(layout.u_min))
        self.title = title
        self.fail_step = np.full(layout.n_props, -1)
        if "fail_prop" in log:
            self.fail_step[int(log["fail_prop"])] = int(log["fail_step"])

    def draw(self, k):
        ax, L, x = self.ax, self.L, self.log["x"][k]
        ax.cla()
        e = self.extent
        ax.set_xlim(-e, e); ax.set_ylim(-e, e); ax.set_zlim(-e, e)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
        p = x[:3]
        Rm = quat_to_rot(x[3:7])
        hx, hy, hz = _sphere(self.R)
        pts = Rm @ np.stack([hx.ravel(), hy.ravel(), hz.ravel()])
        ax.plot_wireframe((pts[0] + p[0]).reshape(hx.shape), (pts[1] + p[1]).reshape(hx.shape),
                          (pts[2] + p[2]).reshape(hx.shape), color="#9fc5e8", lw=0.4, alpha=0.5)
        # gondola, suspension line and a faint keep-out cone
        g = Rm @ np.array([0, 0, -self.R - self.drop]) + p
        b = Rm @ np.array([0, 0, -self.R]) + p
        ax.plot(*zip(b, g), color="#555555", lw=1)
        ax.scatter(*g, color="#555555", s=40, marker="s")
        k_r = self.R * np.sin(np.deg2rad(35))
        th = np.linspace(0, 2 * np.pi, 24)
        ring = Rm @ np.stack([k_r * np.cos(th), k_r * np.sin(th), np.full_like(th, -self.R * np.cos(np.deg2rad(35)))])
        ax.plot(ring[0] + p[0], ring[1] + p[1], ring[2] + p[2], color=RED, lw=0.6, alpha=0.35)
        # propellers and thrust arrows
        u = self.log["u"][k]
        P = (Rm @ L.positions_m.T).T + p
        D = (Rm @ L.directions.T).T
        for i in range(L.n_props):
            if self.fail_step[i] >= 0 and k >= self.fail_step[i]:
                ax.scatter(*P[i], color="black", marker="x", s=60, linewidths=2)
                continue
            frac = abs(u[i]) / max(self.umax[i], 1e-9)
            col = RED if frac > 0.95 else ORANGE if frac > 0.5 else GREY
            ax.scatter(*P[i], color=col, s=14)
            vec = D[i] * np.sign(u[i]) * 1.2 * frac
            if frac > 0.02:
                ax.quiver(*P[i], *vec, color=col, lw=1.4, arrow_length_ratio=0.25)
        # trail, target, wind
        tr = self.log["x"][max(0, k - self.trail):k + 1, :3]
        ax.plot(tr[:, 0], tr[:, 1], tr[:, 2], color="#1f4e9c", lw=1.2)
        ax.scatter(0, 0, 0, color="green", marker="+", s=80)
        w = self.log["wind"][k]
        c = np.array([-0.85 * e, -0.85 * e, 0.85 * e])
        ax.quiver(*c, *(w * 0.8), color="#2a9d8f", lw=2, arrow_length_ratio=0.2)
        roll, pitch, yaw = euler_from_quat(x[3:7])
        err = np.linalg.norm(p)
        yerr = np.degrees(self.log.get("yaw_ref", np.zeros(len(self.log["x"])))[k] - yaw)
        yerr = (yerr + 180) % 360 - 180
        ax.set_title(f"{self.title}\nt = {self.log['t'][k]:5.1f} s   |pos err| = {err:4.2f} m   "
                     f"yaw err = {yerr:5.1f}°   Σu² = {np.sum(u ** 2):5.0f} N²   wind {np.linalg.norm(w):3.1f} m/s",
                     fontsize=9)


def animate(runs, path, extent=4.0, stride=2, fps=10, dpi=90):
    """runs: list of (log, layout, title); one or two panels. Writes MP4 (or GIF by suffix)."""
    n = len(runs)
    fig = plt.figure(figsize=(6.4 * n, 6.2))
    panels = [_Panel(fig.add_subplot(1, n, i + 1, projection="3d"), log, lay, title, extent)
              for i, (log, lay, title) in enumerate(runs)]
    frames = range(0, min(len(r[0]["x"]) for r in runs), stride)

    def update(k):
        for pn in panels:
            pn.draw(k)
        return []

    anim = animation.FuncAnimation(fig, update, frames=frames, blit=False)
    path = str(path)
    if path.endswith(".gif"):
        anim.save(path, writer=animation.PillowWriter(fps=fps), dpi=dpi)
    else:
        anim.save(path, writer=animation.FFMpegWriter(fps=fps, bitrate=1800), dpi=dpi)
    plt.close(fig)
    return path
