"""
Slide 6 visual: one propeller on an airship hull showing
    r_i  - mounting location (from the body-frame origin / CG)
    d_i  - unit thrust direction
    F_i  - force  u_i d_i
    tau_i - resulting moment  r_i x u_i d_i  about the origin

Side view (X-Z plane).  The moment sense drawn around the origin is
computed from the actual r_i and d_i, so it is always consistent.

Usage:
    python thrust_allocation_visual.py
        -> saves thrust_allocation_visual.png (transparent) and .svg
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, Polygon, FancyBboxPatch

A, C          = 3.2, 1.35          # hull semi-axes in the drawing (X, Z)
PROP_ANGLE    = 52.0               # propeller location on hull [deg from +X]
THRUST_DIR    = 150.0              # thrust direction d_i [deg from +X]
FORCE_LEN     = 1.55               # drawn length of u_i d_i
UNIT_LEN      = 0.75               # drawn length of the unit vector d_i
TRANSPARENT   = True
OUT_BASENAME  = "thrust_allocation_visual"

C_HULL  = "#cfe6f2"; C_EDGE = "#6a8fa6"
C_R     = "#1f3b73"              # r_i
C_D     = "#111111"              # d_i
C_F     = "#e07b00"              # force
C_TAU   = "#c0392b"              # moment


def arrow(ax, p0, p1, color, lw, ms=22, z=5, ls="-"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms,
                                 color=color, lw=lw, zorder=z, linestyle=ls,
                                 shrinkA=0, shrinkB=0))


def main():
    fig, ax = plt.subplots(figsize=(10, 6.2))
    ax.set_aspect("equal"); ax.axis("off")

    # ── Hull, fins, gondola ─────────────────────────────────────────────
    ax.add_patch(Ellipse((0, 0), 2*A, 2*C, fc=C_HULL, ec=C_EDGE, lw=2, zorder=1))
    for s in (+1, -1):                                       # tail fins
        ax.add_patch(Polygon([(-A*0.78, s*C*0.55), (-A*1.08, s*C*1.05),
                              (-A*1.12, s*C*1.05), (-A*0.98, s*C*0.35)],
                             closed=True, fc="#9aa9b3", ec=C_EDGE, lw=1.2, zorder=0))
    ax.add_patch(FancyBboxPatch((-0.55, -C - 0.42), 1.1, 0.34,
                                boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc="#7b868e", ec="#5d666c", zorder=1))
    ax.plot([-0.35, -0.25], [-C - 0.08, -C + 0.02], color="#5d666c", lw=1.5)
    ax.plot([0.35, 0.25], [-C - 0.08, -C + 0.02], color="#5d666c", lw=1.5)

    # ── Body frame at origin (CG) ──────────────────────────────────────
    O = np.array([0.0, 0.0])
    arrow(ax, O, (0.75, 0), "#555555", 1.6, ms=14, z=3)
    arrow(ax, O, (0, 0.75), "#555555", 1.6, ms=14, z=3)
    ax.text(0.80, -0.05, r"$x$", fontsize=15, color="#555555", va="top")
    ax.text(-0.13, 0.80, r"$z$", fontsize=15, color="#555555")
    ax.plot(*O, "o", ms=8, color="black", zorder=6)
    ax.text(-0.15, 0.10, "CG", fontsize=14, ha="right", va="bottom")

    # ── Propeller location and thrust ──────────────────────────────────
    t  = np.deg2rad(PROP_ANGLE)
    P  = np.array([A*np.cos(t), C*np.sin(t)])
    d  = np.array([np.cos(np.deg2rad(THRUST_DIR)), np.sin(np.deg2rad(THRUST_DIR))])

    # line of action (faint) and moment arm
    L = np.array([P - 3.0*d, P + 3.0*d])
    ax.plot(L[:, 0], L[:, 1], ls=(0, (4, 4)), color=C_F, lw=1.0, alpha=0.45, zorder=2)
    foot = P + np.dot(O - P, d) * d                 # perpendicular foot from CG
    ax.plot([O[0], foot[0]], [O[1], foot[1]], ls=":", color="#777777", lw=1.4, zorder=2)

    arrow(ax, O, P - 0.30*P/np.linalg.norm(P), C_R, 3.0, ms=24, z=4)   # r_i (stops short of the prop)
    mid = 0.55*P
    ax.text(mid[0] + 0.12, mid[1] - 0.22, r"$\mathbf{r}_i$", fontsize=24,
            color=C_R, fontweight="bold")

    arrow(ax, P, P + FORCE_LEN*d, C_F, 7.0, ms=34, z=5)      # force u_i d_i
    tip = P + FORCE_LEN*d
    ax.text(tip[0] - 0.05, tip[1] + 0.12,
            r"$\mathbf{F}_i = u_i\,\mathbf{d}_i$", fontsize=22, color=C_F,
            ha="center", va="bottom")

    arrow(ax, P, P + UNIT_LEN*d, C_D, 2.2, ms=18, z=6)       # unit d_i on top
    n = np.array([-d[1], d[0]])
    lab = P + 0.5*UNIT_LEN*d - 0.30*n
    ax.text(lab[0], lab[1], r"$\mathbf{d}_i$", fontsize=22, color=C_D,
            ha="center", va="center")

    # propeller symbol: hub + disc perpendicular to d_i
    ax.plot([P[0] - 0.32*n[0], P[0] + 0.32*n[0]],
            [P[1] - 0.32*n[1], P[1] + 0.32*n[1]], color="#2b2b2b", lw=5,
            solid_capstyle="round", zorder=7)
    ax.plot(*P, "o", ms=10, color="#2b2b2b", zorder=8)

    # ── Resulting moment about the CG (sense from r x d) ───────────────
    cross = P[0]*d[1] - P[1]*d[0]                    # >0: counter-clockwise
    R = 0.55
    a0, a1 = (200, 330) if cross > 0 else (330, 200)
    th = np.deg2rad(np.linspace(a0, a1, 60))
    arc = np.column_stack([R*np.cos(th), R*np.sin(th)])
    ax.plot(arc[:-3, 0], arc[:-3, 1], color=C_TAU, lw=3.2, zorder=6)
    arrow(ax, arc[-4], arc[-1], C_TAU, 3.2, ms=24, z=6)
    ax.text(0.15, -0.95, r"$\boldsymbol{\tau}_i = \mathbf{r}_i \times u_i\,\mathbf{d}_i$",
            fontsize=22, color=C_TAU, ha="left", va="center")

    ax.set_xlim(-A*1.25, A*1.25); ax.set_ylim(-C - 0.75, C + 1.75)
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT_BASENAME}.{ext}", dpi=300, bbox_inches="tight",
                    transparent=TRANSPARENT)
    print(f"Saved {OUT_BASENAME}.png and {OUT_BASENAME}.svg")
    plt.show()


if __name__ == "__main__":
    main()
