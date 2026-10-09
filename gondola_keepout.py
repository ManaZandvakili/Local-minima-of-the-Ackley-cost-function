"""
Gondola keep-out zone on the spherical airship.

A candidate mounting position p_i is rejected when its normalised
direction lies inside a 35 deg cone around the downward gondola axis:
        (p_i / |p_i|) . g  >  cos(35 deg),     g = [0, 0, -1]^T

Left : 3-D sphere with the keep-out cap, the K = 30 Fibonacci candidates
       (allowed vs rejected) and the gondola.
Right: side cross-section (X-Z plane) showing the 35 deg cone around g.

Usage:
    python gondola_keepout.py  -> shows the figure and saves gondola_keepout.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Wedge, Circle, FancyArrowPatch, FancyBboxPatch, Arc
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401

R               = 1.575         # sphere radius [m]
K               = 30
KEEPOUT_DEG     = 35.0
SHOW_CANDIDATES = True          # False: hull + keep-out zone only
OUT_FILE        = "gondola_keepout.png"

C_HULL  = "lightblue"; C_EDGE = "#6a8fa6"
C_KO    = "#e0533d"            # keep-out
C_OK    = "#1f3b73"            # allowed candidate
C_GOND  = "#6f7a82"


def fibonacci_unit_sphere(K):
    i = np.arange(K)
    z = 1.0 - (2.0*i + 1.0) / K
    theta = np.mod(i * np.pi * (3.0 - np.sqrt(5.0)), 2.0*np.pi)
    r = np.sqrt(np.maximum(1.0 - z*z, 0.0))
    return np.column_stack([r*np.cos(theta), r*np.sin(theta), z])


def sphere_patch(theta_from, theta_to, n_u=90, n_v=50):
    """Sphere surface for polar angle theta (from +Z) in [theta_from, theta_to]."""
    u = np.linspace(0, 2*np.pi, n_u)
    v = np.linspace(theta_from, theta_to, n_v)
    return (R*np.outer(np.cos(u), np.sin(v)),
            R*np.outer(np.sin(u), np.sin(v)),
            R*np.outer(np.ones_like(u), np.cos(v)))


def panel_3d(ax):
    ko = np.deg2rad(KEEPOUT_DEG)
    ax.plot_surface(*sphere_patch(0, np.pi - ko), color=C_HULL, alpha=0.30,
                    linewidth=0, antialiased=True)
    ax.plot_surface(*sphere_patch(np.pi - ko, np.pi), color=C_KO, alpha=0.55,
                    linewidth=0, antialiased=True)

    # boundary circle of the keep-out cap
    u = np.linspace(0, 2*np.pi, 200)
    rb, zb = R*np.sin(ko), -R*np.cos(ko)
    ax.plot(rb*np.cos(u), rb*np.sin(u), np.full_like(u, zb), color=C_KO, lw=2)

    # cone edges from the centre to the cap boundary
    for a in np.deg2rad([0, 90, 180, 270]):
        ax.plot([0, rb*np.cos(a)], [0, rb*np.sin(a)], [0, zb],
                color=C_KO, lw=1.0, ls="--", alpha=0.8)
    ax.quiver(0, 0, 0, 0, 0, -1.12*R, color="black", lw=2,
              arrow_length_ratio=0.12)
    ax.text(0.15, 0.15, -1.12*R, r"$\mathbf{g}$", fontsize=16)

    # gondola under the hull (small cylinder)
    u = np.linspace(0, 2*np.pi, 40); h = np.linspace(-R-0.55, -R-0.25, 2)
    U, H = np.meshgrid(u, h)
    ax.plot_surface(0.30*np.cos(U), 0.30*np.sin(U), H, color=C_GOND, alpha=0.9,
                    linewidth=0)

    if SHOW_CANDIDATES:
        p = fibonacci_unit_sphere(K)
        bad = p @ np.array([0, 0, -1.0]) > np.cos(ko)
        P = R*p
        ax.scatter(*P[~bad].T, s=40, color=C_OK, edgecolors="white",
                   linewidths=0.6, depthshade=False)
        ax.scatter(*P[bad].T, s=70, marker="X", color=C_KO,
                   edgecolors="black", linewidths=0.6, depthshade=False)

    lim = R*1.05
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim - 0.25, lim - 0.25)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=-14, azim=-55)
    ax.set_axis_off()


def panel_2d(ax):
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Circle((0, 0), R, fc=C_HULL, ec=C_EDGE, lw=2, alpha=0.6))
    # keep-out wedge: 35 deg either side of straight down (270 deg)
    ax.add_patch(Wedge((0, 0), R, 270 - KEEPOUT_DEG, 270 + KEEPOUT_DEG,
                       fc=C_KO, alpha=0.30, ec="none"))
    ax.add_patch(Arc((0, 0), 2*R, 2*R, theta1=270 - KEEPOUT_DEG,
                     theta2=270 + KEEPOUT_DEG, color=C_KO, lw=5))
    ko = np.deg2rad(KEEPOUT_DEG)
    for s in (-1, 1):
        ax.plot([0, s*R*np.sin(ko)], [0, -R*np.cos(ko)], color=C_KO, lw=1.8, ls="--")

    ax.add_patch(FancyArrowPatch((0, 0), (0, -1.18*R), arrowstyle="-|>",
                                 mutation_scale=22, color="black", lw=2.2, zorder=5))
    ax.text(0.12, -1.12*R, r"$\mathbf{g}=[0,0,-1]^T$", fontsize=15, va="center")
    ax.add_patch(Arc((0, 0), 0.9*R, 0.9*R, theta1=270, theta2=270 + KEEPOUT_DEG,
                     color="black", lw=1.5))
    mid = np.deg2rad(270 + KEEPOUT_DEG/2)
    ax.text(0.58*R*np.cos(mid), 0.58*R*np.sin(mid), f"{KEEPOUT_DEG:.0f}°",
            fontsize=15, ha="left", va="center")
    ax.plot(0, 0, "o", color="black", ms=6)

    # gondola
    ax.add_patch(FancyBboxPatch((-0.40, -R - 0.80), 0.80, 0.30,
                                boxstyle="round,pad=0.02,rounding_size=0.08",
                                fc=C_GOND, ec="#4d565c", zorder=4))
    ax.text(-0.45, -R - 0.65, "Gondola", fontsize=13, ha="right", va="center")

    ax.text(0, 0.35*R, "Allowed mounting region", fontsize=14,
            ha="center", color=C_OK)
    ax.text(-1.05*R, -0.95*R, "Keep-out\nregion", fontsize=14, ha="right",
            va="center", color=C_KO, fontweight="bold")
    ax.text(0, -1.80*R,
            r"reject $\mathbf{p}_i$  if  $\hat{\mathbf{p}}_i\cdot\mathbf{g} > \cos 35^\circ$",
            fontsize=15, ha="center")
    ax.set_xlim(-1.9*R, 1.9*R); ax.set_ylim(-1.95*R, 1.1*R)


def main():
    fig = plt.figure(figsize=(15, 7))
    panel_3d(fig.add_subplot(1, 2, 1, projection="3d"))
    panel_2d(fig.add_subplot(1, 2, 2))

    handles = [Patch(color=C_HULL, alpha=0.6, label="Hull"),
               Patch(color=C_KO, alpha=0.55,
                     label=f"Gondola keep-out ({KEEPOUT_DEG:.0f}° cone)")]
    if SHOW_CANDIDATES:
        handles += [Line2D([], [], marker="o", ls="none", color=C_OK, ms=8,
                           label="Allowed candidate"),
                    Line2D([], [], marker="X", ls="none", color=C_KO,
                           markeredgecolor="black", ms=10,
                           label="Rejected candidate")]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               fontsize=12, frameon=False)
    plt.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUT_FILE, dpi=300, bbox_inches="tight")
    print(f"Saved {OUT_FILE}")
    plt.show()


if __name__ == "__main__":
    main()
