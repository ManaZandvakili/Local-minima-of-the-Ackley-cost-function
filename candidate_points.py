"""
Slide 5 figure: the K = 30 candidate mounting locations (spherical
Fibonacci lattice) on the three reference airship hulls.

Candidate points use the same construction as the paper and Final.py:
    z_i     = 1 - (2i + 1) / K
    theta_i = i * gamma,   gamma ~= 137.5 deg (golden angle), wrapped to [0, 2pi)
    p_i     = [sqrt(1 - z_i^2) cos(theta_i), sqrt(1 - z_i^2) sin(theta_i), z_i]
and for an ellipsoid with semi-axes (a, b, c):
    p_i^ellipsoid = [a p_ix, b p_iy, c p_iz]

Usage:
    python candidate_points.py   -> shows the figure and saves candidate_points.png
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401

K = 30
HULLS = [
    ("Sphere",            (1.575, 1.575, 1.575)),
    ("Prolate (Cigar)",   (1.575, 1.575, 3.938)),
    ("Oblate (Flat Disk)", (3.150, 3.150, 0.788)),
]
SAME_SCALE    = True    # True: all panels share axis limits (true relative size)
SHOW_INDICES  = False   # True: label each point with its index i
SHOW_AXES     = True    # False: hide axes/grid for a cleaner slide image
HULL_COLOR    = "lightblue"
POINT_COLOR   = "#1f3b73"
OUT_FILE      = "candidate_points.png"


def fibonacci_unit_sphere(K):
    i = np.arange(K)
    z = 1.0 - (2.0*i + 1.0) / K
    gamma = np.pi * (3.0 - np.sqrt(5.0))           # golden angle ~= 137.5 deg
    theta = np.mod(i * gamma, 2.0*np.pi)
    r = np.sqrt(np.maximum(1.0 - z*z, 0.0))
    return np.column_stack([r*np.cos(theta), r*np.sin(theta), z])


def ellipsoid_surface(a, b, c, n_u=80, n_v=60):
    u = np.linspace(0.0, 2.0*np.pi, n_u)
    v = np.linspace(0.0, np.pi, n_v)
    x = a * np.outer(np.cos(u), np.sin(v))
    y = b * np.outer(np.sin(u), np.sin(v))
    z = c * np.outer(np.ones_like(u), np.cos(v))
    return x, y, z


def main():
    p_unit = fibonacci_unit_sphere(K)
    global_lim = max(max(s) for _, s in HULLS) * 1.08

    fig = plt.figure(figsize=(18, 6.5))
    for k, (name, (a, b, c)) in enumerate(HULLS):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")

        x, y, z = ellipsoid_surface(a, b, c)
        ax.plot_surface(x, y, z, color=HULL_COLOR, alpha=0.30,
                        linewidth=0, antialiased=True, shade=True)

        pts = p_unit * np.array([a, b, c])
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=45, color=POINT_COLOR,
                   edgecolors="white", linewidths=0.6, depthshade=True)
        if SHOW_INDICES:
            for i, (px, py, pz) in enumerate(pts):
                ax.text(px*1.08, py*1.08, pz*1.08, str(i), fontsize=7)

        lim = global_lim if SAME_SCALE else max(a, b, c) * 1.08
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=20, azim=-60)
        if SHOW_AXES:
            ax.set_xlabel("X [m]"); ax.set_ylabel("Y [m]"); ax.set_zlabel("Z [m]")
        else:
            ax.set_axis_off()
        ax.set_title(f"{name}\nK = {K}", fontsize=14, fontweight="bold")

    fig.legend(handles=[Line2D([], [], marker="o", linestyle="none",
                               markerfacecolor=POINT_COLOR,
                               markeredgecolor="white", markersize=9,
                               label="Candidate mounting location")],
               loc="lower center", fontsize=12, frameon=False)

    plt.tight_layout()
    fig.savefig(OUT_FILE, dpi=300, bbox_inches="tight")
    print(f"Saved {OUT_FILE}")
    plt.show()


if __name__ == "__main__":
    main()
