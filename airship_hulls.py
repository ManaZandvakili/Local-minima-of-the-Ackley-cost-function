"""
Plot the three reference airship hulls (sphere, prolate "cigar", oblate)
side by side in one figure, without propellers.

Semi-axes are the physical values from the paper (metres):
    Sphere : (1.575, 1.575, 1.575)
    Cigar  : (1.575, 1.575, 3.938)
    Oblate : (3.150, 3.150, 0.788)

All three panels share the same axis limits, so their relative sizes are
true to scale.  Set SHOW_KEEPOUT = False to hide the 35° gondola keep-out
region at the bottom of each hull.

Usage:
    python airship_hulls.py          -> shows the figure and saves airship_hulls.png
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401

HULLS = [
    ("(a) Sphere", (1.575, 1.575, 1.575)),
    ("(b) Cigar (prolate)", (1.575, 1.575, 3.938)),
    ("(c) Oblate", (3.150, 3.150, 0.788)),
]
SHOW_KEEPOUT        = True
GONDOLA_KEEPOUT_DEG = 35.0                         # half-angle around -Z
HULL_COLOR          = "lightblue"
KEEPOUT_COLOR       = "salmon"
OUT_FILE            = "airship_hulls.png"


def ellipsoid_surface(a, b, c, theta_max=np.pi, n_u=80, n_v=60):
    """Ellipsoid patch; theta is measured from +Z (0 = top, pi = bottom)."""
    u = np.linspace(0.0, 2.0*np.pi, n_u)
    v = np.linspace(0.0, theta_max, n_v)
    x = a * np.outer(np.cos(u), np.sin(v))
    y = b * np.outer(np.sin(u), np.sin(v))
    z = c * np.outer(np.ones_like(u), np.cos(v))
    return x, y, z


def keepout_cap(a, b, c, half_angle_deg, n_u=80, n_v=20):
    """Bottom cap of the ellipsoid inside the gondola keep-out cone.

    Like the optimisation code, the cone is applied in unit-sphere space
    (positions divided by the semi-axes), then mapped back to the hull.
    """
    u = np.linspace(0.0, 2.0*np.pi, n_u)
    v = np.linspace(np.pi - np.deg2rad(half_angle_deg), np.pi, n_v)
    x = a * np.outer(np.cos(u), np.sin(v))
    y = b * np.outer(np.sin(u), np.sin(v))
    z = c * np.outer(np.ones_like(u), np.cos(v))
    return x, y, z


def main():
    lim = max(max(ax) for _, ax in HULLS) * 1.05
    fig = plt.figure(figsize=(18, 6.5))

    for i, (label, (a, b, c)) in enumerate(HULLS):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")

        x, y, z = ellipsoid_surface(a, b, c)
        ax.plot_surface(x, y, z, color=HULL_COLOR, alpha=0.35,
                        linewidth=0, antialiased=True, shade=True)

        if SHOW_KEEPOUT:
            xk, yk, zk = keepout_cap(a, b, c, GONDOLA_KEEPOUT_DEG)
            ax.plot_surface(xk, yk, zk, color=KEEPOUT_COLOR, alpha=0.45,
                            linewidth=0, antialiased=True)

        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("X [m]"); ax.set_ylabel("Y [m]"); ax.set_zlabel("Z [m]")
        ax.set_title(f"{label}\nsemi-axes = ({a:.3f}, {b:.3f}, {c:.3f}) m",
                     fontsize=12)
        ax.view_init(elev=22, azim=-60)

    if SHOW_KEEPOUT:
        from matplotlib.patches import Patch
        fig.legend(handles=[Patch(color=HULL_COLOR, alpha=0.6, label="Hull"),
                            Patch(color=KEEPOUT_COLOR, alpha=0.6,
                                  label=f"Gondola keep-out ({GONDOLA_KEEPOUT_DEG:.0f}°)")],
                   loc="lower center", ncol=2, fontsize=11, frameon=False)

    plt.tight_layout()
    fig.savefig(OUT_FILE, dpi=300, bbox_inches="tight")
    print(f"Saved {OUT_FILE}")
    plt.show()


if __name__ == "__main__":
    main()
