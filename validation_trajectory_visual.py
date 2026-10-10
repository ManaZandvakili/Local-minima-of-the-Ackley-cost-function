"""
Slide 10 visual: circular path-following test used for closed-loop validation.

Shows the 50 m reference circle, the airship tracking it at 0.5 m/s with
its body x-axis aligned with the direction of travel, and the mean wind
(1.5 m/s) applied from 8 directions 45 deg apart, one direction per run.

Usage:
    python validation_trajectory_visual.py
        -> saves validation_trajectory_visual.png (transparent) and .svg
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch

R            = 50.0        # path radius [m]
N_WIND       = 8
SHIP_ANGLE   = 35.0        # where the airship is drawn on the circle [deg]
TRANSPARENT  = True
OUT_BASENAME = "validation_trajectory_visual"

C_PATH = "#1f3b73"; C_SHIP = "#9fd0e6"; C_SHIP_E = "#4f7f99"
C_VEL  = "#e07b00"; C_WIND = "#2a9d8f"; C_TXT = "#222222"


def arrow(ax, p0, p1, color, lw, ms=20, z=5, style="-|>"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms,
                                 color=color, lw=lw, zorder=z,
                                 shrinkA=0, shrinkB=0))


def main():
    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    ax.set_aspect("equal"); ax.axis("off")

    # ── Reference path ────────────────────────────────────────────────
    ax.add_patch(Circle((0, 0), R, fill=False, ec=C_PATH, lw=3.0,
                        ls=(0, (7, 4)), zorder=2))
    for a in (100, 190, 0):                       # direction-of-travel chevrons
        t = np.deg2rad(a)
        p = R*np.array([np.cos(t), np.sin(t)])
        tang = np.array([-np.sin(t), np.cos(t)])
        arrow(ax, p - 3*tang, p + 3*tang, C_PATH, 2.5, ms=22, z=3)

    # start point (t = 0): bottom of the circle, heading +x
    ax.plot(0, -R, "o", ms=11, color=C_PATH, mec="white", mew=1.5, zorder=4)
    ax.text(0, -R + 4.5, "start", fontsize=13, color=C_PATH, ha="center", va="bottom")

    # radius annotation
    arrow(ax, (0, 0), (R*np.cos(np.deg2rad(-40)), R*np.sin(np.deg2rad(-40))),
          "#555555", 1.6, ms=16, z=3)
    ax.plot(0, 0, "+", color="#555555", ms=12, mew=2)
    ax.text(R*0.50*np.cos(np.deg2rad(-40)) + 2, R*0.50*np.sin(np.deg2rad(-40)) - 5,
            "R = 50 m", fontsize=16, color="#444444", rotation=-40,
            ha="center", va="center")

    # ── Airship on the path, nose along the direction of travel ───────
    t  = np.deg2rad(SHIP_ANGLE)
    P  = R*np.array([np.cos(t), np.sin(t)])
    tg = np.array([-np.sin(t), np.cos(t)])
    heading = np.degrees(np.arctan2(tg[1], tg[0]))
    ax.add_patch(Ellipse(P, 13, 6, angle=heading, fc=C_SHIP, ec=C_SHIP_E,
                         lw=2, zorder=6))
    for s in (+1, -1):                                 # tail fins
        base = P - 6.0*tg
        n = np.array([-tg[1], tg[0]])
        ax.add_patch(plt.Polygon([base + 1.0*s*n, base - 2.2*tg + 3.2*s*n,
                                  base - 2.8*tg + 3.2*s*n, base - 1.8*tg + 0.4*s*n],
                                 closed=True, fc="#8fa3ae", ec=C_SHIP_E, zorder=5))
    arrow(ax, P + 6.5*tg, P + 20*tg, C_VEL, 4, ms=26, z=7)
    tip = P + 21*tg
    ax.text(tip[0] - 1, tip[1] + 2.5, "0.5 m/s", fontsize=16, color=C_VEL,
            ha="center", va="bottom", fontweight="bold")
    ax.text(P[0] + 9, P[1] - 4, r"$+x_{body}$ aligned" "\n" "with travel",
            fontsize=12, color=C_TXT, ha="left", va="top")

    # ── Mean wind from 8 directions ───────────────────────────────────
    r_out, r_in = R + 30, R + 12
    for k in range(N_WIND):
        a = np.deg2rad(k * 360.0 / N_WIND)
        u = np.array([np.cos(a), np.sin(a)])
        # wind blowing TOWARD angle a (as in Final.py), drawn on its upwind side
        arrow(ax, -r_out*u, -r_in*u, C_WIND, 3.2, ms=22, z=4)
        lab = -(r_out + 6)*u
        ax.text(lab[0], lab[1], f"{k*360//N_WIND}°", fontsize=13, color=C_WIND,
                ha="center", va="center")

    ax.text(0, -R - 46,
            "Mean wind 1.5 m/s from 8 directions (45° apart)\n"
            r"+ Gaussian gusts, $\sigma$ = 0.1 m/s",
            fontsize=15, color=C_WIND, ha="center", va="top")

    lim = R + 42
    ax.set_xlim(-lim, lim); ax.set_ylim(-R - 62, lim)
    plt.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT_BASENAME}.{ext}", dpi=300, bbox_inches="tight",
                    transparent=TRANSPARENT)
    print(f"Saved {OUT_BASENAME}.png and {OUT_BASENAME}.svg")
    plt.show()


if __name__ == "__main__":
    main()
