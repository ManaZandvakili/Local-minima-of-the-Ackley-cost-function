"""
Refinement (post-evolutionary-sweep + user pruning):
The cost function now excludes:
    J1 — authority (sum of singular values)
    J2 — residual force
    J6 — direction symmetry
    J7 — fault tolerance
    J9 — aerodynamic tangential
Active: J3 (#props), J4 (effort, L∞-norm peak per-motor demand),
J5 (position symmetry), J8 (geometric feasibility soft penalty).
The metrics are still computed and printed for diagnostic purposes —
only the cost-function weighting is affected.

IMPORTANT: the gondola hard-gate constraints (F10/F12/F13) remain fully
active, so physical safety is preserved even without J8's soft penalty.
See ACTIVE_OBJ_MASK and EXCLUDED_OBJ_INDICES near the top of the file.

Also disabled in this version (for runtime):
    • Mission-energy independent metric
    • Monte Carlo robustness
    • Weight sensitivity sweep
The closed-loop path-following validation is now the sole external
performance test.

Fixes over v7:
    F13  gondola_position_hard_check now accepts shape + shape_params
         and maps positions to unit-sphere space (divide by semi-axes)
         before computing the keepout angle.  On the sphere this is
         identical to v7.  On the oblate flat disk (a=b=2, c=0.5) the
         south-pole region covers a MUCH larger footprint in Cartesian
         space: a position at z≈−0.47 near the disk centre has raw
         angle ≈56° from −Z (passed v7 check) but unit-sphere angle
         ≈20° (correctly rejected by v8).  evaluate() now passes
         self.shape and self.shape_params to the gate.

Retained from v7:
    F12  Hard position rejection added alongside F10's hard direction
         rejection.  If ANY propeller mounting point falls within
         GONDOLA_KEEPOUT_DEG of the bottom pole (-Z), evaluate()
         immediately returns penalty_inactive.  Previously only
         directions were hard-gated; positions only had a soft J8
         penalty which the optimizer overrode whenever J1 benefit was
         high, causing propellers to physically appear inside the
         gondola keep-out zone (visible in the 3-D plot).

Retained from v6:
    F10  Hard direction rejection — any thrust vector within 35° of
         −Z → 1e6.
    F11  boxplot() 'labels' → 'tick_labels' (Matplotlib ≥3.9).

Retained from v5:
    F6   Gondola keep-out soft penalty on positions AND directions
         (now serves as a secondary guard; F10 is the primary hard
         gate).
    F7   J1 base weight raised to 2× (weight-sensitivity result).
    F8   set_xticks fix + ASCII Delta in plot_weight_sensitivity.
    F9   plt.close('all') after every plt.show().

Retained from v4:
    F2   Monte Carlo re-projects onto the correct surface
         (sphere/ellipsoid).
    F3   PyVista gondola placed at −Z to match GONDOLA_DOWN = [0,0,−1].
    F4   Second Pareto plot uses log₁₀ scale on the J7 axis.
    F5   MOEA/D per-trial weight vectors for genuine Pareto coverage.

Physical reference
    Airship  : STARLab, spherical, 3.15 m diameter (radius 1.575 m)
    Propeller: Kavan 2826-1200, 7×5 inch blades (diameter 0.178 m)
"""

# ============================================================
# Part 0 — Physical constants
# ============================================================
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D          # noqa
from dataclasses import dataclass, field

AIRSHIP_DIAMETER_M  = 3.15
AIRSHIP_RADIUS_M    = AIRSHIP_DIAMETER_M / 2.0           # 1.575 m
PROP_DIAMETER_M     = 7 * 0.0254                         # 0.1778 m (7-inch)
D_MIN_PHYSICAL_M    = 2.5 * PROP_DIAMETER_M              # 0.445 m
D_MIN_SPHERE        = D_MIN_PHYSICAL_M / AIRSHIP_RADIUS_M  # ≈ 0.282

ANG_MIN_DEG         = 40.0        # I4 — minimum thrust-direction diversity
COND_THRESHOLD      = 100.0       # I5 — condition number hard threshold
GONDOLA_DOWN        = np.array([0.0, 0.0, -1.0])  # I6 — bottom pole (-Z)
GONDOLA_KEEPOUT_DEG = 35.0        # I6 — half-angle of keep-out cone
DOWNWASH_CONE_DEG   = 15.0        # I7 — thrust-wake cone half-angle

# Monte Carlo tolerance (I8)
MC_POS_SIGMA_M      = 0.020       # ±2 cm position tolerance
MC_POS_SIGMA        = MC_POS_SIGMA_M / AIRSHIP_RADIUS_M   # normalised ≈ 0.013
MC_ANG_SIGMA_DEG    = 2.0         # ±2° direction tolerance
MC_N_TRIALS         = 200

# ── Path-following dynamic validation (independent metric, NOT in cost) ──
# A closed-loop simulation in which the airship must track a 5 m circular
# trajectory at 1 m/s against a 3 m/s constant wind plus gusts.  PSO and the
# best-authority random configuration face identical conditions; the only
# difference is the wrench matrix M they were dealt.
SIM_MASS_KG     = 10.0          # airship mass [kg]
SIM_INERTIA     = 9.0           # yaw moment of inertia [kg·m²]
SIM_U_MAX       = 2.5           # max thrust per propeller [N]   (Kavan 7-inch)
SIM_DRAG        = 0.5           # quadratic drag coefficient
SIM_PATH_R      = 50.0          # circular path radius [m]
SIM_CRUISE_V    = 0.5           # cruise speed [m/s]
SIM_WIND        = np.array([1.5, 0.0])     # constant wind [m/s] (+X)
SIM_GUST_SIGMA  = 1.0           # gust RMS [m/s]
SIM_T_TOTAL     = 650.0         # simulation duration [s] (~1 full lap on 50 m circle at 0.5 m/s)
SIM_DT          = 0.05          # timestep [s] (20 Hz — plenty for a 0.5 m/s airship)
SIM_KP_POS, SIM_KD_POS = 2.0, 1.5     # outer position-PD gains
SIM_KP_YAW, SIM_KD_YAW = 3.0, 1.0     # yaw-PD gains
SIM_PATH_SEED   = 7             # fixed wind/gust seed

try:
    import pyvista as pv
    pv.global_theme.font.family = "arial"
    PYVISTA_OK = True
except ImportError:
    PYVISTA_OK = False
    print("[Warning] PyVista not found — pip install pyvista")

# ============================================================
# Active objective set (cost-function refinement)
# ============================================================
# Empirical findings from the J1→J9 evolutionary sweep
# (Final_evolutionary.py) plus current pruning:
#   • J2 (residual force) — L1 and L2 produced bit-for-bit identical
#     configurations on all 3 shapes.  J2 is redundant with the rank-6
#     hard gate, which already drives ‖I − M·M⁺‖ to zero.
#   • J6 (direction symmetry) — currently commented out per user request.
#   • J7 (fault tolerance) — only weak individual effect across L6→L7;
#     does not meaningfully change the best configurations.
#   • J8 (geometric feasibility soft penalty) — currently commented out
#     per user request.  Note: the gondola keepout HARD gates
#     (F10/F12/F13) remain active, so physical safety is still
#     guaranteed.
# All four are excluded from the cost calculation.  The metrics are still
# computed and printed in Pareto reports for diagnostic purposes.
EXCLUDED_OBJ_INDICES = (0, 1, 6, 8)                  # J1, J2, J6, J7, J9 (0-idx)
ACTIVE_OBJ_MASK      = np.ones(9, dtype=bool)
ACTIVE_OBJ_MASK[list(EXCLUDED_OBJ_INDICES)] = False     # False at idx 0, 1, 5, 6, 8

# ============================================================
# Part 1 — Utilities & geometry
# ============================================================
def unit(v, eps=1e-12):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, eps)

def spherical_to_cart(theta, phi):
    return np.stack([np.cos(theta)*np.sin(phi),
                     np.sin(theta)*np.sin(phi),
                     np.cos(phi)], axis=-1)

def fibonacci_directions(K):
    i     = np.arange(K)
    ga    = 2 - (1 + 5**0.5) / 2
    theta = 2 * np.pi * ((i * ga) % 1.0)
    z     = 1 - (2*i + 1) / K
    r     = np.sqrt(np.maximum(1 - z*z, 0))
    return unit(np.stack([r*np.cos(theta), r*np.sin(theta), z], axis=-1))

def fibonacci_sphere(K, radius=1.0):
    return radius * fibonacci_directions(K)

def fibonacci_ellipsoid(K, a=1.0, b=1.0, c=2.5):
    d = fibonacci_directions(K)
    return np.column_stack([a*d[:,0], b*d[:,1], c*d[:,2]])

def surface_normals(positions, shape, shape_params):
    if shape == "sphere":
        return positions / np.linalg.norm(positions, axis=1, keepdims=True)
    a, b, c = shape_params
    g = np.column_stack([positions[:,0]/(a*a),
                         positions[:,1]/(b*b),
                         positions[:,2]/(c*c)])
    return g / np.linalg.norm(g, axis=1, keepdims=True)

def project_onto_surface(positions, shape, shape_params):
    """Radial projection of positions onto the correct surface."""
    if shape == "sphere":
        return positions / np.linalg.norm(positions, axis=1, keepdims=True)
    a, b, c = shape_params
    abc = np.array([a, b, c])
    # map to unit sphere, normalise, map back
    scaled = positions / abc
    scaled = scaled / np.linalg.norm(scaled, axis=1, keepdims=True)
    return scaled * abc

def d_min_for_shape(shape, shape_params):
    """I6 — scale d_min by geometric mean of ellipsoid semi-axes."""
    if shape == "sphere":
        return D_MIN_SPHERE
    a, b, c = shape_params
    return D_MIN_SPHERE * (a * b * c) ** (1.0/3.0)

def pareto_filter(points):
    N = points.shape[0]
    dominated = np.zeros(N, dtype=bool)
    for i in range(N):
        if dominated[i]:
            continue
        for j in range(N):
            if i == j:
                continue
            if np.all(points[j] <= points[i]) and np.any(points[j] < points[i]):
                dominated[i] = True
                break
    idx = np.where(~dominated)[0]
    return points[idx], idx

# F4: log_axes=(axis_indices,...) applies log10 transform before plotting
def plot_pareto_3d(points, labels=("J1","J2","J3"),
                   title="Pareto Front 3D", log_axes=()):
    fig = plt.figure(figsize=(8, 6))
    ax  = fig.add_subplot(111, projection="3d")
    plot_pts = points.astype(float).copy()
    axis_labels = list(labels)
    for i in log_axes:
        plot_pts[:, i] = np.log10(np.maximum(plot_pts[:, i], 1e-9))
        axis_labels[i] = f"log₁₀({labels[i]})"
    ax.scatter(plot_pts[:,0], plot_pts[:,1], plot_pts[:,2], s=40)
    ax.set_xlabel(axis_labels[0], fontsize=12)
    ax.set_ylabel(axis_labels[1], fontsize=12)
    ax.set_zlabel(axis_labels[2], fontsize=12)
    ax.set_title(title, fontsize=14)
    plt.tight_layout(); plt.show(); plt.close('all')

# ============================================================
# Part 2 — Metrics J1–J9
# ============================================================
def build_M_matrix(positions, directions):
    U   = directions.T
    RxU = np.cross(positions, directions, axis=1).T
    return np.vstack([U, RxU])

def controllability_metrics(M, dof_weights, positions, directions):
    _, S, _ = np.linalg.svd(M + 1e-9*np.random.randn(*M.shape), full_matrices=False)
    rank   = int(np.sum(S > 1e-6))
    sv_sum = float(np.sum(S))
    cond   = float(S.max() / max(S.min(), 1e-9))
    try:
        Minv = np.linalg.pinv(M, rcond=1e-9)
        I6 = np.eye(6)
        A    = Minv @ I6
        res_norms = np.linalg.norm(I6 - M @ A, axis=0)
        residual  = float(np.sum(dof_weights * res_norms))
        # J4 = peak per-motor demand per DOF command (L∞ norm of allocation
        # columns), not the L2 norm.  This is what the saturation constraint
        # actually penalises: one motor at the limit fails the sim, no matter
        # how relaxed the other five are.
        effort    = float(np.sum(dof_weights * np.linalg.norm(A, axis=0, ord=np.inf)))
    except Exception:
        residual = effort = 1e6
    sym_pos = float(np.linalg.norm(np.mean(positions,  axis=0))) if positions.size > 0 else 1e6
    sym_dir = float(np.linalg.norm(np.mean(directions, axis=0))) if positions.size > 0 else 1e6
    return dict(rank=rank, sv_sum=sv_sum, cond=cond,
                residual=residual, effort=effort,
                sym_pos=sym_pos, sym_dir=sym_dir, S=S)

def fault_tolerance_metric(positions, directions):
    """N-1 + N-2 redundancy (C5 from v2)."""
    N = positions.shape[0]
    if N <= 6:
        return 1e6
    penalty = 0.0
    for i in range(N):
        pi = np.delete(positions, i, 0)
        di = np.delete(directions, i, 0)
        Si = np.linalg.svd(build_M_matrix(pi, di), full_matrices=False)[1]
        penalty += 1e4 if np.sum(Si > 1e-6) < 6 else 1.0 / max(Si.min(), 1e-9)
    if N > 8:
        for i in range(N):
            for j in range(i+1, N):
                k = [x for x in range(N) if x != i and x != j]
                Sij = np.linalg.svd(build_M_matrix(positions[k], directions[k]),
                                    full_matrices=False)[1]
                penalty += 0.5*1e4 if np.sum(Sij > 1e-6) < 6 else 0.5 / max(Sij.min(), 1e-9)
    return float(penalty)

# ── gondola keep-out (I6) ────────────────────────────────────
def gondola_keepout_penalty(positions, directions=None):
    """
    I6 / F6: two-part keep-out penalty.

    Position part  — penalise props whose mounting point is within
    GONDOLA_KEEPOUT_DEG of the bottom pole (GONDOLA_DOWN = -Z).

    Direction part — penalise props whose thrust vector points within
    GONDOLA_KEEPOUT_DEG of straight down, regardless of where they sit.
    This eliminates axis-snapped [0,0,-1] directions that the position
    check alone cannot catch.
    """
    norms     = positions / np.linalg.norm(positions, axis=1, keepdims=True)
    cos_ang   = norms @ GONDOLA_DOWN
    angles    = np.arccos(np.clip(cos_ang, -1, 1))
    keepout_r = np.deg2rad(GONDOLA_KEEPOUT_DEG)
    penalty   = float(np.sum(np.maximum(keepout_r - angles, 0)))

    if directions is not None:
        # cos(GONDOLA_KEEPOUT_DEG) ≈ 0.819 — threshold for "pointing down"
        cos_thresh = np.cos(keepout_r)
        cos_dir    = directions @ GONDOLA_DOWN   # 1.0 = straight down
        penalty   += float(np.sum(np.maximum(cos_dir - cos_thresh, 0)))

    return penalty

def gondola_position_hard_check(positions, shape="sphere", shape_params=None):
    """
    F12/F13: hard rejection gate — returns True if ANY propeller mounting
    position falls within GONDOLA_KEEPOUT_DEG of the bottom pole (-Z).

    F13 fix: positions are first mapped to unit-sphere coordinates by
    dividing by the ellipsoid semi-axes before computing the angle.
    On a sphere (a=b=c=1) this is identical to the old behaviour.
    On a flat disk (a=b=2, c=0.5) the bottom-centre region has a large
    XY footprint: raw Cartesian normalisation gave ~56° (passed), but
    unit-sphere normalisation gives ~20° (correctly rejected).
    """
    if len(positions) == 0:
        return False
    if shape == "sphere" or shape_params is None:
        unit_pos = positions / np.linalg.norm(positions, axis=1, keepdims=True)
    else:
        a, b, c  = shape_params
        scaled   = positions / np.array([a, b, c])           # → unit sphere
        unit_pos = scaled / np.linalg.norm(scaled, axis=1, keepdims=True)
    cos_ang    = unit_pos @ GONDOLA_DOWN          # 1.0 = at south pole
    cos_thresh = np.cos(np.deg2rad(GONDOLA_KEEPOUT_DEG))   # ≈ 0.819
    return bool(np.any(cos_ang > cos_thresh))

def gondola_direction_hard_check(directions):
    """
    F10: hard rejection gate — returns True if ANY propeller direction
    points within GONDOLA_KEEPOUT_DEG of straight down (GONDOLA_DOWN).

    Called inside evaluate() before any matrix work.  Violating
    solutions are immediately returned with penalty_inactive so they can
    never become gbest or pbest regardless of how good J1 would
    otherwise be.

    The existing soft penalty in gondola_keepout_penalty (J8) is kept as
    a secondary gradient signal to steer particles away from the
    boundary.
    """
    if len(directions) == 0:
        return False
    cos_thresh = np.cos(np.deg2rad(GONDOLA_KEEPOUT_DEG))   # ≈ 0.819
    cos_dir    = directions @ GONDOLA_DOWN                  # 1.0 = straight down
    return bool(np.any(cos_dir > cos_thresh))

# ── downwash interference (I7) ───────────────────────────────
def downwash_interference_metric(positions, directions):
    """
    I7: penalise pairs where one propeller sits inside the thrust wake
    cone (half-angle DOWNWASH_CONE_DEG) of another.
    """
    N = positions.shape[0]
    if N < 2:
        return 0.0
    cos_thresh = np.cos(np.deg2rad(DOWNWASH_CONE_DEG))
    penalty    = 0.0
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            r_ij  = positions[j] - positions[i]
            dist  = np.linalg.norm(r_ij)
            if dist < 1e-9:
                continue
            cos_a = np.dot(directions[i], r_ij / dist)
            if cos_a > cos_thresh:
                penalty += (cos_a - cos_thresh) / (1.0 - cos_thresh)
    return float(penalty)

# ── J8 combined (C3,C4,I6,I7) ────────────────────────────────
def geometric_constraints_metric(positions, directions,
                                 shape="sphere", shape_params=None):
    N = positions.shape[0]
    if N < 2:
        return 0.0
    d_min   = d_min_for_shape(shape, shape_params)
    ang_min = np.deg2rad(ANG_MIN_DEG)
    penalty = 0.0
    for i in range(N):
        for j in range(i+1, N):
            dij = np.linalg.norm(positions[i] - positions[j])
            if dij < d_min:
                penalty += (d_min - dij)
            dot = np.clip(np.dot(directions[i], directions[j]), -1, 1)
            ang = np.arccos(dot)
            if ang < ang_min:
                penalty += (ang_min - ang)
    penalty += gondola_keepout_penalty(positions, directions)   # F6
    penalty += downwash_interference_metric(positions, directions)
    return float(penalty)

def aerodynamic_tangential_metric(positions, directions,
                                  shape="sphere", shape_params=None):
    if positions.size == 0:
        return 0.0
    norms = surface_normals(positions, shape, shape_params)
    return float(np.sum(np.abs(np.sum(directions * norms, axis=1))))

# ============================================================
# Part 3 — InnerPSO
# ============================================================
@dataclass
class InnerPSOConfig:
    swarm_size       : int   = 30
    iters            : int   = 100
    w_start          : float = 0.90
    w_end            : float = 0.40
    c_cog            : float = 1.6
    c_soc            : float = 1.6
    sel_min          : float = -5.0
    sel_max          : float =  5.0
    theta_min        : float =  0.0
    theta_max        : float =  2 * np.pi
    phi_min          : float =  0.0
    phi_max          : float =  np.pi
    threshold        : float =  0.5
    penalty_inactive : float =  1e6
    penalty_rank     : float =  1e4
    # Option B: dof_weights re-balanced toward planar control
    # [Fx Fy Fz Mx My Mz].  Airships do not need much vertical-force
    # (Fz) or pitch/roll (Mx, My) authority — buoyancy carries weight
    # and passive aerodynamic stability handles pitch/roll.
    # Concentrating PSO's authority budget on the three DOFs actually
    # exercised in flight (Fx, Fy, Mz) yields better trajectory tracking
    # without changing the metric we evaluate against.
    dof_weights      : np.ndarray = field(
        default_factory=lambda: np.array([2.0, 2.0, 0.3, 0.3, 0.3, 2.0])
    )

class InnerPSO:
    def __init__(self, cfg, positions, weights,
                 shape="sphere", shape_params=None, seed=42):
        self.cfg          = cfg
        self.positions    = positions
        self.weights      = weights
        self.shape        = shape
        self.shape_params = shape_params
        self.dim          = positions.shape[0] * 3
        self.rng          = np.random.default_rng(seed)
        self.K            = positions.shape[0]

        self.lb = np.empty(self.dim); self.ub = np.empty(self.dim)
        self.lb[0::3] = cfg.sel_min;   self.ub[0::3] = cfg.sel_max
        self.lb[1::3] = cfg.theta_min; self.ub[1::3] = cfg.theta_max
        self.lb[2::3] = cfg.phi_min;   self.ub[2::3] = cfg.phi_max

        self.V_max = 0.2 * (self.ub - self.lb)

        self.X         = self.rng.uniform(self.lb, self.ub,
                                          size=(cfg.swarm_size, self.dim))
        self.V         = np.zeros_like(self.X)
        self.pbest     = self.X.copy()
        self.pbest_val = np.full(cfg.swarm_size, np.inf)
        self.gbest     = None
        self.gbest_val = np.inf

        self.obj_history     = []
        self.eval_positions  = []
        self.eval_directions = []
        self.eval_masks      = []

    @staticmethod
    def _sigmoid(x):
        return 1.0 / (1.0 + np.exp(-x))

    def decode(self, x):
        logits = x[0::3]; theta = x[1::3]; phi = x[2::3]
        mask  = self._sigmoid(logits) >= self.cfg.threshold
        dirs  = spherical_to_cart(theta, phi)
        if int(np.sum(mask)) < 6:
            mask[np.argsort(logits)[-6:]] = True
        pos_a  = self.positions[mask]
        dirs_a = dirs[mask].copy()
        if pos_a.shape[0] > 0:
            norms = surface_normals(pos_a, self.shape, self.shape_params)
            dirs_a[np.sum(dirs_a * norms, axis=1) < 0] *= -1
        return mask, pos_a, dirs_a

    def evaluate(self, x):
        cfg = self.cfg; w = self.weights
        mask, pos, dirs = self.decode(x)
        if pos.shape[0] == 0:
            return cfg.penalty_inactive, np.zeros(9), mask, np.zeros((0,3)), np.zeros((0,3))

        # F12/F13: hard gate — reject any position inside the gondola cone
        #   F13: pass shape so check is done in unit-sphere space (oblate fix)
        if gondola_position_hard_check(pos, self.shape, self.shape_params):
            return cfg.penalty_inactive, np.zeros(9), mask, pos, dirs

        # F10: hard gate — reject any direction pointing into the gondola cone
        if gondola_direction_hard_check(dirs):
            return cfg.penalty_inactive, np.zeros(9), mask, pos, dirs

        M   = build_M_matrix(pos, dirs)
        met = controllability_metrics(M, cfg.dof_weights, pos, dirs)
        N   = float(pos.shape[0]); K = float(self.K)

        J1 = -met["sv_sum"]; J2 = met["residual"]; J3 = N; J4 = met["effort"]
        J5 = met["sym_pos"]; J6 = met["sym_dir"]
        J7 = fault_tolerance_metric(pos, dirs)
        J8 = geometric_constraints_metric(pos, dirs,
                shape=self.shape, shape_params=self.shape_params)
        J9 = aerodynamic_tangential_metric(pos, dirs,
                shape=self.shape, shape_params=self.shape_params)

        n2  = max(N*(N-1)/2.0, 1.0) if N > 8 else 0.0
        j7s = max(N*5.0 + n2*2.5, 1.0)
        scales = np.array([10., 1., K/2., 1., 1., 1., j7s, 1., max(N,1.)])
        J_norm = np.array([J1, J2, J3, J4, J5, J6, J7, J8, J9]) / scales

        # Apply active-objective mask: zero-weight excluded Js, then
        # renormalise so the remaining weights still sum to 1 (cost
        # magnitudes stay comparable to the pre-exclusion baseline).
        w_active = w * ACTIVE_OBJ_MASK
        w_active = w_active / max(w_active.sum(), 1e-12)
        score    = float(np.dot(w_active, J_norm))

        if met["rank"] < 6:
            score += cfg.penalty_rank * (6 - met["rank"])

        if met["cond"] <= COND_THRESHOLD:
            score += 0.01 * met["cond"] / COND_THRESHOLD
        else:
            score += 0.5 * np.log10(met["cond"] / COND_THRESHOLD)

        return float(score), np.array([J1,J2,J3,J4,J5,J6,J7,J8,J9]), mask, pos, dirs

    def _eval_and_update_bests(self, i):
        val, Js, mask, pos, dirs = self.evaluate(self.X[i])
        self.obj_history.append(Js); self.eval_positions.append(pos)
        self.eval_directions.append(dirs); self.eval_masks.append(mask)
        if val < self.pbest_val[i]:
            self.pbest_val[i] = val; self.pbest[i] = self.X[i].copy()
            if val < self.gbest_val:
                self.gbest_val = val; self.gbest = self.X[i].copy()

    def step(self, iter_idx):
        w = (self.cfg.w_start
             - (self.cfg.w_start - self.cfg.w_end)
             * iter_idx / max(self.cfg.iters - 1, 1))

        for i in range(self.cfg.swarm_size):
            self._eval_and_update_bests(i)

        rp = self.rng.random(self.X.shape); rg = self.rng.random(self.X.shape)
        self.V = (w * self.V
                  + self.cfg.c_cog * rp * (self.pbest - self.X)
                  + self.cfg.c_soc * rg * (self.gbest  - self.X))
        self.V = np.clip(self.V, -self.V_max, self.V_max)
        self.X = np.clip(self.X + self.V, self.lb, self.ub)

    def run(self):
        convergence = []
        for it in range(self.cfg.iters):
            self.step(it)
            convergence.append(self.gbest_val)
        for i in range(self.cfg.swarm_size):
            self._eval_and_update_bests(i)
        return (self.gbest_val, np.array(convergence), self.gbest,
                self.eval_positions, self.eval_directions,
                self.eval_masks, self.obj_history)

# ============================================================
# Part 4 — run_multiple_trials  (F5: MOEA/D per-trial weights)
# ============================================================
def _scalarisation_weights(num_trials, n_obj=9, seed=0):
    """
    F5 — MOEA/D-style weight vectors.  Dirichlet(alpha=0.5) draws cluster
    near the simplex boundaries, giving good coverage of trade-off
    extremes so each PSO trial genuinely targets a different region of
    the Pareto front.

    Excluded objectives (ACTIVE_OBJ_MASK is False at their indices) are
    forced to zero weight and the remaining columns renormalised, so
    each trial only explores trade-offs across the active subset.
    """
    rng = np.random.default_rng(seed)
    w   = rng.dirichlet(np.full(n_obj, 0.5), size=num_trials)
    if n_obj == ACTIVE_OBJ_MASK.size:
        w = w * ACTIVE_OBJ_MASK                              # zero excluded
        w = w / np.maximum(w.sum(axis=1, keepdims=True), 1e-12)
    return w

def run_multiple_trials(positions, weights, num_trials=20,
                        shape="sphere", shape_params=None):
    all_Js = []; all_best = []; all_pos = []; all_dirs = []
    all_masks = []; all_convergences = []; trial_costs = []
    best_pos_list = []; best_dirs_list = []

    # F5: generate diverse per-trial weight vectors instead of reusing
    # the same equal-weight vector for every trial
    trial_weights = _scalarisation_weights(num_trials, n_obj=len(weights))

    for t in range(num_trials):
        print(f"  Trial {t+1}/{num_trials}")
        cfg = InnerPSOConfig(iters=80)
        pso = InnerPSO(cfg, positions, trial_weights[t],
                       shape=shape, shape_params=shape_params, seed=t*17+3)
        (best_cost, conv, best_x,
         pos_list, dir_list, mask_list, Js_list) = pso.run()
        all_convergences.append(conv); trial_costs.append(float(best_cost))
        all_pos.extend(pos_list); all_dirs.extend(dir_list)
        all_masks.extend(mask_list); all_Js.extend(Js_list)

        mask, pos, dirs = pso.decode(best_x)
        M = build_M_matrix(pos, dirs)
        met = controllability_metrics(M, cfg.dof_weights, pos, dirs)
        J_best = np.array([
            -met["sv_sum"], met["residual"], float(pos.shape[0]),
             met["effort"], met["sym_pos"], met["sym_dir"],
             fault_tolerance_metric(pos, dirs),
             geometric_constraints_metric(pos, dirs, shape=shape, shape_params=shape_params),
             aerodynamic_tangential_metric(pos, dirs, shape=shape, shape_params=shape_params),
        ])
        all_best.append(J_best)
        best_pos_list.append(pos); best_dirs_list.append(dirs)

    return (np.array(all_Js), np.array(all_best),
            np.array(all_pos, dtype=object), np.array(all_dirs, dtype=object),
            np.array(all_masks, dtype=object),
            best_pos_list, best_dirs_list,
            all_convergences, np.array(trial_costs))

# ============================================================
# Part 5 — Ablation helpers
# ============================================================
def run_full_inner_pso(positions, weights, shape="sphere", shape_params=None):
    cfg = InnerPSOConfig(iters=100)
    pso = InnerPSO(cfg, positions, weights, shape=shape, shape_params=shape_params)
    best_cost, conv, best_x, *_ = pso.run()
    return best_cost, conv, best_x

def run_ablation(positions, weights, shape="sphere", shape_params=None):
    labels = ["J1","J2","J3","J4","J5","J6","J7","J8","J9"]
    results = {}
    print("  Running full ...")
    bf, cf, _ = run_full_inner_pso(positions, weights, shape=shape, shape_params=shape_params)
    results["Full"] = (bf, cf)
    for k in range(9):
        lbl = f"no{labels[k]}"
        w_new = weights.copy(); w_new[k] = 0.0
        s = np.sum(w_new)
        if s > 1e-12:
            w_new /= s
        print(f"  Running {lbl} ...")
        b, c, _ = run_full_inner_pso(positions, w_new, shape=shape, shape_params=shape_params)
        results[lbl] = (b, c)
    return results

def plot_convergence_ablation(results):
    plt.figure(figsize=(10,6))
    for key, (_, conv) in results.items():
        plt.plot(conv, label=key)
    plt.xlabel("Iteration"); plt.ylabel("Best Cost")
    plt.title("Convergence — Full vs Ablations")
    plt.legend(fontsize=7); plt.grid(True); plt.tight_layout()
    plt.show(); plt.close('all')

def print_ablation_summary(results):
    print("\n=== Ablation Summary ===")
    fc = results["Full"][0]
    print(f"  {'Run':<10}  {'Best Cost':>12}  {'Delta':>10}")
    print("  " + "-"*36)
    for key, (cost, _) in results.items():
        print(f"  {key:<10}  {cost:>12.4f}  {cost-fc:>+10.4f}")

# ============================================================
# Part 6 — Result visualizations
# ============================================================
def plot_config_3d(positions, directions, shape, shape_params, title="Best Config"):
    fig = plt.figure(figsize=(9,7)); ax = fig.add_subplot(111, projection="3d")
    u = np.linspace(0, 2*np.pi, 40); v = np.linspace(0, np.pi, 20)
    if shape == "sphere":
        a = b = c = 1.0
        xs = np.outer(np.cos(u), np.sin(v))
        ys = np.outer(np.sin(u), np.sin(v))
        zs = np.outer(np.ones_like(u), np.cos(v))
        arrow_scale = 0.30
    else:
        a, b, c = shape_params
        xs = a*np.outer(np.cos(u), np.sin(v))
        ys = b*np.outer(np.sin(u), np.sin(v))
        zs = c*np.outer(np.ones_like(u), np.cos(v))
        arrow_scale = 0.55
    ax.plot_surface(xs, ys, zs, alpha=0.10, color="skyblue", linewidth=0)
    # shade gondola keep-out zone at -Z pole (matches GONDOLA_DOWN)
    phi_ko = np.deg2rad(GONDOLA_KEEPOUT_DEG)
    u2 = np.linspace(0, 2*np.pi, 30)
    p2 = np.linspace(np.pi - phi_ko, np.pi, 10)
    ax.plot_surface(a*np.outer(np.cos(u2), np.sin(p2)),
                    b*np.outer(np.sin(u2), np.sin(p2)),
                    c*np.outer(np.ones_like(u2), np.cos(p2)),
                    alpha=0.25, color="salmon", linewidth=0)
    ax.scatter(positions[:,0], positions[:,1], positions[:,2],
               s=90, c="red", zorder=5, label=f"{len(positions)} props")
    for p, d in zip(positions, directions):
        ax.quiver(p[0], p[1], p[2],
                  d[0]*arrow_scale, d[1]*arrow_scale, d[2]*arrow_scale,
                  color="darkorange", linewidth=1.5, arrow_length_ratio=0.35)
    ax.set_box_aspect([a, b, c])
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title(title, fontsize=13)
    ax.legend(handles=[
        plt.scatter([], [], c="red", s=60, label=f"{len(positions)} props"),
        plt.Rectangle((0,0), 1, 1, fc="salmon", alpha=0.5, label="Gondola keep-out")
    ], fontsize=9)
    plt.tight_layout(); plt.show(); plt.close('all')

def plot_dof_authority(positions, directions, title="6-DOF Authority"):
    M = build_M_matrix(positions, directions)
    Minv = np.linalg.pinv(M, rcond=1e-9)
    dof_labels = ["Fx\n(surge)", "Fy\n(sway)", "Fz\n(heave)",
                  "Roll\n(Mx)", "Pitch\n(My)", "Yaw\n(Mz)"]
    auth = np.array([1.0/max(np.linalg.norm(Minv @ np.eye(6)[:, d]), 1e-9)
                     for d in range(6)])
    from matplotlib.patches import Patch
    fig, ax_bar = plt.subplots(figsize=(8,5))
    colors = ["steelblue"]*3 + ["coral"]*3
    bars = ax_bar.bar(dof_labels, auth, color=colors, edgecolor="white", linewidth=0.8)
    ax_bar.set_ylabel("Authority (1/effort)", fontsize=11)
    ax_bar.grid(axis="y", alpha=0.4)
    ax_bar.set_title(title, fontsize=13)
    for bar, val in zip(bars, auth):
        ax_bar.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.02,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=9)
    ax_bar.legend(handles=[Patch(color="steelblue", label="Force"),
                           Patch(color="coral", label="Torque")], fontsize=9)
    plt.tight_layout(); plt.show(); plt.close('all')

def plot_convergence_trials(convergences, label):
    arr = np.array(convergences)
    fig, ax = plt.subplots(figsize=(10,5))
    for row in arr:
        ax.plot(row, alpha=0.35, linewidth=1.0, color="steelblue")
    ax.plot(arr.mean(axis=0), color="red", linewidth=2.5, label="Mean")
    ax.set_xlabel("Iteration", fontsize=11); ax.set_ylabel("Best Cost", fontsize=11)
    ax.set_title(f"{label} — Convergence across {len(convergences)} trials", fontsize=13)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3); plt.tight_layout()
    plt.show(); plt.close('all')

def plot_cost_distribution(costs_dict):
    labels = list(costs_dict.keys()); data = [costs_dict[k] for k in labels]
    fig, ax = plt.subplots(figsize=(9,5))
    palette = ["steelblue", "coral", "mediumseagreen"]
    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.5)
    for patch, col in zip(bp["boxes"], palette):
        patch.set_facecolor(col); patch.set_alpha(0.7)
    for med in bp["medians"]:
        med.set_color("black"); med.set_linewidth(2)
    for i, d in enumerate(data):
        ax.text(i+1, np.median(d), f" {np.median(d):.2f}", va="center", fontsize=9)
    ax.set_ylabel("Final Best Cost", fontsize=11)
    ax.set_title("PSO Consistency — Final Cost Distribution per Shape", fontsize=12)
    ax.grid(axis="y", alpha=0.4); plt.tight_layout()
    plt.show(); plt.close('all')

# ============================================================
# Independent performance metric — Mission Energy (NOT in cost)
# ============================================================
# Benchmark mission profile: a fixed sequence of wrench commands every
# configuration must execute.  For each segment, we compute the minimum-
# norm thrust allocation u = M⁺·τ and integrate the quadratic motor
# power Σuᵢ² over the segment duration.
MISSION_PROFILE = [
    # (wrench command [Fx, Fy, Fz, Mx, My, Mz],  duration_seconds)
    (np.array([5.0,  0.0,  0.0,  0.0,  0.0,  0.0]), 60.0),  # hover vs side wind
    (np.array([2.0,  0.0,  0.0,  0.0,  0.0,  0.0]),  5.0),  # forward cruise
    (np.array([0.0,  0.0,  0.0,  0.0,  0.0,  1.0]),  2.0),  # yaw maneuver
    (np.array([0.0,  0.0, -1.0,  0.0,  0.5,  0.0]),  3.0),  # descent + pitch
]

def mission_energy(M, motor_const=1.0):
    """Independent metric: total energy to execute the benchmark mission."""
    if M.size == 0 or M.shape[1] < 6:
        return np.inf
    try:
        Minv = np.linalg.pinv(M, rcond=1e-9)
        total_E = 0.0
        for tau, dt in MISSION_PROFILE:
            u = Minv @ tau
            total_E += float(np.sum(u**2)) * motor_const * dt
        return float(total_E)
    except Exception:
        return np.inf

def run_random_baseline(positions, weights, shape, shape_params,
                        n_target=1500, max_attempts=200000):
    """
    Constrained random baseline (smart random).

    Provides a statistically meaningful baseline by sampling ONLY from
    the feasible region:

      1) Positions inside the gondola keepout cone are pre-filtered out
         (their selection logits are forced to sel_min so they cannot be
         chosen by the sigmoid threshold).  Mapping to unit-sphere space
         first (F13) so the test is correct on oblate/prolate shapes.

      2) The polar angle phi of every direction is sampled from
         [0, arccos(-cos(35°))] ≈ [0, 145°], so the raw thrust vector
         cannot point within 35° of -Z.

      3) Any residual infeasible draws (rare) are rejected and
         resampled, continuing until n_target valid configurations are
         collected (or max_attempts is reached as a safety limit).
    """
    cfg = InnerPSOConfig()
    pso = InnerPSO(cfg, positions, weights,
                   shape=shape, shape_params=shape_params)
    rng = np.random.default_rng(99)

    # ── 1) forbidden-position mask (unit-sphere space, matches F13) ─────
    if shape == "sphere" or shape_params is None:
        unit_pos = positions / np.linalg.norm(positions, axis=1, keepdims=True)
    else:
        a, b, c = shape_params
        scaled   = positions / np.array([a, b, c])
        unit_pos = scaled / np.linalg.norm(scaled, axis=1, keepdims=True)
    cos_thresh    = np.cos(np.deg2rad(GONDOLA_KEEPOUT_DEG))    # ≈ 0.819
    forbidden_pos = (unit_pos @ GONDOLA_DOWN) > cos_thresh

    # ── 2) safe phi bound for raw directions ────────────────────────────
    phi_safe_max = np.arccos(-cos_thresh)        # ≈ 2.531 rad ≈ 145°
    lb_safe      = pso.lb.copy()
    ub_safe      = pso.ub.copy()
    ub_safe[2::3] = phi_safe_max                  # restrict every phi

    # ── 3) sample loop with rejection of residual edge cases ────────────
    RAND_SAMPLES_PER_N_CAP = 25
    costs       = []
    energies    = []
    attempts    = 0
    best_auth   = -np.inf
    best_pos    = None
    best_dirs   = None
    best_by_N   = {}    # N → (auth, pos, dirs)  — best per N (oracle)
    all_by_N    = {}    # N → list of (auth, pos, dirs) — capped, for median
    while len(costs) < n_target and attempts < max_attempts:
        x = rng.uniform(lb_safe, ub_safe)
        x[0::3][forbidden_pos] = cfg.sel_min      # never select keepout positions
        cost, _, _, pos, dirs = pso.evaluate(x)
        if pos.shape[0] > 0 and cost < 1e5:
            costs.append(cost)
            M_rand = build_M_matrix(pos, dirs)
            energies.append(mission_energy(M_rand))
            auth = float(np.linalg.svd(M_rand, full_matrices=False)[1].sum())
            n_props = pos.shape[0]
            if n_props not in best_by_N or auth > best_by_N[n_props][0]:
                best_by_N[n_props] = (auth, pos.copy(), dirs.copy())
            bucket = all_by_N.setdefault(n_props, [])
            if len(bucket) < RAND_SAMPLES_PER_N_CAP:
                bucket.append((auth, pos.copy(), dirs.copy()))
            if auth > best_auth:
                best_auth, best_pos, best_dirs = auth, pos.copy(), dirs.copy()
        attempts += 1
    return (np.array(costs), np.array(energies), attempts,
            best_pos, best_dirs, best_by_N, all_by_N)

def plot_mission_energy(pso_energy, pareto_energies, random_energies, label,
                        random_attempts=None):
    """Independent-metric comparison plot: Mission energy distribution."""
    finite_rand = random_energies[np.isfinite(random_energies)]
    x_max = max(np.max(pareto_energies)*2.0, np.percentile(finite_rand, 70)*1.5)
    rz = finite_rand[finite_rand <= x_max]

    fig, ax = plt.subplots(figsize=(10,5))

    if random_attempts is not None and random_attempts > 0:
        rate = 100.0 * len(finite_rand) / random_attempts
        rand_label = (f"Smart Random (n={len(finite_rand)} valid / "
                      f"{random_attempts} tries, {rate:.1f}%, ≤{x_max:.1f})")
    else:
        rand_label = f"Random (n={len(finite_rand)}, ≤{x_max:.1f})"

    ax.hist(rz, bins=60, alpha=0.55, color="gray", label=rand_label)
    ax.axvline(np.min(finite_rand), color="dimgray", linestyle="--", linewidth=1.8,
               label=f"Random best={np.min(finite_rand):.2f}")
    ax.axvline(np.median(finite_rand), color="black", linestyle=":", linewidth=1.5,
               label=f"Random median={np.median(finite_rand):.2f}")
    for e in pareto_energies:
        ax.axvline(e, color="steelblue", alpha=0.4, linewidth=1.2)
    ax.axvline(pso_energy, color="red", linewidth=2.5,
               label=f"PSO best energy={pso_energy:.2f}")
    ax.axvline(np.mean(pareto_energies), color="orange", linewidth=2, linestyle="--",
               label=f"PSO Pareto mean={np.mean(pareto_energies):.2f}")
    ax.set_xlim(0, x_max)
    ax.set_xlabel("Mission energy  Σuᵢ²·dt  (dimensionless)", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title(f"{label} — Independent Metric: Mission Energy "
                 f"(NOT used in optimization cost)", fontsize=12)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    plt.tight_layout(); plt.show(); plt.close('all')

def plot_pso_vs_random(trial_costs, random_costs, label, random_attempts=None):
    x_max = max(np.max(trial_costs)*3.0, np.percentile(random_costs, 30)*3.0)
    rz    = random_costs[random_costs <= x_max]
    fig, ax = plt.subplots(figsize=(10,5))
    if random_attempts is not None and random_attempts > 0:
        rate = 100.0 * len(random_costs) / random_attempts
        rand_label = (f"Smart Random (n={len(random_costs)} valid / "
                      f"{random_attempts} tries, {rate:.1f}%, ≤{x_max:.1f})")
    else:
        rand_label = f"Random (n={len(random_costs)}, ≤{x_max:.1f})"
    ax.hist(rz, bins=60, alpha=0.55, color="gray", label=rand_label)
    ax.axvline(np.min(random_costs), color="dimgray", linestyle="--", linewidth=1.8,
               label=f"Random best={np.min(random_costs):.2f}")
    ax.axvline(np.median(random_costs), color="black", linestyle=":", linewidth=1.5,
               label=f"Random median={np.median(random_costs):.2f}")
    for c in trial_costs:
        ax.axvline(c, color="steelblue", alpha=0.4, linewidth=1.2)
    ax.axvline(np.min(trial_costs), color="red", linewidth=2.5,
               label=f"PSO best={np.min(trial_costs):.2f}")
    ax.axvline(np.mean(trial_costs), color="orange", linewidth=2, linestyle="--",
               label=f"PSO mean={np.mean(trial_costs):.2f}")
    ax.set_xlim(0, x_max)
    ax.set_xlabel("Cost", fontsize=11); ax.set_ylabel("Count", fontsize=11)
    ax.set_title(f"{label} — PSO vs Smart-Random (constrained to feasible region)",
                 fontsize=12)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    plt.tight_layout(); plt.show(); plt.close('all')

# ============================================================
# Part 6c — Monte Carlo robustness (I8) & Weight sensitivity (I10)
# ============================================================
def monte_carlo_robustness(positions, directions,
                           pos_sigma=MC_POS_SIGMA,
                           ang_sigma_deg=MC_ANG_SIGMA_DEG,
                           n_trials=MC_N_TRIALS,
                           shape="sphere", shape_params=None):
    """
    F2 — I8: 200 random perturbations of ±2 cm (normalised) on
    positions and ±2° on thrust directions.  Perturbed positions are
    re-projected onto the correct surface (sphere or ellipsoid) instead
    of always normalising to the unit sphere.
    """
    rng = np.random.default_rng(42)
    M0  = build_M_matrix(positions, directions)
    S0  = np.linalg.svd(M0, full_matrices=False)[1]
    sv0 = float(np.sum(S0))
    ang_sigma = np.deg2rad(ang_sigma_deg)

    sv_sums = []; ranks = []; min_svs = []
    for _ in range(n_trials):
        pos_p = positions + rng.normal(0, pos_sigma, positions.shape)
        pos_p = project_onto_surface(pos_p, shape, shape_params)

        dirs_p = directions.copy()
        for k in range(len(dirs_p)):
            perp = rng.normal(0, 1, 3)
            perp -= np.dot(perp, dirs_p[k]) * dirs_p[k]
            pn   = np.linalg.norm(perp)
            if pn > 1e-9:
                perp /= pn
                a = rng.normal(0, ang_sigma)
                dirs_p[k] = (dirs_p[k]*np.cos(a)
                             + np.cross(perp, dirs_p[k])*np.sin(a)
                             + perp*np.dot(perp, dirs_p[k])*(1-np.cos(a)))
                dirs_p[k] /= np.linalg.norm(dirs_p[k])

        Sp = np.linalg.svd(build_M_matrix(pos_p, dirs_p), full_matrices=False)[1]
        sv_sums.append(float(np.sum(Sp)))
        ranks.append(int(np.sum(Sp > 1e-6)))
        min_svs.append(float(Sp.min()))

    sv_sums = np.array(sv_sums); ranks = np.array(ranks); min_svs = np.array(min_svs)
    return dict(sv_sum_0=sv0, sv_mean=sv_sums.mean(), sv_std=sv_sums.std(),
                degradation_pct=(sv0 - sv_sums.mean())/sv0*100,
                rank6_frac=float(np.mean(ranks == 6)),
                min_sv_mean=min_svs.mean(), min_sv_std=min_svs.std(),
                sv_sums=sv_sums, ranks=ranks)

def plot_monte_carlo(mc, label):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].hist(mc["sv_sums"], bins=30, color="steelblue", alpha=0.7, edgecolor="white")
    axes[0].axvline(mc["sv_sum_0"], color="red", linewidth=2,
                    label=f"Nominal = {mc['sv_sum_0']:.3f}")
    axes[0].axvline(mc["sv_mean"], color="orange", linewidth=2, linestyle="--",
                    label=f"MC mean = {mc['sv_mean']:.3f}  "
                          f"({mc['degradation_pct']:.1f}% drop)")
    axes[0].set_xlabel("Sum of Singular Values", fontsize=11)
    axes[0].set_ylabel("Count", fontsize=11)
    axes[0].set_title("Controllability under Manufacturing Tolerances\n"
                      f"(±{MC_POS_SIGMA_M*100:.0f} cm pos, ±{MC_ANG_SIGMA_DEG}° dir)",
                      fontsize=11)
    axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

    vals = [mc["rank6_frac"]*100, (1-mc["rank6_frac"])*100]
    bars = axes[1].bar(["Rank 6\n(Full DOF)", "Rank < 6\n(Degraded)"],
                       vals, color=["steelblue", "coral"], alpha=0.8)
    axes[1].set_ylabel("% of perturbations", fontsize=11)
    axes[1].set_ylim(0, 110)
    axes[1].set_title("Full 6-DOF Rank Retention", fontsize=12)
    for bar, v in zip(bars, vals):
        axes[1].text(bar.get_x()+bar.get_width()/2, v+1,
                     f"{v:.1f}%", ha="center", fontsize=12)
    axes[1].grid(axis="y", alpha=0.3)

    fig.suptitle(f"{label} — Monte Carlo Robustness (n={MC_N_TRIALS})", fontsize=13)
    plt.tight_layout(); plt.show(); plt.close('all')

    print(f"\n  Monte Carlo Robustness [{label}]:")
    print(f"    Nominal sv_sum      : {mc['sv_sum_0']:.4f}")
    print(f"    MC mean sv_sum      : {mc['sv_mean']:.4f}  (±{mc['sv_std']:.4f})")
    print(f"    Degradation         : {mc['degradation_pct']:.2f}%")
    print(f"    Rank-6 retention    : {mc['rank6_frac']*100:.1f}%")
    print(f"    Min SV mean ± std   : {mc['min_sv_mean']:.4f} ± {mc['min_sv_std']:.4f}")

def run_weight_sensitivity(positions, shape, shape_params, n_trials=4):
    """
    I10: re-run PSO with each ACTIVE J-weight doubled in turn.
    Skips excluded objectives (J2, J7) — doubling a zero weight is a
    no-op.  Uses fewer iterations (60) and trials (4) to keep runtime
    manageable.
    """
    labels  = ["J1","J2","J3","J4","J5","J6","J7","J8","J9"]
    base_w  = np.ones(9) / 9.0
    results = {}
    cfg     = InnerPSOConfig(iters=60)

    print("  Weight sensitivity: baseline ...")
    costs = []
    for t in range(n_trials):
        pso = InnerPSO(cfg, positions, base_w,
                       shape=shape, shape_params=shape_params, seed=t*7+1)
        c, _, _, *_ = pso.run(); costs.append(c)
    results["Baseline"] = np.array(costs)

    for k in range(9):
        if not ACTIVE_OBJ_MASK[k]:
            continue                       # skip excluded J2, J7
        w_new = base_w.copy(); w_new[k] *= 2.0; w_new /= w_new.sum()
        print(f"  Weight sensitivity: 2x {labels[k]} ...")
        costs = []
        for t in range(n_trials):
            pso = InnerPSO(cfg, positions, w_new,
                           shape=shape, shape_params=shape_params, seed=t*7+1)
            c, _, _, *_ = pso.run(); costs.append(c)
        results[f"2x {labels[k]}"] = np.array(costs)

    return results

def plot_weight_sensitivity(results, label):
    keys  = list(results.keys())
    means = np.array([results[k].mean() for k in keys])
    stds  = np.array([results[k].std()  for k in keys])
    bm    = means[0]
    delta = (means - bm) / max(abs(bm), 1e-9) * 100

    colors = ["gray"] + ["coral" if d > 0 else "steelblue" for d in delta[1:]]
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.bar(keys, delta, color=colors, alpha=0.82, edgecolor="white")
    ax.errorbar(range(len(keys)), delta,
                yerr=stds/max(abs(bm), 1e-9)*100,
                fmt="none", color="black", capsize=4)
    ax.axhline(0, color="black", linewidth=1)
    ax.set_ylabel("Cost change vs baseline (%)", fontsize=11)
    ax.set_title(f"{label} — Weight Sensitivity (I10)\n"
                 "Red = worse when that J is up-weighted  "
                 "Blue = better (that J was being over-penalised)", fontsize=11)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys, rotation=30, ha="right", fontsize=9)
    ax.grid(axis="y", alpha=0.3); plt.tight_layout()
    plt.show(); plt.close('all')

    print(f"\n  Weight Sensitivity [{label}]:")
    print(f"  {'Config':<15}  {'Mean cost':>10}  {'Delta vs base':>14}")
    print("  " + "-"*43)
    for k, m, d in zip(keys, means, delta):
        print(f"  {k:<15}  {m:>10.4f}  {d:>+14.2f}%")

# ============================================================
# Path-following validation (closed-loop dynamic simulation)
# ============================================================
def _path_reference(t):
    """Smooth circular reference: (x_ref, y_ref, yaw_ref) at time t."""
    omega = SIM_CRUISE_V / SIM_PATH_R
    return (SIM_PATH_R * np.sin(omega * t),
            SIM_PATH_R * (1.0 - np.cos(omega * t)),
            omega * t)

def _wrap_pi(a):
    return (a + np.pi) % (2 * np.pi) - np.pi

def simulate_path_following(positions, directions,
                            seed=SIM_PATH_SEED,
                            wind_angle_rad=0.0,
                            kp_pos=None, kd_pos=None,
                            kp_yaw=None, kd_yaw=None,
                            gains_per_axis=None,
                            wind_mag_override=None,
                            t_total=None):
    """
    3-DOF planar closed-loop simulation.  The airship tracks a circular
    reference with a cascade PD controller; the controller produces a
    desired body-frame wrench [Fx_b, Fy_b, Mz] which is allocated to the
    propellers via M⁺ and clipped to [0, u_max].
    """
    KP_POS  = SIM_KP_POS  if kp_pos  is None else float(kp_pos)
    KD_POS  = SIM_KD_POS  if kd_pos  is None else float(kd_pos)
    KP_YAW  = SIM_KP_YAW  if kp_yaw  is None else float(kp_yaw)
    KD_YAW  = SIM_KD_YAW  if kd_yaw  is None else float(kd_yaw)
    WIND_M  = float(np.linalg.norm(SIM_WIND)) if wind_mag_override is None \
              else float(wind_mag_override)
    T_TOT   = SIM_T_TOTAL if t_total is None else float(t_total)

    if gains_per_axis is not None:
        g           = np.asarray(gains_per_axis, dtype=float).ravel()
        KP_XY       = np.array([g[0], g[2]])
        KD_XY       = np.array([g[1], g[3]])
        KP_YAW      = float(g[4])
        KD_YAW      = float(g[5])
        PER_AXIS    = True
    else:
        PER_AXIS    = False

    # Scale torque rows to physical N·m (positions are dimensionless × R)
    M_phys = build_M_matrix(positions, directions).copy()
    M_phys[3:] *= AIRSHIP_RADIUS_M
    A      = M_phys[[0, 1, 5], :]                  # Fx, Fy, Mz rows
    A_pinv = np.linalg.pinv(A, rcond=1e-9)

    rng   = np.random.default_rng(seed)
    state = np.zeros(6)                             # [x, y, vx, vy, yaw, yaw_rate]
    N     = int(T_TOT / SIM_DT)
    pos_h = np.zeros((N, 2))
    err_h = np.zeros(N)
    nrg_h = np.zeros(N)
    sat_n = 0
    e_run = 0.0

    for i in range(N):
        t = i * SIM_DT
        x_r, y_r, yaw_r = _path_reference(t)
        pos_err = np.array([x_r - state[0], y_r - state[1]])
        yaw_err = _wrap_pi(yaw_r - state[4])
        vel_err = -state[2:4]
        w_err   = -state[5]

        if PER_AXIS:
            F_des_w = SIM_MASS_KG * (KP_XY * pos_err + KD_XY * vel_err)
        else:
            F_des_w = SIM_MASS_KG * (KP_POS * pos_err + KD_POS * vel_err)
        Tau_des = SIM_INERTIA * (KP_YAW * yaw_err + KD_YAW * w_err)
        c, s    = np.cos(state[4]), np.sin(state[4])
        R_w2b   = np.array([[ c,  s], [-s,  c]])
        F_des_b = R_w2b @ F_des_w
        wrench  = np.array([F_des_b[0], F_des_b[1], Tau_des])

        u      = A_pinv @ wrench
        u_clip = np.clip(u, 0.0, SIM_U_MAX)
        sat_n += int(np.sum((u < 0.0) | (u > SIM_U_MAX)))

        wrench_act = A @ u_clip
        F_b, Tau_b = wrench_act[:2], wrench_act[2]
        F_w        = R_w2b.T @ F_b

        c_w, s_w   = np.cos(wind_angle_rad), np.sin(wind_angle_rad)
        wind_base  = WIND_M * np.array([c_w, s_w])
        wind  = wind_base + rng.normal(0.0, SIM_GUST_SIGMA * 0.1, 2)
        v_rel = state[2:4] - wind
        drag  = -SIM_DRAG * np.linalg.norm(v_rel) * v_rel

        accel = (F_w + drag) / SIM_MASS_KG
        alpha = Tau_b / SIM_INERTIA
        state[0:2] += state[2:4] * SIM_DT
        state[2:4] += accel * SIM_DT
        state[4]   += state[5] * SIM_DT
        state[5]   += alpha * SIM_DT

        e_run    += float(np.sum(u_clip**2)) * SIM_DT
        pos_h[i]  = state[0:2]
        err_h[i]  = float(np.linalg.norm(pos_err))
        nrg_h[i]  = e_run

    return dict(
        rms_error    = float(np.sqrt(np.mean(err_h**2))),
        max_error    = float(np.max(err_h)),
        total_energy = float(e_run),
        sat_pct      = float(100.0 * sat_n / max(N * A.shape[1], 1)),
        completed    = bool(np.max(err_h) < 3.0),
        traj         = pos_h,
        error_curve  = err_h,
        energy_curve = nrg_h,
    )

def plot_path_comparison(label, pso_sim, rand_sim, n_pso, n_rand):
    """Three-panel comparison: trajectory, tracking error, energy."""
    t     = np.arange(len(pso_sim["traj"])) * SIM_DT
    omega = SIM_CRUISE_V / SIM_PATH_R
    x_ref = SIM_PATH_R * np.sin(omega * t)
    y_ref = SIM_PATH_R * (1.0 - np.cos(omega * t))

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    ax = axes[0]
    ax.plot(x_ref, y_ref, "k--", linewidth=1.5, label="Reference")
    ax.plot(pso_sim["traj"][:,0], pso_sim["traj"][:,1], "b-", linewidth=2,
            label=f"PSO  ({n_pso} props, RMS={pso_sim['rms_error']:.2f} m)")
    ax.plot(rand_sim["traj"][:,0], rand_sim["traj"][:,1], "r-", linewidth=2,
            label=f"Random ({n_rand} props, RMS={rand_sim['rms_error']:.2f} m)")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_title("Trajectory")
    ax.axis("equal"); ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

    ax = axes[1]
    ax.plot(t, pso_sim["error_curve"], "b-", linewidth=1.5, label="PSO")
    ax.plot(t, rand_sim["error_curve"], "r-", linewidth=1.5, label="Random")
    ax.set_xlabel("time [s]"); ax.set_ylabel("position error [m]")
    ax.set_title(f"Tracking error  (PSO max={pso_sim['max_error']:.2f}, "
                 f"Random max={rand_sim['max_error']:.2f})")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

    ax = axes[2]
    ax.plot(t, pso_sim["energy_curve"], "b-", linewidth=1.5, label="PSO")
    ax.plot(t, rand_sim["energy_curve"], "r-", linewidth=1.5, label="Random")
    ax.set_xlabel("time [s]"); ax.set_ylabel("cumulative ∫u² dt  [N²·s]")
    ax.set_title(f"Energy  (PSO={pso_sim['total_energy']:.1f}, "
                 f"Random={rand_sim['total_energy']:.1f})")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=9)

    matched_tag = (f"matched N={n_pso}" if n_pso == n_rand
                   else f"PSO N={n_pso}, Random N={n_rand}")
    fig.suptitle(f"{label} — Matched-N closed-loop validation  ({matched_tag})\n"
                 f"{SIM_PATH_R:.0f} m circle, {SIM_CRUISE_V:.1f} m/s cruise, "
                 f"{SIM_WIND[0]:.1f} m/s wind  (only placement quality differs)",
                 fontsize=12)
    plt.tight_layout(); plt.show(); plt.close('all')

def simulate_multi_wind(positions, directions, n_angles=8,
                        seed=SIM_PATH_SEED, gains=None,
                        per_angle_seeds=None):
    """
    Run the path-following simulation across n_angles equally-spaced
    wind directions.
    """
    angles = np.linspace(0.0, 2.0*np.pi, n_angles, endpoint=False)

    if per_angle_seeds is None:
        per_angle_seeds = [int(seed) + i for i in range(n_angles)]

    kp_p = kd_p = kp_y = kd_y = None
    g_axis = None
    if gains is not None:
        gains_arr = np.asarray(gains).ravel()
        if gains_arr.size == 6:
            g_axis = gains_arr
        elif gains_arr.size == 4:
            kp_p, kd_p, kp_y, kd_y = [float(v) for v in gains_arr]

    runs = [simulate_path_following(
                positions, directions,
                seed=int(per_angle_seeds[i]),
                wind_angle_rad=a,
                kp_pos=kp_p, kd_pos=kd_p, kp_yaw=kp_y, kd_yaw=kd_y,
                gains_per_axis=g_axis)
            for i, a in enumerate(angles)]
    rms = np.array([r['rms_error']    for r in runs])
    mxe = np.array([r['max_error']    for r in runs])
    eng = np.array([r['total_energy'] for r in runs])
    sat = np.array([r['sat_pct']      for r in runs])
    return dict(
        angles_deg          = np.degrees(angles),
        rms_error_per_angle = rms,
        max_error_per_angle = mxe,
        energy_per_angle    = eng,
        sat_per_angle       = sat,
        rms_mean            = float(rms.mean()),
        rms_worst           = float(rms.max()),
        rms_best            = float(rms.min()),
        rms_std             = float(rms.std()),
        max_mean            = float(mxe.mean()),
        max_worst           = float(mxe.max()),
        energy_mean         = float(eng.mean()),
        energy_worst        = float(eng.max()),
        sat_mean            = float(sat.mean()),
        all_completed       = all(r['completed'] for r in runs),
        runs                = runs,
    )

# ────────────────────────────────────────────────────────────────────────
# Patel-style PSO controller tuning  (integer-order PD)
# ────────────────────────────────────────────────────────────────────────
TUNE_T_TOTAL     = 100.0    # raised from 20 s — long enough to see curvature on the 50 m circle
TUNE_SEED        = 17

TUNE_LB          = np.array([1.0, 0.4, 1.0, 0.4, 1.5, 0.3])
TUNE_UB          = np.array([5.0, 3.0, 5.0, 3.0, 6.0, 2.5])

TUNE_W_TRACK     = 5.0
TUNE_W_ENERGY    = 1e-3

TUNE_PSO_SWARM   = 12
TUNE_PSO_ITERS   = 30       # reduced from 50 — Patel early-stop usually fires well before this
TUNE_PSO_W       = 0.7
TUNE_PSO_C1      = 1.5
TUNE_PSO_C2      = 1.5
TUNE_PSO_VMAX_K  = 0.2
TUNE_CONV_TOL    = 1e-3
TUNE_CONV_WIN    = 20

def _airship_controller_cost(positions, directions, gains_6d):
    """
    Compute the airship controller cost for one gain vector.
    Two-term objective on a wind-free 20-s reference run.
    """
    res = simulate_path_following(
        positions, directions,
        seed=TUNE_SEED, wind_angle_rad=0.0,
        gains_per_axis=gains_6d,
        wind_mag_override=0.0,
        t_total=TUNE_T_TOTAL,
    )
    cost = (TUNE_W_TRACK  * res['rms_error']
          + TUNE_W_ENERGY * res['total_energy'])
    return float(cost)

def tune_pd_gains(positions, directions, return_history=False):
    """
    Patel-style PSO controller tuning.  Searches a 6-D space (per-axis
    PD gains) using PSO with airship-specific objective.  Returns a
    6-tuple of gains.
    """
    rng = np.random.default_rng(TUNE_SEED)
    n_dim = 6

    X      = rng.uniform(TUNE_LB, TUNE_UB, size=(TUNE_PSO_SWARM, n_dim))
    V      = np.zeros_like(X)
    V_max  = TUNE_PSO_VMAX_K * (TUNE_UB - TUNE_LB)
    pbest      = X.copy()
    pbest_val  = np.array([_airship_controller_cost(positions, directions, x)
                            for x in X])
    g_i        = int(np.argmin(pbest_val))
    gbest      = pbest[g_i].copy()
    gbest_val  = float(pbest_val[g_i])

    history = [gbest_val]
    for it in range(TUNE_PSO_ITERS):
        r1 = rng.random((TUNE_PSO_SWARM, n_dim))
        r2 = rng.random((TUNE_PSO_SWARM, n_dim))
        V  = (TUNE_PSO_W*V
              + TUNE_PSO_C1*r1*(pbest - X)
              + TUNE_PSO_C2*r2*(gbest - X))
        V  = np.clip(V, -V_max, V_max)
        X  = np.clip(X + V, TUNE_LB, TUNE_UB)

        for i in range(TUNE_PSO_SWARM):
            val = _airship_controller_cost(positions, directions, X[i])
            if val < pbest_val[i]:
                pbest_val[i] = val
                pbest[i]     = X[i].copy()
                if val < gbest_val:
                    gbest_val = val
                    gbest     = X[i].copy()

        history.append(gbest_val)

        if len(history) > TUNE_CONV_WIN:
            window = history[-TUNE_CONV_WIN:]
            if (max(window) - min(window)) < TUNE_CONV_TOL:
                break

    gains_tuple = tuple(float(v) for v in gbest)
    if return_history:
        return gains_tuple, gbest_val, history
    return gains_tuple

# ── Constrained random sampling at EXACT N (same feasibility as PSO) ─────
def sample_random_at_n_feasible(target_n, positions, shape, shape_params,
                                n_target=20, max_attempts=20000, seed=42):
    """
    Generate up to n_target feasible random configurations with EXACTLY
    target_n propellers.  Each sample passes the same hard constraints
    PSO solutions must pass (gondola position/direction gates, min
    spacing, rank-6).
    """
    rng = np.random.default_rng(seed)

    if shape == "sphere" or shape_params is None:
        unit_pos = positions / np.linalg.norm(positions, axis=1, keepdims=True)
    else:
        a, b, c = shape_params
        scaled   = positions / np.array([a, b, c])
        unit_pos = scaled / np.linalg.norm(scaled, axis=1, keepdims=True)
    cos_thresh   = np.cos(np.deg2rad(GONDOLA_KEEPOUT_DEG))
    safe_pos_idx = np.where((unit_pos @ GONDOLA_DOWN) <= cos_thresh)[0]
    if len(safe_pos_idx) < target_n:
        return []

    phi_safe_max = np.arccos(-cos_thresh)
    d_min        = d_min_for_shape(shape, shape_params)

    samples  = []
    attempts = 0
    while len(samples) < n_target and attempts < max_attempts:
        attempts += 1
        idx = rng.choice(safe_pos_idx, size=target_n, replace=False)
        pos_sub = positions[idx]

        diffs = pos_sub[:, None, :] - pos_sub[None, :, :]
        d     = np.linalg.norm(diffs, axis=2)
        np.fill_diagonal(d, np.inf)
        if d.min() < d_min:
            continue

        theta = rng.uniform(0.0, 2.0*np.pi, target_n)
        phi   = rng.uniform(0.0, phi_safe_max, target_n)
        dirs  = np.stack([np.cos(theta)*np.sin(phi),
                          np.sin(theta)*np.sin(phi),
                          np.cos(phi)], axis=-1)

        norms = surface_normals(pos_sub, shape, shape_params)
        flip  = (np.sum(dirs * norms, axis=1) < 0)
        dirs[flip] = -dirs[flip]

        if gondola_direction_hard_check(dirs):
            continue

        S = np.linalg.svd(build_M_matrix(pos_sub, dirs),
                          full_matrices=False)[1]
        if S.min() < 1e-6:
            continue

        samples.append((pos_sub.copy(), dirs.copy()))

    return samples

# ── Percentile rank helper ───────────────────────────────────────────────
def _percentile_rank(value, distribution):
    """PSO's percentile within the distribution (0 = best, 100 = worst)."""
    dist = np.asarray(distribution, dtype=float)
    if dist.size == 0:
        return float('nan')
    return float(100.0 * np.mean(dist <= value))

def plot_multi_wind_comparison(label, pso_res, rand_median, rand_population,
                                r_rms, r_eng, r_sat,
                                pct_rms, pct_eng, pct_sat,
                                n_pso, n_rand):
    """Clean publication-quality 3-panel comparison."""
    import matplotlib as _mpl
    _mpl.rcParams['figure.dpi']       = 120
    _mpl.rcParams['savefig.dpi']      = 150
    _mpl.rcParams['axes.titleweight'] = 'bold'
    _mpl.rcParams['axes.labelweight'] = 'semibold'

    C_PSO   = "#1F4E79"
    C_RAND  = "#C8102E"
    C_REF   = "#000000"
    C_POP   = "#9AA3AB"

    PSO_LS  = "-"  ; PSO_MK  = "o"
    RAND_LS = "--" ; RAND_MK = "s"
    REF_LS  = ":"

    ang = pso_res['angles_deg']
    wind_str = (f"{SIM_PATH_R:.0f} m circle, {SIM_CRUISE_V:.1f} m/s cruise, "
                f"{np.linalg.norm(SIM_WIND):.1f} m/s wind from {len(ang)} angles")

    # ═══════════════════════════════════════════════════════════════════
    # Figure 1 — XY Trajectories
    # ═══════════════════════════════════════════════════════════════════
    fig1, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    n_steps = len(pso_res['runs'][0]['traj'])
    omega   = SIM_CRUISE_V / SIM_PATH_R
    t_ref   = np.arange(n_steps) * SIM_DT
    x_ref   = SIM_PATH_R * np.sin(omega * t_ref)
    y_ref   = SIM_PATH_R * (1.0 - np.cos(omega * t_ref))

    if rand_population is not None and len(rand_population) > 0:
        for r_mw in rand_population:
            for run in r_mw['runs']:
                ax.plot(run['traj'][:, 0], run['traj'][:, 1],
                         color=C_POP, alpha=0.08, linewidth=0.7, zorder=1)
        ax.plot([], [], color=C_POP, linewidth=2.0, alpha=0.5,
                 label=f"Random population  ({len(rand_population)} configs × {len(ang)} winds)")

    for i, run in enumerate(pso_res['runs']):
        lab = (f"PSO  (N={n_pso}, mean RMS={pso_res['rms_mean']:.3f} m)"
               if i == 0 else None)
        ax.plot(run['traj'][:, 0], run['traj'][:, 1],
                 color=C_PSO, linestyle=PSO_LS, linewidth=1.6,
                 alpha=0.85, label=lab, zorder=3)

    for i, run in enumerate(rand_median['runs']):
        lab = (f"Random median  (N={n_rand}, mean RMS={rand_median['rms_mean']:.3f} m)"
               if i == 0 else None)
        ax.plot(run['traj'][:, 0], run['traj'][:, 1],
                 color=C_RAND, linestyle=RAND_LS, linewidth=1.6,
                 alpha=0.85, label=lab, zorder=2)

    ax.plot(x_ref, y_ref, linestyle=REF_LS, color=C_REF, linewidth=2.2,
             label=f"Reference path  ({SIM_PATH_R:.0f} m circle)", zorder=4)

    ax.set_xlabel("x [m]", fontsize=13)
    ax.set_ylabel("y [m]", fontsize=13)
    ax.set_title(f"{label} — Trajectories  (all {len(ang)} wind directions overlaid)\n{wind_str}",
                  fontsize=13, pad=10, fontweight="bold")
    ax.axis("equal"); ax.grid(True, alpha=0.30, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=11)
    ax.legend(fontsize=10, loc="upper left", framealpha=0.95, edgecolor="black")
    ax.set_axisbelow(True)
    plt.show(); plt.close(fig1)

    # ═══════════════════════════════════════════════════════════════════
    # Figure 2 — Tracking error vs wind direction
    # ═══════════════════════════════════════════════════════════════════
    fig2, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax.plot(ang, pso_res['rms_error_per_angle'],
             color=C_PSO, linestyle=PSO_LS, marker=PSO_MK, markersize=8,
             linewidth=2.2, label="PSO", zorder=3)
    ax.plot(ang, rand_median['rms_error_per_angle'],
             color=C_RAND, linestyle=RAND_LS, marker=RAND_MK, markersize=8,
             linewidth=2.2, label="Random median", zorder=2)
    ax.axhline(pso_res['rms_mean'], color=C_PSO, linestyle=":",
                alpha=0.55, linewidth=1.4,
                label=f"PSO mean = {pso_res['rms_mean']:.3f} m")
    ax.axhline(rand_median['rms_mean'], color=C_RAND, linestyle=":",
                alpha=0.55, linewidth=1.4,
                label=f"Random mean = {rand_median['rms_mean']:.3f} m")
    ax.set_xlabel("Wind direction  [°]", fontsize=13)
    ax.set_ylabel("RMS tracking error  [m]", fontsize=13)
    ax.set_title(f"{label} — Tracking error vs wind direction\n{wind_str}",
                 fontsize=13, pad=10, fontweight="bold")
    ax.set_xticks(np.linspace(0, 360, 9))
    ax.grid(True, alpha=0.30, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=11)
    ax.legend(fontsize=10, loc="best", framealpha=0.95, edgecolor="black")
    ax.set_axisbelow(True)
    plt.show(); plt.close(fig2)

    # ═══════════════════════════════════════════════════════════════════
    # Figure 3 — Motor saturation vs wind direction
    # ═══════════════════════════════════════════════════════════════════
    fig3, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax.plot(ang, pso_res['sat_per_angle'],
             color=C_PSO, linestyle=PSO_LS, marker=PSO_MK, markersize=8,
             linewidth=2.2, label="PSO", zorder=3)
    ax.plot(ang, rand_median['sat_per_angle'],
             color=C_RAND, linestyle=RAND_LS, marker=RAND_MK, markersize=8,
             linewidth=2.2, label="Random median", zorder=2)
    ax.axhline(pso_res['sat_mean'], color=C_PSO, linestyle=":",
                alpha=0.55, linewidth=1.4,
                label=f"PSO mean = {pso_res['sat_mean']:.2f}%")
    ax.axhline(rand_median['sat_mean'], color=C_RAND, linestyle=":",
                alpha=0.55, linewidth=1.4,
                label=f"Random mean = {rand_median['sat_mean']:.2f}%")
    ax.set_xlabel("Wind direction  [°]", fontsize=13)
    ax.set_ylabel("Motor saturation  [%]", fontsize=13)
    ax.set_title(f"{label} — Motor saturation vs wind direction\n{wind_str}",
                 fontsize=13, pad=10, fontweight="bold")
    ax.set_xticks(np.linspace(0, 360, 9))
    ax.grid(True, alpha=0.30, linestyle=":")
    ax.tick_params(axis="both", which="major", labelsize=11)
    ax.legend(fontsize=10, loc="best", framealpha=0.95, edgecolor="black")
    ax.set_axisbelow(True)
    plt.show(); plt.close(fig3)

# ============================================================
# Part 7 — Main pipeline
# ============================================================
def _run_geometry(label, positions, weights, shape, shape_params, num_trials=20):
    print(f"\n{'='*60}\n  Geometry : {label}\n{'='*60}")
    dm = d_min_for_shape(shape, shape_params)
    print(f"  d_min    = {dm:.4f} normalised  "
          f"({dm*AIRSHIP_RADIUS_M*100:.1f} cm physical)")

    (all_Js, all_best, all_pos, all_dirs, all_masks,
     best_pos_list, best_dirs_list,
     all_convergences, trial_costs) = run_multiple_trials(
        positions, weights, num_trials=num_trials,
        shape=shape, shape_params=shape_params)

    # I3a — primary Pareto: authority vs effort vs count
    sel1       = all_best[:, [0, 3, 2]]
    pts1, idx1 = pareto_filter(sel1)
    print(f"\n  Pareto [J1×J4×J3] : {pts1.shape[0]} solutions from {len(all_best)} trials")
    plot_pareto_3d(pts1, labels=("J1 Authority", "J4 Effort", "J3 #Props"),
                   title=f"{label} — Pareto: Authority × Effort × Count")

    # I3b — secondary Pareto: authority vs the two symmetry metrics (J5, J6)
    sel2       = all_best[:, [0, 4, 5]]
    pts2, idx2 = pareto_filter(sel2)
    print(f"  Pareto [J1×J5×J6] : {pts2.shape[0]} solutions from {len(all_best)} trials")
    plot_pareto_3d(pts2, labels=("J1 Authority", "J5 Pos Symmetry", "J6 Dir Symmetry"),
                   title=f"{label} — Pareto: Authority × Position-Symmetry × Direction-Symmetry")

    # Best-authority solution
    best_k    = idx1[int(np.argmin(pts1[:, 0]))]
    best_pos  = best_pos_list[best_k]
    best_dirs = best_dirs_list[best_k]
    print(f"\n  >>> Best-authority solution: {len(best_pos)} props")

    plot_config_3d(best_pos, best_dirs, shape, shape_params,
                   title=f"{label} — Best Config ({len(best_pos)} props)")

    # Smart-random baseline collection (data still needed for the
    # closed-loop path-following test) — plot itself is disabled.
    print(f"\n  Computing smart-random baseline (constrained sampling) ...")
    (rc, re, rc_attempts,
     rand_best_pos, rand_best_dirs, rand_by_N, rand_all_by_N) = \
        run_random_baseline(positions, weights, shape, shape_params)
    print(f"    Collected {len(rc)} valid samples from {rc_attempts} attempts "
          f"({100.0*len(rc)/max(rc_attempts,1):.1f}% feasibility)")

    # ── FAIR closed-loop path-following validation ───────────────────────
    def _planar_auth(p, d):
        M = build_M_matrix(p, d)
        return float(np.linalg.svd(M[[0,1,5], :], full_matrices=False)[1].sum())

    N_RAND_TARGET = 20
    N_RAND_MIN    = 15

    n_angles      = 8
    fixed_seeds   = [SIM_PATH_SEED + i*1009 for i in range(n_angles)]

    pareto_with_n = sorted(
        ((k, len(best_pos_list[k]),
          _planar_auth(best_pos_list[k], best_dirs_list[k])) for k in idx1),
        key=lambda x: -x[2])

    chosen_k, matched_n, rand_samples = None, 0, []
    for k, n_k, _ in pareto_with_n:
        print(f"    trying matched N = {n_k} ...", end="")
        cand = sample_random_at_n_feasible(
            n_k, positions, shape, shape_params,
            n_target=N_RAND_TARGET, seed=42 + n_k)
        print(f"  generated {len(cand)} feasible random configs")
        if len(cand) >= N_RAND_MIN:
            chosen_k, matched_n, rand_samples = k, n_k, cand
            break

    if chosen_k is None:
        print("  [Skip] path-following: cannot generate ≥"
              f"{N_RAND_MIN} feasible random configs at any Pareto N.")
    else:
        pso_pos    = best_pos_list[chosen_k]
        pso_dirs   = best_dirs_list[chosen_k]
        n_pso      = matched_n
        pso_planar = _planar_auth(pso_pos, pso_dirs)
        pso_auth_full = float(np.linalg.svd(
            build_M_matrix(pso_pos, pso_dirs), full_matrices=False)[1].sum())

        print(f"\n  Path-following  (MATCHED N = {n_pso} for both sides)")
        print(f"    PSO Pareto config: planar-σ-sum={pso_planar:.3f}, "
              f"6-DOF σ-sum={pso_auth_full:.3f}")
        print(f"    Random population: {len(rand_samples)} feasible configs "
              f"(same gondola gates, same min-spacing, rank-6)")

        print(f"    Auto-tuning PD gains (same procedure) ...")
        pso_gains = tune_pd_gains(pso_pos, pso_dirs)
        rand_gains_list = [tune_pd_gains(rp, rd) for rp, rd in rand_samples]
        def _gfmt(g):
            return (f"Kp_x={g[0]:.2f}  Kd_x={g[1]:.2f}  "
                    f"Kp_y={g[2]:.2f}  Kd_y={g[3]:.2f}  "
                    f"Kp_yaw={g[4]:.2f}  Kd_yaw={g[5]:.2f}")
        print(f"      PSO gains  : {_gfmt(pso_gains)}")
        median_gain_idx = len(rand_gains_list) // 2
        sample_g        = rand_gains_list[median_gain_idx]
        print(f"      sample rand: {_gfmt(sample_g)}  (one example)")

        print(f"    Running {n_angles}-direction path test for PSO + "
              f"{len(rand_samples)} random ...")
        pso_res = simulate_multi_wind(pso_pos, pso_dirs, n_angles=n_angles,
                                       gains=pso_gains,
                                       per_angle_seeds=fixed_seeds)
        rand_results = [
            simulate_multi_wind(rp, rd, n_angles=n_angles,
                                gains=g, per_angle_seeds=fixed_seeds)
            for (rp, rd), g in zip(rand_samples, rand_gains_list)
        ]

        r_rms = np.array([r['rms_mean']    for r in rand_results])
        r_eng = np.array([r['energy_mean'] for r in rand_results])
        r_sat = np.array([r['sat_mean']    for r in rand_results])

        med_idx     = int(np.argsort(r_rms)[len(r_rms)//2])
        rand_median = rand_results[med_idx]

        pct_rms = _percentile_rank(pso_res['rms_mean'],    r_rms)
        pct_eng = _percentile_rank(pso_res['energy_mean'], r_eng)
        pct_sat = _percentile_rank(pso_res['sat_mean'],    r_sat)

        print(f"\n  ╔════════════════════════════════════════════════════════════════╗")
        print(f"  ║  FAIR Path-following verdict  [{label}]                            ║")
        print(f"  ║  Matched N = {n_pso} (both sides, no fallback)                     ║")
        print(f"  ║  Random population = {len(rand_samples)} feasible configs "
              f"× {n_angles} wind directions    ║")
        print(f"  ║  PD gains auto-tuned per config (identical procedure)            ║")
        print(f"  ║  Wind, gust seeds, controller, dynamics — all identical          ║")
        print(f"  ╚════════════════════════════════════════════════════════════════╝")
        print(f"    {'Metric':<22} {'PSO':>10}  {'Rand med':>10}  "
              f"{'Rand IQR':>22}  {'PSO pct':>8}")
        print(f"    {'-'*22} {'-'*10}  {'-'*10}  {'-'*22}  {'-'*8}")
        print(f"    {'RMS error [m]':<22} "
              f"{pso_res['rms_mean']:>10.3f}  "
              f"{np.median(r_rms):>10.3f}  "
              f"[{np.percentile(r_rms,25):>5.3f}, {np.percentile(r_rms,75):>5.3f}]    "
              f"{pct_rms:>7.1f}%")
        print(f"    {'Energy [N²·s]':<22} "
              f"{pso_res['energy_mean']:>10.1f}  "
              f"{np.median(r_eng):>10.1f}  "
              f"[{np.percentile(r_eng,25):>5.1f}, {np.percentile(r_eng,75):>5.1f}]    "
              f"{pct_eng:>7.1f}%")
        print(f"    {'Saturation [%]':<22} "
              f"{pso_res['sat_mean']:>10.2f}  "
              f"{np.median(r_sat):>10.2f}  "
              f"[{np.percentile(r_sat,25):>5.2f}, {np.percentile(r_sat,75):>5.2f}]    "
              f"{pct_sat:>7.1f}%")
        print(f"\n    PSO percentile interpretation: 0% = PSO best, 100% = PSO worst.")
        print(f"    A PSO at <50% means it BEATS the random median.")
        if max(pct_rms, pct_eng, pct_sat) < 50:
            print(f"    >>> PSO beats random median on ALL three metrics.")
        elif min(pct_rms, pct_eng, pct_sat) < 50:
            print(f"    >>> PSO beats random median on some metrics.")
        else:
            print(f"    >>> PSO is at or above random median — placement quality"
                  f" not detectable at this N.")

        plot_multi_wind_comparison(label, pso_res, rand_median,
                                     rand_population=rand_results,
                                     r_rms=r_rms, r_eng=r_eng, r_sat=r_sat,
                                     pct_rms=pct_rms, pct_eng=pct_eng,
                                     pct_sat=pct_sat,
                                     n_pso=n_pso, n_rand=matched_n)

    print(f"\n  Trial costs:  mean={trial_costs.mean():.4f}  "
          f"std={trial_costs.std():.4f}  "
          f"cv={trial_costs.std()/abs(trial_costs.mean())*100:.1f}%")

    return (all_Js, all_best, all_pos, all_dirs, all_masks,
            idx1, trial_costs, best_pos, best_dirs)

def main():
    _J_LABELS = ("J1","J2","J3","J4","J5","J6","J7","J8","J9")
    _on  = [j for j, m in zip(_J_LABELS, ACTIVE_OBJ_MASK) if m]
    _off = [j for j, m in zip(_J_LABELS, ACTIVE_OBJ_MASK) if not m]
    print(f"  [Cost configuration]  active = {_on}   excluded = {_off}")

    K = 30
    # F7: J1 (authority) raised to 2x
    weights    = np.ones(9)
    weights[0] = 2.0
    weights   /= weights.sum()
    costs_per_shape = {}

    # Sphere (STARLab reference)
    pos_sphere = fibonacci_sphere(K, radius=1.0)
    *_, tc_sphere, bp_sphere, bd_sphere = _run_geometry(
        "Sphere", pos_sphere, weights, "sphere", None, num_trials=20)
    costs_per_shape["Sphere"] = tc_sphere

    # Cigar (prolate)
    a, b, c = 1.0, 1.0, 2.5
    pos_cig = fibonacci_ellipsoid(K, a, b, c)
    *_, tc_cig, bp_cig, bd_cig = _run_geometry(
        "Cigar (ellipsoid)", pos_cig, weights, "ellipsoid", (a,b,c), num_trials=20)
    costs_per_shape["Cigar"] = tc_cig

    # Flat disk (oblate)
    a2, b2, c2 = 2.0, 2.0, 0.5
    pos_disk   = fibonacci_ellipsoid(K, a2, b2, c2)
    *_, tc_disk, bp_disk, bd_disk = _run_geometry(
        "Flat Disk (oblate)", pos_disk, weights, "ellipsoid", (a2,b2,c2), num_trials=20)
    costs_per_shape["Flat Disk"] = tc_disk

# ============================================================
# Part 8 — PyVista
# ============================================================
def _ellipsoid_r_at_z(z, a, c):
    return a * np.sqrt(max(1.0 - (z/c)**2, 0.0))

def make_airship_pyvista(shape, shape_params):
    if not PYVISTA_OK:
        return None, None, []
    a = b = c = 1.0
    if shape != "sphere":
        a, b, c = shape_params
    hull = pv.Sphere(radius=1.0, theta_resolution=80, phi_resolution=60)
    hull.scale([a, b, c], inplace=True)
    gond_r = max(0.09*(a+b)*0.5, 0.06)
    gond_h = max(0.28*c*0.5, 0.22)
    # F3: gondola placed at -Z to match GONDOLA_DOWN = [0, 0, -1]
    gondola = pv.Cylinder(center=(0., 0., -(c+0.22)), direction=(0., 0., 1.),
                          radius=gond_r, height=gond_h, resolution=20, capping=True)
    z_att   = -c*0.78
    r_hull  = _ellipsoid_r_at_z(z_att, a, c)
    fin_r   = 0.55*b
    fin_axial = max(0.30*c, 0.20)
    fin_t   = max(0.042*a, 0.025)
    fins = []
    for ang in [0., 90., 180., 270.]:
        fin = pv.Box(bounds=(-fin_t/2, fin_t/2, r_hull, r_hull+fin_r,
                              z_att-fin_axial*0.35, z_att+fin_axial*0.65))
        fin.rotate_z(ang, inplace=True); fins.append(fin)
    return hull, gondola, fins

def _dof_authority(positions, directions):
    if len(positions) == 0:
        return np.zeros(6)
    M    = build_M_matrix(np.asarray(positions), np.asarray(directions))
    Minv = np.linalg.pinv(M, rcond=1e-9)
    return np.array([1./max(np.linalg.norm(Minv @ np.eye(6)[:, d]), 1e-9)
                     for d in range(6)])

def _add_airship_to_plotter(pl, hull, gondola, fins, positions, directions,
                            shape, shape_params):
    a = b = c = 1.0
    if shape != "sphere":
        a, b, c = shape_params
    mx      = max(a, b, c)
    prop_r  = 0.12*max(a, b)
    hub_r   = 0.045*max(a, b)
    arr_mag = 0.42*mx
    pl.add_mesh(hull, color="lightblue", opacity=0.38, smooth_shading=True,
                show_edges=False)
    pl.add_mesh(gondola, color="#505050", opacity=1.0, smooth_shading=True)
    for fin in fins:
        pl.add_mesh(fin, color="#505050", opacity=1.0, smooth_shading=True)
    pos  = np.asarray(positions, dtype=float)
    dirs = np.asarray(directions, dtype=float)
    for p, d in zip(pos, dirs):
        pl.add_mesh(pv.Disc(center=p.tolist(), normal=d.tolist(),
                            inner=0., outer=prop_r, r_res=2, c_res=32),
                    color="crimson", opacity=0.82)
        pl.add_mesh(pv.Sphere(radius=hub_r, center=p.tolist()), color="#1a1a1a")
    if len(pos) > 0:
        pl.add_arrows(pos, dirs, mag=arr_mag, color="darkorange")

def compare_airships_pyvista(configs, off_screen=False, save_path=None):
    if not PYVISTA_OK:
        print("[Skip] PyVista not available."); return
    n  = len(configs)
    pl = pv.Plotter(shape=(1, n), window_size=(720*n, 760), off_screen=off_screen)
    for i, cfg in enumerate(configs):
        pl.subplot(0, i); pl.set_background("white")
        _add_airship_to_plotter(pl, cfg["hull"], cfg["gondola"], cfg["fins"],
                                cfg["positions"], cfg["directions"],
                                cfg["shape"], cfg["shape_params"])
        auth = _dof_authority(np.asarray(cfg["positions"]),
                              np.asarray(cfg["directions"]))
        short = ["Fx","Fy","Fz","Roll","Pitch","Yaw"]
        row = "   ".join(f"{l}:{v:.2f}" for l, v in zip(short, auth))
        pl.add_text(f"{cfg['label']}  |  {len(cfg['positions'])} props\n{row}",
                    position="upper_left", font_size=9, color="black")
    pl.reset_camera()
    if save_path:
        pl.screenshot(save_path, transparent_background=False)
    if not off_screen:
        pl.show()

if __name__ == "__main__":
    main()
