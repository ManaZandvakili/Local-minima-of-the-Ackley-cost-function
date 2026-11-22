# Propeller layout optimization and learning-based control of a spherical airship

Code for the journal paper (extension of the conference paper by Zandvakili et al.).
It contains everything needed to reproduce the results:

1. **Layout optimizer**: particle swarm optimization of where the propellers go on the sphere,
   with a bounded-thrust capacity test, and 20 random layouts per optimized layout for comparison.
2. **Flight simulator**: 6-degree-of-freedom rigid-body model of the 3.15 m sphere (added mass,
   drag on the air-relative velocity, buoyancy and weight, motor lag), wind with turbulence and gusts.
3. **Controllers**: thrust allocation (bounded weighted least squares), a cascaded PID tuned per
   layout, and two reinforcement-learning controllers trained with PPO: plain (the network gives
   the whole command) and residual (the network corrects the PID).
4. **Tests**: station keeping in gusty wind and the conference paper's 50 m circle.
5. **MATLAB/Simulink** version of the full closed loop (`matlab/`).

## Install

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
python -m pytest -q            # 67 tests, about 1 minute
```

## Reproduce the results

All commands are run from this folder. Each step writes a new time-stamped folder under
`results/`; put those folder names into the next command. Propellers are reversible
(`reverse_thrust: true` in `configs/sphere.yaml`). Times are for 4 CPU cores.

```bash
# 1. reference layouts (conference paper, hand-designed)              ~1 min
python scripts/export_reference_layouts.py
python scripts/make_hand_layouts.py

# 2. layout optimizer: 5 seeds x 20 runs, then the run limited to 11 propellers   ~1 h each
python scripts/run_optimizer.py --workers 4
python scripts/run_optimizer.py --workers 4 --max-props 11

# 3. tune the PID of every layout                                    ~1 h
python scripts/tune_pid.py --run results/layouts/phase1_<run> \
    --capped-run results/layouts/phase1_maxN11_<run> --workers 4

# 4. train the learning controllers (5 million steps each)           ~75 min each
python scripts/train_rl.py --seed 0 --tag _plain --layout-runs phase1_<run> phase1_maxN11_<run>
python scripts/train_rl.py --seed 0 --residual --pid-run results/pid/<pid run> \
    --layout-runs phase1_<run> phase1_maxN11_<run>

# 5. station keeping in gusty wind (3 m/s, three directions)          ~20 min per run
python scripts/evaluate_rl.py --rl results/rl/<residual run> --pid results/pid/<pid run> --workers 4 --no-videos
python scripts/evaluate_rl.py --rl results/rl/<plain run> --pid results/pid/<pid run> --workers 4 --no-videos

# 6. the 50 m circle (1.5 and 3 m/s wind, three directions)           ~25 min
python scripts/path_following.py --case reversible --pid results/pid/<pid run> \
    --residual results/rl/<residual run> --plain results/rl/<plain run>

# 7. the main figure (optimized vs random 8-propeller layout, three controllers, two tests)
python scripts/main_result_figure.py --residual results/rl/<residual run> --plain results/rl/<plain run> \
    --circle results/path/<run>
```

Reports: `results/runs/<run>/report.md` (optimizer), `results/pid/<run>/report.md` (PID),
`results/rl/<run>/eval/report.md` (station keeping), `results/path/<run>/report.md` (circle).

### Shortcut: the runs used in the paper

`results/` already holds the layouts, tuned gains and trained networks used in the paper,
so steps 2 to 4 can be skipped. For example:

```bash
python scripts/evaluate_rl.py --rl results/rl/20261117_155720_bidir_residual_seed0 \
    --pid results/pid/20261117_153423 --workers 4 --no-videos
python scripts/path_following.py --case reversible      # uses these runs by default
```

Results do not depend on the computer, except for tiny floating-point differences:
every random number (optimizer, random layouts, wind, training) comes from a fixed seed.

## MATLAB and Simulink

`matlab/` holds the same simulator, thrust allocation, PID and trained networks.
See `matlab/README.md`. In MATLAB, from `matlab/`:

```matlab
validate_closed_loop      % checks against the Python flights (all PASS)
compare_controllers       % optimized vs random layout, three controllers, two tests
run_airship_closed_loop   % builds and runs the Simulink model airship_closed_loop.slx
```

## Contents

| Folder / file | What it is |
|---|---|
| `airship/` | the Python package: geometry, force matrix, feasibility, optimizer, simulator, wind, allocation, controllers, learning environment |
| `configs/` | airship, optimizer, PID tuning and training settings (every physical value states its source) |
| `scripts/` | the programs in the steps above |
| `tests/` | automatic checks, including agreement with the conference-paper code |
| `reference/Final.py` | the conference-paper code (used by the tests and for the paper's layouts) |
| `matlab/` | MATLAB/Simulink version |
| `results/` | layouts, tuned gains and trained networks used in the paper |
