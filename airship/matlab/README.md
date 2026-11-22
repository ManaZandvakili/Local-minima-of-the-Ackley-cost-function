# Sphere airship simulator in MATLAB and Simulink

The same simulator as the Python code (`airship/airship/dynamics.py` and
`wind.py`), ported line by line. Same equations, same numbers, same time step
(0.01 s, 4th-order Runge-Kutta). Replaying a recorded Python flight gives the same
states to within 1e-14 m.

## Full closed loop (airship + wind + target + controller + thrust sharing)

The complete flight simulator with the three controllers of the journal paper:
the classic controller (PID), residual learning and plain learning (the trained
networks, run directly in MATLAB). In MATLAB, go to this folder and run:

```matlab
validate_closed_loop          % plain MATLAB: 6 checks against Python, all should say PASS
compare_controllers           % plain MATLAB: optimized vs typical random 8-propeller layout, 3 controllers, table + chart
run_airship_closed_loop       % needs Simulink: builds airship_closed_loop.slx and flies it
run_airship_closed_loop('random8_reversible_typical', 'pid', 'station')   % any layout, controller, test
```

Controllers: `'pid'`, `'residual'`, `'plain'`. Tests: `'station'` (hold position for
60 s in gusty 3 m/s wind, heading step of 45 degrees at 30 s) or `'circle'` (the
paper's 50 m circle at 0.5 m/s, one lap, 3 m/s gusty wind). The exact Python wind is
used, so the errors are the same as in the Python results. For a new wind:
`run_closed_loop(closed_loop_setup('knee8_reversible_capped', 'residual', 'circle', 42))`.

The Simulink model (every block is a MATLAB Function block calling the file in brackets):

```
[Clock] -> [Target (reference_signal.m)] ----------------------\
       \-> [Wind (environment) (wind_lookup.m)] -------------\  |
                                                             v  v
 [State x] -> [Controller (controller_step.m)] -> [Thrust sharing (allocation_step.m)] -> [Airship (plant) (airship_plant_step.m)] -+
     ^              ^ memory (delay)        ^ last result (delay)                                                                 |
     +----------------------------------------------------------------------------------------------------------------------------+
```

It runs every 0.1 s (the controller rate); the airship block steps the physics 10 times
at 0.01 s inside, exactly as the Python simulator. Logs: `x_log`, `u_log`, `ref_log`, `tau_log`.

| File | What it does |
|---|---|
| `closed_loop_setup.m` | Loads a layout, its tuned PID gains, the trained networks, the target and the wind into one struct |
| `reference_signal.m` | Target: fixed point with a heading step, or the paper's circle |
| `controller_step.m` | One 0.1 s step of the classic, residual or plain controller |
| `pid_step.m` | The cascaded PID (position, roll and pitch, heading) with anti-windup |
| `make_observation.m`, `policy_action.m` | What the network sees, and the network itself (2 x 256 tanh) |
| `allocation_step.m` | Thrust sharing: bounded weighted least squares (active-set method) |
| `airship_plant_step.m` | The airship over 0.1 s (10 Runge-Kutta steps) |
| `wind_lookup.m` | Wind of the current step from the wind table |
| `run_closed_loop.m` | The same loop in plain MATLAB, with the scores |
| `validate_closed_loop.m` | Checks against the Python flights |
| `compare_controllers.m` | The main comparison: optimized vs random, three controllers, two tests |
| `build_airship_closed_loop.m`, `run_airship_closed_loop.m` | Build and run the Simulink model |
| `data/controllers/` | Tuned PID gains of each layout; the residual and plain networks (reversible and one-way) |
| `data/validation/` | Python flights and their wind, for the checks |

## Replay only (the airship and wind, with recorded thrust commands)

In MATLAB, go to this folder and run:

```matlab
validate_against_python   % plain MATLAB: 4 checks, all should say PASS
plot_reference_flight     % plain MATLAB: flies the reference flight and plots it (no Simulink needed)
run_airship_simulink      % needs Simulink: builds airship_sim.slx, replays the same flight, plots it
```

To build the Simulink model for another layout:

```matlab
build_airship_simulink('knee9_oneway_capped')
```

## Files

| File | What it does |
|---|---|
| `airship_derivative.m` | The equations of motion: thrust, buoyancy and weight, drag, damping, added mass, Coriolis terms, motor lag |
| `airship_rk4_step.m` | One 0.01 s step: clip the command, Runge-Kutta step, re-normalise the attitude, clip the motor thrusts |
| `quat_to_rot.m` | Attitude quaternion to rotation matrix |
| `airship_load.m` | Loads the physical model and one layout into the parameter struct `P` |
| `plot_reference_flight.m` | Flies the recorded gusty flight in plain MATLAB and plots position, path, wind, thrusts and heading |
| `run_airship.m` | Flies a list of thrust commands through the simulator in plain MATLAB |
| `wind_model.m` | Mean wind + turbulence + gusts, same model as Python (MATLAB's own random numbers) |
| `validate_against_python.m` | Replay of a Python flight, pendulum period, energy and wind checks |
| `build_airship_simulink.m` | Builds `airship_sim.slx` |
| `run_airship_simulink.m` | Replays the Python flight through the Simulink model and compares |
| `data/airship_params.mat` | Mass matrices, weight, buoyancy, centre of gravity, drag, damping, motor time constant |
| `data/layouts/*.mat` | Nine layouts: the one-way and reversible knees, the capped knees, the typical random 8- and 9-propeller layouts, the paper layout both ways, STARLab |
| `data/reference_flight.mat` | A 60 s gusty flight from Python: commands, wind and states |

The `.mat` files are written by `python scripts/export_matlab.py` (from the
`airship/` folder). Run it again after changing `configs/sphere.yaml`, so MATLAB
uses the new numbers.

## The state vector

| Index | Meaning | Unit |
|---|---|---|
| `x(1:3)` | position, world frame (z up) | m |
| `x(4:7)` | attitude quaternion `[w x y z]`, body to world | - |
| `x(8:10)` | velocity, body frame | m/s |
| `x(11:13)` | angular rate, body frame | rad/s |
| `x(14:end)` | actual thrust of each propeller | N |

Body frame: origin at the hull centre, z up, gondola below (-z).

## The Simulink model

```
[Thrust command] --\
[Wind] -------------> [Airship step] --> [State x (unit delay)] --+--> [Log x]
                   ^                                               +--> [Position scope]
                   +-----------------------------------------------+
```

A fixed-step discrete loop at 0.01 s. "Airship step" is a MATLAB Function block
that calls `airship_rk4_step.m`. The inputs come from the workspace (`cmd_data`,
`wind_data`, `x0`, `P`). To fly your own controller, replace "Thrust command"
with your controller, fed from the "State x" signal, and give it a sample time of
0.1 s.

## What is and is not included

- Included: the airship physics, the motors, the wind model, nine layouts (including the
  typical random 8- and 9-propeller layouts), the thrust sharing, the tuned PID of every
  layout and the trained residual and plain networks.
- Not included: training (the networks are trained in Python and exported), motor
  failures, and the layout optimizer.
- The plain-MATLAB files were tested in GNU Octave 8.4: `validate_closed_loop` passes
  all 6 checks (largest difference from Python 1e-6 m over a full 650 s lap). The
  Simulink builders could not be run there (Simulink needs MATLAB); the blocks call the
  same tested files, in the same order as `run_closed_loop.m`.
