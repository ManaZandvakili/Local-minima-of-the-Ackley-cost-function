# Sphere airship simulator in MATLAB and Simulink

The same simulator as the Python code (`airship/airship/dynamics.py` and
`wind.py`), ported line by line. Same equations, same numbers, same time step
(0.01 s, 4th-order Runge-Kutta). Replaying a recorded Python flight gives the same
states to within 1e-14 m.

## Quick start

In MATLAB, go to this folder and run:

```matlab
validate_against_python   % plain MATLAB: 4 checks, all should say PASS
run_airship_simulink      % builds airship_sim.slx, replays the same flight in Simulink, plots it
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
| `run_airship.m` | Flies a list of thrust commands through the simulator in plain MATLAB |
| `wind_model.m` | Mean wind + turbulence + gusts, same model as Python (MATLAB's own random numbers) |
| `validate_against_python.m` | Replay of a Python flight, pendulum period, energy and wind checks |
| `build_airship_simulink.m` | Builds `airship_sim.slx` |
| `run_airship_simulink.m` | Replays the Python flight through the Simulink model and compares |
| `data/airship_params.mat` | Mass matrices, weight, buoyancy, centre of gravity, drag, damping, motor time constant |
| `data/layouts/*.mat` | Seven layouts: the one-way and reversible knees, the capped knees, the paper layout both ways, STARLab |
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

- Included: the airship physics, the motors, the wind model, seven layouts.
- Not yet ported: the thrust-sharing solver and the tuned PID controller (the
  reference flight carries their commands, so the replay does not need them).
  They can be ported the same way if you want to close the loop in Simulink.
- The plain-MATLAB files were tested in GNU Octave 8.4 (all checks pass).
  The Simulink builder could not be run there (Simulink needs MATLAB); its timing
  was checked by imitating Simulink's sample-and-hold in Octave, which gives the
  same flight as Python.
