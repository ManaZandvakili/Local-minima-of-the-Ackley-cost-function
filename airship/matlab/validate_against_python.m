%VALIDATE_AGAINST_PYTHON  Check that the MATLAB simulator matches the Python one.
%   Run from the matlab/ folder:  validate_against_python
%
%   1. Replay: the 60 s gusty flight recorded in Python (same thrust commands,
%      same wind) is flown again here; the states must match.
%   2. Pendulum: released at 10 degrees roll in still air with no thrust, the
%      airship must swing with the period the pendulum formula predicts.
%   3. Energy: with drag, damping and wind off, total energy must stay constant.
%   4. Wind: the turbulence must have the requested spread.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
ok = true;

%% 1. Replay the Python reference flight
R = load(fullfile(here, 'data', 'reference_flight.mat'));
ref = R.ref;
P = airship_load(ref.layout);
tic;
X = run_airship(P, ref.x0, ref.cmd, ref.wind);
el = toc;
dpos = max(max(abs(X(:, 1:3) - ref.x_python(:, 1:3))));
dq = max(max(abs(X(:, 4:7) - ref.x_python(:, 4:7))));
dvel = max(max(abs(X(:, 8:13) - ref.x_python(:, 8:13))));
du = max(max(abs(X(:, 14:end) - ref.x_python(:, 14:end))));
fprintf('1. Replay of the Python flight (%d s, %s, %.1f s to run)\n', ...
        round(size(ref.cmd, 1) * 0.1), ref.layout, el);
fprintf('   largest difference: position %.2e m, attitude %.2e, velocity %.2e, thrust %.2e N\n', ...
        dpos, dq, dvel, du);
if dpos < 1e-6 && du < 1e-6
    fprintf('   PASS: MATLAB and Python give the same flight\n');
else
    fprintf('   FAIL\n'); ok = false;
end

%% 2. Pendulum period
P0 = airship_load('knee14_reversible');
roll0 = deg2rad(10);
x = zeros(13 + P0.n_props, 1);
x(4:7) = [cos(roll0 / 2); sin(roll0 / 2); 0; 0];
n = 3000;                                        % 30 s
roll = zeros(n, 1);
for k = 1:n
    x = airship_rk4_step(x, zeros(P0.n_props, 1), [0; 0; 0], P0);
    q = x(4:7);
    roll(k) = atan2(2 * (q(1) * q(2) + q(3) * q(4)), 1 - 2 * (q(2)^2 + q(3)^2));
end
t = (1:n).' * P0.dt;
up = find(roll(1:end - 1) < 0 & roll(2:end) >= 0);       % upward zero crossings
tc = t(up) - roll(up) ./ (roll(up + 1) - roll(up)) * P0.dt;
T_sim = mean(diff(tc));
% small-angle pendulum about the hull centre: inertia incl. added inertia, restoring m*g*h
h = -P0.r_g(3);
Ixx = P0.inertia_kgm2(1, 1) + P0.added_inertia_kgm2(1, 1);
m_eff = P0.mass_kg + P0.added_mass_kg;
% the hull centre also moves sideways; the exact linear period comes from the 2x2 sway-roll system
Mss = [m_eff, P0.mass_kg * h; P0.mass_kg * h, Ixx];
Kss = [0 0; 0 P0.mass_kg * 9.81 * h];
T_theory = 2 * pi / sqrt(max(eig(Mss \ Kss)));
fprintf('2. Pendulum: period %.3f s in the simulator, %.3f s from theory\n', T_sim, T_theory);
if abs(T_sim - T_theory) / T_theory < 0.02
    fprintf('   PASS\n');
else
    fprintf('   FAIL\n'); ok = false;
end

%% 3. Energy with drag, damping and wind turned off
Pe = P0; Pe.k_drag = 0; Pe.d_lin = [0; 0; 0]; Pe.d_quad = [0; 0; 0];
x = zeros(13 + Pe.n_props, 1);
x(4:7) = [cos(roll0 / 2); sin(roll0 / 2); 0; 0];
x(8:13) = [0.3; -0.2; 0.1; 0.05; -0.04; 0.1];
E = @(x) 0.5 * x(8:13).' * (Pe.Mrb + Pe.Ma) * x(8:13) + ...
         Pe.weight * (x(3) + [0 0 1] * (quat_to_rot(x(4:7)) * Pe.r_g)) - Pe.buoy * x(3);
E0 = E(x);
for k = 1:2000
    x = airship_rk4_step(x, zeros(Pe.n_props, 1), [0; 0; 0], Pe);
end
dE = abs(E(x) - E0);
fprintf('3. Energy after 20 s without drag: changed by %.1e J\n', dE);
if dE < 1e-6
    fprintf('   PASS\n');
else
    fprintf('   FAIL\n'); ok = false;
end

%% 4. Wind statistics
cfg = struct('mean_speed_m_s', 3, 'mean_dir_deg', 120, 'turbulence_sigma_m_s', [0.5 0.5 0.25]);
W = wind_model(cfg, 1, 200000, 0.01);
s = std(W);
fprintf('4. Wind: mean (%.2f, %.2f) m/s, turbulence spread (%.2f, %.2f, %.2f) m/s\n', ...
        mean(W(:, 1)), mean(W(:, 2)), s(1), s(2), s(3));
fprintf('   requested mean (%.2f, %.2f), spread (0.50, 0.50, 0.25)\n', 3 * cosd(120), 3 * sind(120));
if all(abs(s - [0.5 0.5 0.25]) < 0.05)
    fprintf('   PASS\n');
else
    fprintf('   FAIL\n'); ok = false;
end

if ok
    fprintf('\nAll checks passed.\n');
else
    fprintf('\nSome checks FAILED.\n');
end
