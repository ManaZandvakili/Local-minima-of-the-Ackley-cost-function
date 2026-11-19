%RUN_AIRSHIP_SIMULINK  Replay the Python reference flight through the Simulink model.
%   Run from the matlab/ folder:  run_airship_simulink
%
%   Builds airship_sim.slx, feeds it the thrust commands and wind recorded in
%   Python (data/reference_flight.mat), and compares the states with Python's.
%   Needs MATLAB with Simulink (and a C compiler for the MATLAB Function block,
%   which MATLAB normally has).

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here);
R = load(fullfile(here, 'data', 'reference_flight.mat'));
ref = R.ref;
mdl = build_airship_simulink(ref.layout);
P = evalin('base', 'P');

K = size(ref.cmd, 1);
n = size(ref.wind, 1);
% times are set half a step early so the held value switches exactly on the sample
tc = (0:K - 1).' * P.n_sub * P.dt - P.dt / 2;  tc(1) = 0;
tw = (0:n - 1).' * P.dt - P.dt / 2;            tw(1) = 0;
assignin('base', 'cmd_data', [tc, ref.cmd]);
assignin('base', 'wind_data', [tw, ref.wind]);
assignin('base', 'x0', ref.x0(:));

out = sim(mdl, 'StopTime', num2str(K * P.n_sub * P.dt));
x_log = out.get('x_log');
if ndims(x_log) == 3
    x_log = squeeze(x_log).';                    % [nx x 1 x T] -> [T x nx]
end
idx = (1:K).' * P.n_sub + 1;                     % rows at t = 0.1, 0.2, ..., 60 s
X = x_log(idx, :);
dpos = max(max(abs(X(:, 1:3) - ref.x_python(:, 1:3))));
du = max(max(abs(X(:, 14:end) - ref.x_python(:, 14:end))));
fprintf('Simulink vs Python: largest difference position %.2e m, thrust %.2e N\n', dpos, du);
if dpos < 1e-6
    fprintf('PASS: the Simulink model gives the same flight as Python\n');
else
    fprintf('FAIL: check the time alignment of cmd_data and wind_data\n');
end

t = ref.t_x(:);
figure('Name', 'Reference flight');
subplot(2, 1, 1);
plot(t, X(:, 1:3)); hold on; plot(t, ref.x_python(:, 1:3), 'k:');
ylabel('Position [m]'); legend('x', 'y', 'z', 'Python'); grid on;
title('Simulink (colour) vs Python (dotted)');
subplot(2, 1, 2);
plot(t, X(:, 14:end));
ylabel('Thrust per propeller [N]'); xlabel('Time [s]'); grid on;
