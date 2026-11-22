%PLOT_REFERENCE_FLIGHT  Fly the recorded reference flight in plain MATLAB and plot it.
%   Run from the matlab/ folder:  plot_reference_flight
%   No Simulink needed. Same flight as run_airship_simulink: the reversible
%   14-propeller knee layout in 3 m/s gusty wind, held by the tuned classic
%   controller (its thrust commands were recorded in Python), with a heading
%   change to 45 degrees at 30 s.

clear;
here = fileparts(mfilename('fullpath'));
addpath(here);
R = load(fullfile(here, 'data', 'reference_flight.mat'));
ref = R.ref;
P = airship_load(ref.layout);
[X, t] = run_airship(P, ref.x0, ref.cmd, ref.wind);

pos = X(:, 1:3);
q = X(:, 4:7);
yaw = atan2(2 * (q(:, 1) .* q(:, 4) + q(:, 2) .* q(:, 3)), 1 - 2 * (q(:, 3).^2 + q(:, 4).^2));
err = sqrt(sum(pos.^2, 2));
tw = (0:size(ref.wind, 1) - 1).' * P.dt;
fprintf('Position error: average (RMS) %.2f m, largest %.2f m\n', sqrt(mean(err.^2)), max(err));
fprintf('Difference from the Python flight: %.1e m\n', max(max(abs(pos - ref.x_python(:, 1:3)))));

figure('Name', 'Reference flight (plain MATLAB)', 'Color', 'w');
subplot(2, 2, 1);
plot(t, pos, 'LineWidth', 1.2); hold on; plot(t, ref.x_python(:, 1:3), 'k:');
ylabel('Position [m]'); legend('x', 'y', 'z', 'Python', 'Location', 'best'); grid on;
title('Position (target is 0, 0, 0)');
subplot(2, 2, 2);
plot(pos(:, 1), pos(:, 2), 'LineWidth', 1.2); hold on; plot(0, 0, 'k+', 'MarkerSize', 12);
plot(pos(1, 1), pos(1, 2), 'go'); axis equal; grid on;
xlabel('x [m]'); ylabel('y [m]'); title('Path seen from above (green = start, + = target)');
subplot(2, 2, 3);
plot(tw, ref.wind(:, 1:2)); grid on; xlabel('Time [s]'); ylabel('Wind [m/s]');
legend('wind x', 'wind y', 'Location', 'best'); title('Wind (mean 3 m/s + turbulence + gusts)');
subplot(2, 2, 4);
plot(t, X(:, 14:end)); hold on; yyaxis right; plot(t, rad2deg(yaw), 'k', 'LineWidth', 1.5);
ylabel('Heading [deg]'); yyaxis left; ylabel('Thrust per propeller [N]');
xlabel('Time [s]'); grid on; title('Propeller thrusts and heading');
