function out = run_airship_closed_loop(layout_name, controller, scenario)
%RUN_AIRSHIP_CLOSED_LOOP  Build and run the full closed-loop Simulink model, then plot it.
%   Run from the matlab/ folder. Needs Simulink.
%   run_airship_closed_loop                                   % optimized 8 propellers, residual, circle
%   run_airship_closed_loop('knee8_reversible_capped', 'pid', 'station')
%   run_airship_closed_loop('random8_reversible_typical', 'plain', 'circle')
%
%   controller: 'pid' (classic), 'residual' (residual learning), 'plain' (plain learning)
%   scenario:   'station' (hold position, 60 s) or 'circle' (the paper's 50 m circle, 650 s)
%   Uses the exact Python wind, so the error printed here equals the Python result.
%   The model file airship_closed_loop.slx stays open: you can look inside every block.

if nargin < 1 || isempty(layout_name), layout_name = 'knee8_reversible_capped'; end
if nargin < 2 || isempty(controller), controller = 'residual'; end
if nargin < 3 || isempty(scenario), scenario = 'circle'; end
here = fileparts(mfilename('fullpath'));
addpath(here);
S = closed_loop_setup(layout_name, controller, scenario);
mdl = build_airship_closed_loop(S);
open_system(mdl);
fprintf('Simulating %.0f s ...\n', S.K * S.dt);
so = sim(mdl);
X = getlog(so, 'x_log');
REF = getlog(so, 'ref_log');
U = getlog(so, 'u_log');
n = size(X, 1);
t = (1:n).' * S.dt;
nxt = zeros(n, 3);
for k = 1:n
    r = reference_signal(k * S.dt, S.RP);
    nxt(k, :) = r(1:3).';
end
err = sqrt(sum((X(:, 1:3) - nxt) .^ 2, 2));
out = struct('X', X, 'U', U, 'REF', REF, 'T', t, 'err', err, 'rms_err_m', sqrt(mean(err .^ 2)), ...
             'max_err_m', max(err), 'effort_N2s', sum(sum(X(:, 14:end) .^ 2)) * S.dt);
names = struct('pid', 'classic controller (PID)', 'residual', 'residual learning', 'plain', 'plain learning');
fprintf('%s, %s, %s: RMS error %.3f m, largest %.2f m, effort %.0f N^2 s\n', layout_name, ...
        names.(controller), scenario, out.rms_err_m, out.max_err_m, out.effort_N2s);

figure('Name', sprintf('%s - %s - %s', layout_name, controller, scenario), 'Color', 'w');
subplot(2, 2, [1 3]);
if strcmp(scenario, 'circle')
    th = linspace(0, 2 * pi, 400);
    plot(S.RP.radius * sin(th), S.RP.radius * (1 - cos(th)), 'k--'); hold on;
else
    plot(0, 0, 'k+', 'MarkerSize', 12, 'LineWidth', 2); hold on;
end
plot(X(:, 1), X(:, 2), 'Color', [0.16 0.47 0.84]);
axis equal; grid on; xlabel('x [m]'); ylabel('y [m]'); title('Seen from above (dashed: target)');
subplot(2, 2, 2);
plot(t, err, 'Color', [0.16 0.47 0.84]); grid on;
ylabel('Distance from target [m]'); title(sprintf('RMS %.3f m', out.rms_err_m));
subplot(2, 2, 4);
plot(t, U); grid on; xlabel('Time [s]'); ylabel('Thrust per propeller [N]');
end

function v = getlog(so, name)
v = so.get(name);
if ndims(v) == 3
    v = squeeze(v).';
end
end
