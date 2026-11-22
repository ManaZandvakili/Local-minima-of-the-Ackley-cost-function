%COMPARE_CONTROLLERS  The main comparison of the journal paper, in plain MATLAB (no Simulink needed).
%   Optimized 8-propeller layout vs the typical random 8-propeller layout (reversible propellers),
%   each flown by the classic controller, residual learning and plain learning, in gusty
%   station keeping and on the paper's 50 m circle (3 m/s wind, exact Python wind).
%   Prints a table and draws a bar chart. Takes a few minutes.

here = fileparts(mfilename('fullpath'));
addpath(here);
layouts = {'knee8_reversible_capped', 'random8_reversible_typical'};
lname = {'Optimized', 'Typical random'};
ctrls = {'pid', 'residual', 'plain'};
cname = {'Classic (PID)', 'Residual learning', 'Plain learning'};
scen = {'station', 'circle'};
sname = {'Holding position', 'Circle'};
E = nan(2, 2, 3); done = false(2, 2, 3);
for s = 1:2
    for l = 1:2
        for c = 1:3
            out = run_closed_loop(closed_loop_setup(layouts{l}, ctrls{c}, scen{s}));
            E(s, l, c) = out.rms_err_m;
            done(s, l, c) = out.completed || strcmp(scen{s}, 'station') && out.survived;
            note = '';
            if ~out.survived
                note = sprintf(' (lost control after %.0f s)', out.T(end));
            elseif ~done(s, l, c)
                note = ' (lap not finished: fell more than 3 m behind)';
            end
            fprintf('%-17s %-15s %-18s RMS error %.3f m%s\n', sname{s}, lname{l}, cname{c}, out.rms_err_m, note);
        end
    end
end
figure('Name', 'Optimized vs random, three controllers', 'Color', 'w');
for s = 1:2
    subplot(1, 2, s);
    b = bar(squeeze(E(s, :, :)).');
    b(1).FaceColor = [0.16 0.47 0.84]; b(2).FaceColor = [0.65 0.65 0.65];
    set(gca, 'XTickLabel', cname); ylabel('RMS distance from target [m]'); grid on;
    legend(lname, 'Location', 'northwest'); title([sname{s} ', 3 m/s gusty wind']);
end
