%VALIDATE_CLOSED_LOOP  Check the MATLAB controllers and closed loop against Python.
%   Run from the matlab/ folder (plain MATLAB, no Simulink needed).
%   Flies classic, residual and plain on the optimized 8-propeller reversible layout,
%   in gusty station keeping and on the 50 m circle, with the exact Python wind, and
%   compares every state with the Python flight (data/validation/python_*.mat).

here = fileparts(mfilename('fullpath'));
addpath(here);
names = {'pid', 'residual', 'plain'};
label = struct('pid', 'Classic (PID)', 'residual', 'Residual learning', 'plain', 'Plain learning');
all_ok = true;
for sc = {'station', 'circle'}
    for c = names
        S = closed_loop_setup('knee8_reversible_capped', c{1}, sc{1});
        tic; out = run_closed_loop(S); el = toc;
        Q = load(fullfile(here, 'data', 'validation', ['python_' c{1} '_' sc{1} '.mat']));
        n = min(size(out.X, 1), size(Q.py.X, 1));
        dpos = max(max(abs(out.X(1:n, 1:3) - Q.py.X(1:n, 1:3))));
        if strcmp(sc{1}, 'circle')
            w = S.RP.speed / S.RP.radius; tt = (1:n).' * S.dt;
            nxt = [S.RP.radius * sin(w * tt), S.RP.radius * (1 - cos(w * tt)), zeros(n, 1)];
        else
            nxt = zeros(n, 3);
        end
        py_rms = sqrt(mean(sum((Q.py.X(1:n, 1:3) - nxt) .^ 2, 2)));
        ok = dpos < 1e-3 && n == size(Q.py.X, 1);
        all_ok = all_ok && ok;
        verdict = 'FAIL';
        if ok, verdict = 'PASS'; end
        fprintf('%-8s %-18s RMS error MATLAB %.4f m, Python %.4f m, largest position difference %.1e m (%.0f s to run)  %s\n', ...
                sc{1}, label.(c{1}), out.rms_err_m, py_rms, dpos, el, verdict);
    end
end
if all_ok
    fprintf('All checks PASS: the MATLAB controllers fly exactly like the Python ones.\n');
end

