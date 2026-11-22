function out = run_closed_loop(S)
%RUN_CLOSED_LOOP  Fly one closed-loop flight in plain MATLAB (no Simulink needed).
%   out = run_closed_loop(closed_loop_setup('knee8_reversible_capped', 'residual', 'circle'))
%
%   Every 0.1 s: target -> controller -> thrust sharing -> airship (10 steps of 0.01 s).
%   Exactly the loop of the Simulink model airship_closed_loop.slx, and of
%   run_episode / run_path_episode in the Python code.
%   out: X (state after each step), U (thrust commands), TAU (wanted force and moment),
%   REF (target), T (time), and the scores: rms_err_m, max_err_m, effort_N2s, survived,
%   completed (circle: error under 3 m the whole lap).

K = S.K; dt = S.dt; N = S.P.n_props;
x = S.x0; cs = zeros(20, 1); fb = zeros(12, 1); u_prev = zeros(N, 1);
X = zeros(K, 13 + N); U = zeros(K, N); TAU = zeros(K, 6); REF = zeros(K, 11); err = zeros(K, 1);
n = K; survived = true;
for k = 1:K
    t = (k - 1) * dt;
    ref = reference_signal(t, S.RP);
    [tau, cs] = controller_step(x, ref, fb, cs, S.CP);
    [u, achieved, gap] = allocation_step(tau, u_prev, S.AP);
    u_prev = u; fb = [achieved; gap];
    x = airship_plant_step(x, u, wind_lookup(t, S.WIND), S.P);
    X(k, :) = x.'; U(k, :) = u.'; TAU(k, :) = tau.'; REF(k, :) = ref.';
    nxt = reference_signal(t + dt, S.RP);
    err(k) = norm(x(1:3) - nxt(1:3));
    Rm = quat_to_rot(x(4:7));
    tilt = acosd(max(min(Rm(3, 3), 1), -1));
    if any(~isfinite(x)) || err(k) > 10 || tilt > 60
        survived = false; n = k; break
    end
end
out.X = X(1:n, :); out.U = U(1:n, :); out.TAU = TAU(1:n, :); out.REF = REF(1:n, :);
out.T = (1:n).' * dt; out.err = err(1:n);
out.rms_err_m = sqrt(mean(err(1:n) .^ 2));
out.max_err_m = max(err(1:n));
out.effort_N2s = sum(sum(out.X(:, 14:end) .^ 2)) * dt;
out.survived = survived;
out.completed = survived && out.max_err_m < 3;
end
