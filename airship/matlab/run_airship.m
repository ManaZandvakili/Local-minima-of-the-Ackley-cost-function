function [X, T] = run_airship(P, x0, cmd, wind)
%RUN_AIRSHIP  Fly the simulator in plain MATLAB (no Simulink needed).
%   [X, T] = run_airship(P, x0, cmd, wind)
%
%   cmd   K x N  thrust commands, one row per control step (held 0.1 s)
%   wind  (K*P.n_sub) x 3  world wind, one row per physics step (held 0.01 s)
%   X     K x (13+N) state after each control step, T its time [s]
%
%   Same loop as Simulator.step() in airship/airship/dynamics.py.

K = size(cmd, 1);
x = x0(:);
X = zeros(K, numel(x));
for k = 1:K
    u = cmd(k, :).';
    for j = 1:P.n_sub
        x = airship_rk4_step(x, u, wind((k - 1) * P.n_sub + j, :).', P);
    end
    X(k, :) = x.';
end
T = (1:K).' * P.n_sub * P.dt;
end
