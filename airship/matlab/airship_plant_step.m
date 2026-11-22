function x = airship_plant_step(x, u_cmd, wind_block, P)
%AIRSHIP_PLANT_STEP  The airship over one 0.1 s control step: P.n_sub (10) Runge-Kutta steps of 0.01 s.
%   x = airship_plant_step(x, u_cmd, wind_block, P)
%   wind_block = the world wind for the 10 physics steps, [w1; w2; ...; w10] (30 x 1).
%   Same as Simulator.step in airship/airship/dynamics.py.
%#codegen

for j = 1:P.n_sub
    x = airship_rk4_step(x, u_cmd, wind_block(3 * j - 2:3 * j), P);
end
end
