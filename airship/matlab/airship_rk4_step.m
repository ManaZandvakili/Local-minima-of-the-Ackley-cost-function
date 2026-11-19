function x = airship_rk4_step(x, u_cmd, wind_w, P)
%AIRSHIP_RK4_STEP  Advance the airship by one physics step of P.dt seconds (0.01 s).
%   x = airship_rk4_step(x, u_cmd, wind_w, P)
%
%   Same as one pass of integrate() in airship/airship/dynamics.py:
%   the command is clipped to the thrust limits, the state is advanced with
%   4th-order Runge-Kutta, the quaternion is re-normalised and the motor
%   thrusts are clipped to their limits. The wind is held constant over the step.
%   This is the function the Simulink model calls every 0.01 s.
%#codegen

x = x(:);
lo = P.u_lo(:);
hi = P.u_hi(:);
u_cmd = min(max(u_cmd(:), lo), hi);
wind_w = wind_w(:);
dt = P.dt;
k1 = airship_derivative(x, u_cmd, wind_w, P);
k2 = airship_derivative(x + 0.5 * dt * k1, u_cmd, wind_w, P);
k3 = airship_derivative(x + 0.5 * dt * k2, u_cmd, wind_w, P);
k4 = airship_derivative(x + dt * k3, u_cmd, wind_w, P);
x = x + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
x(4:7) = x(4:7) / sqrt(sum(x(4:7) .^ 2));
x(14:end) = min(max(x(14:end), lo), hi);
end
