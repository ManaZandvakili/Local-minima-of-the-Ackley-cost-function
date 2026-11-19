function dx = airship_derivative(x, u_cmd, wind_w, P)
%AIRSHIP_DERIVATIVE  Time derivative of the sphere airship state.
%   dx = airship_derivative(x, u_cmd, wind_w, P)
%
%   Line-by-line port of derivative() in airship/airship/dynamics.py.
%
%   State x (column, 13 + N):
%     x(1:3)    position in the world frame [m] (z up)
%     x(4:7)    attitude quaternion [w x y z], body -> world
%     x(8:10)   linear velocity in the body frame [m/s]
%     x(11:13)  angular velocity in the body frame [rad/s]
%     x(14:end) actual thrust of each propeller [N]
%   u_cmd  commanded thrust per propeller [N] (column, N)
%   wind_w wind velocity in the world frame [m/s] (column, 3)
%   P      parameter struct from airship_load()
%
%   Equations (body frame, nu = [v; w], nu_r = nu - [R'*wind; 0]):
%     (Mrb + Ma) nu_dot = A*u + buoyancy/weight + drag + damping
%                         - C_rb(nu) nu - C_a(nu_r) nu_r
%     u_dot = (u_cmd - u) / tau_m        (motor lag)
%#codegen

q = x(4:7);
v = x(8:10);
w = x(11:13);
u = x(14:end);
R = quat_to_rot(q);
nu = x(8:13);
v_air = R.' * wind_w(:);
nu_r = nu;
nu_r(1:3) = v - v_air;

% propeller thrust and a constant disturbance torque
tau = P.A * u;
tau(4:6) = tau(4:6) + P.tau_dist(:);

% buoyancy at the hull centre, weight at the centre of gravity
up = R.' * [0; 0; 1];
tau(1:3) = tau(1:3) + (P.buoy - P.weight) * up;
tau(4:6) = tau(4:6) + cross3(P.r_g(:), -P.weight * up);

% aerodynamic drag on the air-relative velocity, rotational damping
sp = sqrt(nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2);
tau(1:3) = tau(1:3) - P.k_drag * sp * nu_r(1:3);
tau(4:6) = tau(4:6) - P.d_lin(:) .* w - P.d_quad(:) .* abs(w) .* w;

% Coriolis and centripetal terms (rigid body + added mass)
tau = tau - coriolis(P.Mrb, nu) - coriolis(P.Ma, nu_r);

dx = zeros(size(x));
dx(1:3) = R * v;
dx(4) = -0.5 * (q(2) * w(1) + q(3) * w(2) + q(4) * w(3));
dx(5) =  0.5 * (q(1) * w(1) + q(3) * w(3) - q(4) * w(2));
dx(6) =  0.5 * (q(1) * w(2) + q(4) * w(1) - q(2) * w(3));
dx(7) =  0.5 * (q(1) * w(3) + q(2) * w(2) - q(3) * w(1));
dx(8:13) = P.Minv * tau;
dx(14:end) = (u_cmd(:) - u) / P.tau_m;
end

function c = coriolis(M, nu)
% C(nu) nu for a symmetric 6x6 mass matrix: [w x n1; v x n1 + w x n2], n = M nu
n = M * nu;
c = zeros(6, 1);
c(1:3) = cross3(nu(4:6), n(1:3));
c(4:6) = cross3(nu(1:3), n(1:3)) + cross3(nu(4:6), n(4:6));
end

function c = cross3(a, b)
c = [a(2) * b(3) - a(3) * b(2); a(3) * b(1) - a(1) * b(3); a(1) * b(2) - a(2) * b(1)];
end
