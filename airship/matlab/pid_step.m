function [tau, s] = pid_step(x, ref, s, G, saturated)
%PID_STEP  Classic controller: cascaded PID -> wanted force and moment (body frame).
%   [tau, s] = pid_step(x, ref, s, G, saturated)
%   x    state (13+N), ref from reference_signal, s = [position integral (3); heading integral]
%   G    gains from data/controllers/pid_<layout>.mat (struct pid)
%   saturated  true while the thrust sharing could not deliver the force (anti-windup)
%   Same as PIDController.act in airship/airship/controllers/pid.py.
%#codegen

p_ref = ref(1:3); v_ref = ref(4:6); a_ref = ref(7:9); yaw_ref = ref(10); yaw_rate_ref = ref(11);
R = quat_to_rot(x(4:7));
e = p_ref - x(1:3);
v_world = R * x(8:10);
if ~saturated
    s(1:3) = min(max(s(1:3) + e * G.dt, -G.i_max_pos), G.i_max_pos);
end
f_world = G.mass_eff * (G.kp .* e + G.ki .* s(1:3) + G.kd .* (v_ref - v_world) + a_ref);
w = x(11:13);
ea = attitude_error(x(4:7), yaw_ref);
if ~saturated
    s(4) = min(max(s(4) + ea(3) * G.dt, -G.i_max_yaw), G.i_max_yaw);
end
torque = G.inertia .* [G.kp_rp * ea(1) - G.kd_rp * w(1);
                       G.kp_rp * ea(2) - G.kd_rp * w(2);
                       G.kp_yaw * ea(3) + G.ki_yaw * s(4) - G.kd_yaw * (w(3) - yaw_rate_ref)];
tau = [R.' * f_world; torque];
end

function ea = attitude_error(q, yaw_ref)
% small-angle rotation vector (body frame) from the attitude q to level at yaw_ref
q_ref = [cos(yaw_ref / 2); 0; 0; sin(yaw_ref / 2)];
a = [q(1); -q(2); -q(3); -q(4)];
b = q_ref;
qe = [a(1) * b(1) - a(2) * b(2) - a(3) * b(3) - a(4) * b(4);
      a(1) * b(2) + a(2) * b(1) + a(3) * b(4) - a(4) * b(3);
      a(1) * b(3) - a(2) * b(4) + a(3) * b(1) + a(4) * b(2);
      a(1) * b(4) + a(2) * b(3) - a(3) * b(2) + a(4) * b(1)];
if qe(1) < 0
    qe = -qe;
end
sn = sqrt(sum(qe(2:4) .^ 2));
if sn < 1e-12
    ea = zeros(3, 1);
else
    ea = 2.0 * atan2(sn, qe(1)) * qe(2:4) / sn;
end
end
