function obs = make_observation(x, ref, integral, a_prev, gap, caps, NN, extra)
%MAKE_OBSERVATION  What the learning controller sees (41 numbers, 47 with the PID's wrench).
%   Same as make_obs in airship/airship/envs/station_keeping.py, called the way
%   RLController.act (airship/airship/controllers/rl_policy.py) calls it: speed and
%   turn rate relative to the (moving) target. Rounded to single precision, as in Python.
%#codegen

R = quat_to_rot(x(4:7));
p_err = x(1:3) - ref(1:3);
v_rel = x(8:10) - R.' * ref(4:6);
w_rel = x(11:13) - R.' * [0; 0; ref(11)];
yaw = atan2(R(2, 1), R(1, 1));
ey = mod(ref(10) - yaw + pi, 2 * pi) - pi;
scale = NN.action_scale;
o = [R.' * (-p_err) / NN.pos_scale_m; v_rel / NN.vel_scale_m_s; R.' * [0; 0; -1]; sin(ey); cos(ey);
     w_rel / NN.rate_scale_rad_s; R.' * integral / NN.integral_clip_m_s; a_prev; gap ./ scale;
     caps ./ kron(scale, [1; 1])];
if nargin > 7
    o = [o; extra];
end
obs = double(single(o));
end
