function [tau, cs] = controller_step(x, ref, fb, cs, CP)
%CONTROLLER_STEP  One 0.1 s step of the chosen controller -> wanted force and moment (6 x 1).
%   [tau, cs] = controller_step(x, ref, fb, cs, CP)
%   x    state (13+N);  ref  target from reference_signal (11 x 1)
%   fb   last thrust-sharing result [achieved (6); shortfall (6)] (zeros at the start)
%   cs   controller memory (20 x 1, zeros at the start):
%        PID integrals (1:4), learning integral (5:7), last action (8:13),
%        last PID wrench (14:19), first-step flag (20, 0 = first step not done yet)
%   CP.mode  1 = classic (PID), 2 = residual learning (PID + correction), 3 = plain learning
%   Same as PIDController and RLController in airship/airship/controllers/.
%#codegen

achieved = fb(1:6);
gap = fb(7:12);
% anti-windup flag from the last thrust sharing (PIDController.notify)
saturated = sqrt(sum(gap(1:3) .^ 2)) > 0.05 * sqrt(sum((achieved(1:3) + gap(1:3)) .^ 2)) + 0.05;
if CP.mode == 1
    [tau, s] = pid_step(x, ref, cs(1:4), CP.pid, saturated);
    cs(1:4) = s;
    return
end
% learning controllers (RLController.act)
if CP.mode == 2
    NN = CP.residual;
else
    NN = CP.plain;
end
integral = cs(5:7);
if cs(20) > 0
    integral = min(max(integral + (ref(1:3) - x(1:3)) * CP.dt, -NN.integral_clip_m_s), NN.integral_clip_m_s);
end
a_prev = cs(8:13);
if CP.mode == 2
    obs = make_observation(x, ref, integral, a_prev, -gap, CP.caps, NN, cs(14:19) ./ NN.action_scale);
    a = policy_action(obs, CP.residual);
    [pid_w, s] = pid_step(x, ref, cs(1:4), CP.pid, saturated);
    tau = pid_w + a .* NN.residual_scale;
    cs(1:4) = s;
    cs(14:19) = pid_w;
else
    obs = make_observation(x, ref, integral, a_prev, -gap, CP.caps, NN);
    a = policy_action(obs, CP.plain);
    tau = a .* NN.action_scale;
end
cs(5:7) = integral;
cs(8:13) = a;
cs(20) = 1;
end
