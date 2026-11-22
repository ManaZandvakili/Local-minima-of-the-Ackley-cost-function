function ref = reference_signal(t, RP)
%REFERENCE_SIGNAL  Target for the controller at time t (11 x 1).
%   ref = [position (3); velocity (3); acceleration (3); heading; heading rate]
%   RP.kind = 1: station keeping at the origin; the heading command steps from 0 to
%                RP.yaw_step_to [rad] at RP.yaw_step_time [s].
%   RP.kind = 2: the conference paper's circle (Final.py, _path_reference):
%                radius RP.radius [m], speed RP.speed [m/s], nose along the path.
%   Same as airship/closed_loop.py (Scenario.yaw_at) and airship/path_following.py.
%#codegen

k = round(t / RP.dt);              % control step number (avoids rounding in t)
t = k * RP.dt;
ref = zeros(11, 1);
if RP.kind == 1
    if t >= RP.yaw_step_time
        ref(10) = RP.yaw_step_to;
    end
else
    R = RP.radius; V = RP.speed; w = V / R;
    s = sin(w * t); c = cos(w * t);
    ref(1:3) = [R * s; R * (1 - c); 0];
    ref(4:6) = [V * c; V * s; 0];
    ref(7:9) = [-V * w * s; V * w * c; 0];
    ref(10) = w * t;
    ref(11) = w;
end
end
