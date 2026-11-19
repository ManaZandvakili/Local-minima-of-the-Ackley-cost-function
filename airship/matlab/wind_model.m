function W = wind_model(cfg, seed, n, dt)
%WIND_MODEL  Wind samples: mean wind + turbulence + 1-cosine gusts.
%   W = wind_model(cfg, seed, n, dt)   ->  n x 3 world wind [m/s], one row per step
%
%   cfg fields (defaults in brackets):
%     mean_speed_m_s [0], mean_dir_deg [0] (direction the wind blows toward, from +x, CCW),
%     turbulence_sigma_m_s [0.5 0.5 0.25], turbulence_tau_s [2],
%     gust_rate_per_min [0], gust_amplitude_m_s [1 3], gust_duration_s [2 6]
%
%   Same model as airship/airship/wind.py (turbulence = Ornstein-Uhlenbeck
%   process with exact discrete update; gusts drawn up front). MATLAB's random
%   numbers differ from Python's, so the wind is statistically the same but not
%   sample for sample. To replay a Python flight exactly, use the wind stored in
%   data/reference_flight.mat instead.

d = struct('mean_speed_m_s', 0, 'mean_dir_deg', 0, 'turbulence_sigma_m_s', [0.5 0.5 0.25], ...
           'turbulence_tau_s', 2, 'gust_rate_per_min', 0, 'gust_amplitude_m_s', [1 3], ...
           'gust_duration_s', [2 6]);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}), cfg.(f{i}) = d.(f{i}); end
end
rng(seed);
a = deg2rad(cfg.mean_dir_deg);
mean_w = cfg.mean_speed_m_s * [cos(a) sin(a) 0];
sigma = cfg.turbulence_sigma_m_s(:).';
turb = randn(1, 3) .* sigma;                      % start in the stationary distribution

% gust schedule: start time, duration, peak vector
gusts = zeros(0, 5);
if cfg.gust_rate_per_min > 0
    t = 0;
    while true
        t = t - log(rand) * 60 / cfg.gust_rate_per_min;   % exponential waiting time
        if t >= n * dt, break; end
        dur = cfg.gust_duration_s(1) + diff(cfg.gust_duration_s) * rand;
        amp = cfg.gust_amplitude_m_s(1) + diff(cfg.gust_amplitude_m_s) * rand;
        ang = 2 * pi * rand;
        gusts(end + 1, :) = [t dur amp * [cos(ang) sin(ang) 0]]; %#ok<AGROW>
    end
end

aa = exp(-dt / cfg.turbulence_tau_s);
bb = sigma * sqrt(1 - aa^2);
W = zeros(n, 3);
for k = 1:n
    t = (k - 1) * dt;
    g = [0 0 0];
    for i = 1:size(gusts, 1)
        if t >= gusts(i, 1) && t < gusts(i, 1) + gusts(i, 2)
            g = g + 0.5 * gusts(i, 3:5) * (1 - cos(2 * pi * (t - gusts(i, 1)) / gusts(i, 2)));
        end
    end
    W(k, :) = mean_w + turb + g;
    turb = aa * turb + bb .* randn(1, 3);          % exact OU update
end
end
