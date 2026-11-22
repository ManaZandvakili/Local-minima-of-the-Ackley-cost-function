function S = closed_loop_setup(layout_name, controller, scenario, wind_seed)
%CLOSED_LOOP_SETUP  Everything one closed-loop flight needs, in one struct.
%   S = closed_loop_setup(layout_name, controller, scenario)
%   S = closed_loop_setup(layout_name, controller, scenario, wind_seed)
%
%   layout_name  'knee8_reversible_capped' (optimized, 8 propellers), 'random8_reversible_typical'
%                (typical random layout, 8 propellers), or any layout in data/layouts/
%   controller   'pid' (classic), 'residual' (residual learning) or 'plain' (plain learning)
%   scenario     'station' (hold position 60 s, heading step 45 deg at 30 s, start 1.2 m off)
%                'circle'  (the paper's 50 m circle at 0.5 m/s, 650 s)
%   Wind: 3 m/s toward 0 deg, turbulence 0.5 m/s, 2 gusts per minute. Without wind_seed the
%   exact wind of the Python run is used (data/validation/wind_<scenario>.mat), so the numbers
%   match Python; with wind_seed, wind_model.m makes a new wind with MATLAB's random numbers.
%
%   Fields: P (airship + layout), CP (controller), AP (thrust sharing), RP (target),
%   WIND (wind table), x0 (initial state), K (number of 0.1 s steps), dt.

if nargin < 1 || isempty(layout_name), layout_name = 'knee8_reversible_capped'; end
if nargin < 2 || isempty(controller), controller = 'residual'; end
if nargin < 3 || isempty(scenario), scenario = 'station'; end
here = fileparts(mfilename('fullpath'));
P = airship_load(layout_name);
N = P.n_props;
dt = P.n_sub * P.dt;                               % 0.1 s control step
if ~isempty(strfind(layout_name, 'oneway')), kind = 'oneway'; else, kind = 'reversible'; end

% thrust sharing (Allocator in airship/allocation.py)
L = load(fullfile(here, 'data', 'layouts', [layout_name '.mat']));
AP.A = P.A;
AP.W2 = ([2; 2; 1; 0.5; 0.5; 2]) .^ 2;
AP.H = P.A.' * (AP.W2 .* P.A) + 1e-8 * eye(N);
AP.lo = P.u_lo;
AP.hi = P.u_hi;

% controller
modes = {'pid', 'residual', 'plain'};
CP.mode = find(strcmp(modes, controller));
if isempty(CP.mode), error('controller must be pid, residual or plain'); end
CP.dt = dt;
G = load(fullfile(here, 'data', 'controllers', ['pid_' layout_name '.mat']));
CP.pid = rmfield(G.pid, {'gains', 'source_file'});
CP.caps = L.layout.caps(:);
R = load(fullfile(here, 'data', 'controllers', ['residual_' kind '.mat']));
Q = load(fullfile(here, 'data', 'controllers', ['plain_' kind '.mat']));
CP.residual = rmfield(R.net, 'source_run');
CP.plain = rmfield(Q.net, 'source_run');

% target and duration
RP.dt = dt;
RP.yaw_step_time = 30.0; RP.yaw_step_to = deg2rad(45);
RP.radius = 50.0; RP.speed = 0.5;
x0 = zeros(13 + N, 1); x0(4) = 1;
if strcmp(scenario, 'station')
    RP.kind = 1; K = 600; x0(1:3) = [1.0; -0.6; 0.3]; py_seed = 500;
else
    RP.kind = 2; K = 6500; py_seed = 7;
end

% wind
if nargin < 4 || isempty(wind_seed)
    W = load(fullfile(here, 'data', 'validation', ['wind_' scenario '.mat']));
    w = W.wind;
    wind_note = sprintf('Python wind (seed %d)', py_seed);
else
    cfg = struct('mean_speed_m_s', 3.0, 'mean_dir_deg', 0.0, 'turbulence_sigma_m_s', [0.5 0.5 0.25], ...
                 'gust_rate_per_min', 2.0);
    w = wind_model(cfg, wind_seed, K * P.n_sub, P.dt);
    wind_note = sprintf('MATLAB wind (seed %d)', wind_seed);
end
WIND.dt = dt;
WIND.table = reshape(w(1:K * P.n_sub, :).', 3 * P.n_sub, K).';   % row k = wind of control step k

S = struct('P', P, 'CP', CP, 'AP', AP, 'RP', RP, 'WIND', WIND, 'x0', x0, 'K', K, 'dt', dt, ...
           'layout', layout_name, 'controller', controller, 'scenario', scenario, 'wind', wind_note);
end
