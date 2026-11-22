function mdl = build_airship_closed_loop(S)
%BUILD_AIRSHIP_CLOSED_LOOP  Create the full closed-loop Simulink model airship_closed_loop.slx.
%   build_airship_closed_loop(closed_loop_setup('knee8_reversible_capped', 'residual', 'circle'))
%   (run_airship_closed_loop does this for you.)
%
%   Everything runs every 0.1 s (the controller rate); inside the airship block the
%   physics is stepped 10 times at 0.01 s, exactly as in Python:
%
%   [Clock] -> [Target] ---------------------------\
%          \-> [Wind (environment)] ------------\   v
%                                     [State x]-+->[Controller]->[Thrust sharing]->[Airship (plant)]-> x
%                                         ^      ^   memory ^          | thrust, result      |
%                                         |      +---(delay)+----------+ (delay)             |
%                                         +--------------------------------------------------+
%
%   Blocks (each a MATLAB Function block calling the file in brackets):
%     Target            (reference_signal.m)  fixed point with a heading step, or the 50 m circle
%     Wind (environment)(wind_lookup.m)       mean wind + turbulence + gusts, 10 samples per step
%     Controller        (controller_step.m)   classic PID, residual learning or plain learning
%     Thrust sharing    (allocation_step.m)   force and moment -> thrust of each propeller
%     Airship (plant)   (airship_plant_step.m) 6-direction rigid-body model with motors, 10 x 0.01 s
%   Unit Delay blocks hold the state, the controller memory and the last thrust-sharing result.
%   Logged to the workspace: x_log, u_log, ref_log, tau_log.

if isempty(ver('simulink'))
    error('airship:noSimulink', ['Simulink is not installed in this MATLAB. Install it from ' ...
          'Home > Add-Ons, or run run_closed_loop(closed_loop_setup(...)) for the same flight without Simulink.']);
end
here = fileparts(mfilename('fullpath'));
addpath(here);
mdl = 'airship_closed_loop';
if bdIsLoaded(mdl)
    close_system(mdl, 0);
end
N = S.P.n_props; nx = 13 + N;
assignin('base', 'P', S.P); assignin('base', 'CP', S.CP); assignin('base', 'AP', S.AP);
assignin('base', 'RP', S.RP); assignin('base', 'WIND', S.WIND); assignin('base', 'x0', S.x0);
assignin('base', 'u0', zeros(N, 1)); assignin('base', 'cs0', zeros(20, 1)); assignin('base', 'fb0', zeros(12, 1));

new_system(mdl);
set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', num2str(S.dt), ...
          'StopTime', num2str((S.K - 1) * S.dt), 'SaveTime', 'on', 'SaveOutput', 'off', ...
          'ReturnWorkspaceOutputs', 'on');

ts = num2str(S.dt);
add_block('simulink/Sources/Digital Clock', [mdl '/Clock'], 'SampleTime', ts, 'Position', [30 40 80 70]);
fn(mdl, 'Target', [150 20 290 80], ...
   {'function ref = target(t, RP)', '% fixed point with a heading step, or the 50 m circle', ...
    'ref = reference_signal(t, RP);'}, {'RP'}, {'ref', '[11 1]'});
fn(mdl, 'Wind (environment)', [150 330 290 390], ...
   {'function wind_block = wind(t, WIND)', '% world wind for the 10 physics steps of this control step', ...
    'wind_block = wind_lookup(t, WIND);'}, {'WIND'}, {'wind_block', sprintf('[%d 1]', 3 * S.P.n_sub)});
fn(mdl, 'Controller', [380 60 540 200], ...
   {'function [tau, cs_next] = controller(x, ref, fb, cs, CP)', ...
    '% CP.mode: 1 classic (PID), 2 residual learning, 3 plain learning', ...
    '[tau, cs_next] = controller_step(x, ref, fb, cs, CP);'}, {'CP'}, {'tau', '[6 1]'; 'cs_next', '[20 1]'});
fn(mdl, 'Thrust sharing', [620 70 760 170], ...
   {'function [u, fb] = thrust_sharing(tau, u_prev, AP)', ...
    '% force and moment -> thrust of each propeller, within the thrust limits', ...
    '[u, achieved, gap] = allocation_step(tau, u_prev, AP);', 'fb = [achieved; gap];'}, ...
   {'AP'}, {'u', sprintf('[%d 1]', N); 'fb', '[12 1]'});
fn(mdl, 'Airship (plant)', [840 80 1000 200], ...
   {'function x_next = airship(x, u, wind_block, P)', ...
    '% 6-direction rigid-body airship with motors: 10 Runge-Kutta steps of 0.01 s', ...
    'x_next = airship_plant_step(x, u, wind_block, P);'}, {'P'}, {'x_next', sprintf('[%d 1]', nx)});

add_block('simulink/Discrete/Unit Delay', [mdl '/State x'], 'InitialCondition', 'x0', ...
          'SampleTime', ts, 'Position', [1060 230 1110 270], 'Orientation', 'left');
add_block('simulink/Discrete/Unit Delay', [mdl '/Controller memory'], 'InitialCondition', 'cs0', ...
          'SampleTime', ts, 'Position', [430 250 480 290], 'Orientation', 'left');
add_block('simulink/Discrete/Unit Delay', [mdl '/Last thrust'], 'InitialCondition', 'u0', ...
          'SampleTime', ts, 'Position', [660 230 710 270], 'Orientation', 'left');
add_block('simulink/Discrete/Unit Delay', [mdl '/Last sharing result'], 'InitialCondition', 'fb0', ...
          'SampleTime', ts, 'Position', [560 300 610 340], 'Orientation', 'left');
logs = {'x_log', [1120 20 1200 50]; 'u_log', [840 10 920 40]; 'ref_log', [380 0 460 30]; 'tau_log', [620 0 700 30]};
for i = 1:size(logs, 1)
    add_block('simulink/Sinks/To Workspace', [mdl '/Log ' logs{i, 1}], 'VariableName', logs{i, 1}, ...
              'SaveFormat', 'Array', 'SampleTime', ts, 'Position', logs{i, 2});
end

L = @(a, b) add_line(mdl, a, b, 'autorouting', 'on');
L('Clock/1', 'Target/1');
L('Clock/1', 'Wind (environment)/1');
L('State x/1', 'Controller/1');
L('Target/1', 'Controller/2');
L('Last sharing result/1', 'Controller/3');
L('Controller memory/1', 'Controller/4');
L('Controller/1', 'Thrust sharing/1');
L('Controller/2', 'Controller memory/1');
L('Last thrust/1', 'Thrust sharing/2');
L('Thrust sharing/1', 'Airship (plant)/2');
L('Thrust sharing/1', 'Last thrust/1');
L('Thrust sharing/2', 'Last sharing result/1');
L('State x/1', 'Airship (plant)/1');
L('Wind (environment)/1', 'Airship (plant)/3');
L('Airship (plant)/1', 'State x/1');
L('Airship (plant)/1', 'Log x_log/1');
L('Thrust sharing/1', 'Log u_log/1');
L('Target/1', 'Log ref_log/1');
L('Controller/1', 'Log tau_log/1');

% scope: position (x, y, z) and the target position
try
    add_block('simulink/Signal Routing/Selector', [mdl '/Position'], 'NumberOfDimensions', '1', ...
              'IndexOptions', 'Index vector (dialog)', 'Indices', '[1 2 3]', 'InputPortWidth', num2str(nx), ...
              'Position', [1120 120 1150 150]);
    add_block('simulink/Signal Routing/Selector', [mdl '/Target position'], 'NumberOfDimensions', '1', ...
              'IndexOptions', 'Index vector (dialog)', 'Indices', '[1 2 3]', 'InputPortWidth', '11', ...
              'Position', [1120 170 1150 200]);
    add_block('simulink/Signal Routing/Mux', [mdl '/Mux'], 'Inputs', '2', 'Position', [1180 120 1185 200]);
    add_block('simulink/Sinks/Scope', [mdl '/Position vs target'], 'Position', [1220 140 1260 180]);
    L('Airship (plant)/1', 'Position/1');
    L('Target/1', 'Target position/1');
    L('Position/1', 'Mux/1');
    L('Target position/1', 'Mux/2');
    L('Mux/1', 'Position vs target/1');
catch err
    warning('Scope not added (%s); the rest of the model is complete.', err.message);
end

save_system(mdl, fullfile(here, [mdl '.slx']));
fprintf('Built %s.slx: %s, %s, %s (%d propellers)\n', mdl, S.layout, S.controller, S.scenario, N);
end

function fn(mdl, name, pos, code, params, outs)
% add a MATLAB Function block with the given code, parameter inputs and output sizes
blk = [mdl '/' name];
add_block('simulink/User-Defined Functions/MATLAB Function', blk, 'Position', pos);
ch = find(sfroot, '-isa', 'Stateflow.EMChart', 'Path', blk);
ch.Script = sprintf('%s\n', code{:});
for i = 1:numel(params)
    d = find(ch, '-isa', 'Stateflow.Data', 'Name', params{i});
    d.Scope = 'Parameter';
end
for i = 1:size(outs, 1)
    d = find(ch, '-isa', 'Stateflow.Data', 'Name', outs{i, 1});
    d.Props.Array.Size = outs{i, 2};
end
end
