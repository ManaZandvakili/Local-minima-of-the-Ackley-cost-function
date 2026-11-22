function mdl = build_airship_simulink(layout_name)
%BUILD_AIRSHIP_SIMULINK  Create the Simulink model airship_sim.slx.
%   build_airship_simulink                      % 14-propeller reversible knee layout
%   build_airship_simulink('knee9_oneway_capped')
%
%   The model is a fixed-step discrete loop at 0.01 s (100 Hz), the same as the
%   Python simulator:
%
%     [Thrust command] ---\
%     [Wind] -------------> [Airship step: one 0.01 s RK4 step] --> [State x (unit delay)] --+--> [Log x]
%                        ^                                                                    |
%                        +--------------------------------------------------------------------+
%                                                                                             +--> [Position scope]
%
%   "Airship step" is a MATLAB Function block that calls airship_rk4_step.m,
%   so the Simulink model runs exactly the same equations as run_airship.m
%   and the Python code.
%
%   Inputs (base workspace, [time, values] matrices; values are held between samples):
%     cmd_data   [t, u_1 ... u_N]  thrust commands [N] (one row per 0.1 s control step)
%     wind_data  [t, wx wy wz]     world wind [m/s]
%     x0         initial state (13+N): position, quaternion [w x y z], body velocity,
%                body angular rate, propeller thrusts
%     P          parameters from airship_load()
%   Output: x_log, the state at every 0.01 s.
%
%   The model starts with zero thrust and no wind for 60 s; run_airship_simulink.m
%   replays the Python reference flight through it. To fly your own controller,
%   replace the "Thrust command" block with your controller and thrust sharing,
%   fed from the "State x" signal.

if nargin < 1
    layout_name = 'knee14_reversible';
end
if isempty(ver('simulink'))
    error('airship:noSimulink', ['Simulink is not installed in this MATLAB. Install it from ' ...
          'Home > Add-Ons > Get Add-Ons (search "Simulink"), then run this again. ' ...
          'Without Simulink, run plot_reference_flight to see the same flight in plain MATLAB.']);
end
here = fileparts(mfilename('fullpath'));
addpath(here);
mdl = 'airship_sim';
if bdIsLoaded(mdl)
    close_system(mdl, 0);
end
P = airship_load(layout_name);
N = P.n_props;
nx = 13 + N;
x0 = zeros(nx, 1);
x0(4) = 1;                                       % level, heading 0
assignin('base', 'P', P);
assignin('base', 'x0', x0);
assignin('base', 'cmd_data', [0, zeros(1, N)]);
assignin('base', 'wind_data', [0, 0, 0, 0]);

new_system(mdl);
set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', '0.01', 'StopTime', '60', ...
          'SaveTime', 'on', 'SaveOutput', 'off', 'ReturnWorkspaceOutputs', 'on');

add_block('simulink/Sources/From Workspace', [mdl '/Thrust command'], ...
          'VariableName', 'cmd_data', 'SampleTime', '0.01', 'Interpolate', 'off', ...
          'OutputAfterFinalValue', 'Holding final value', 'Position', [40 95 150 125]);
add_block('simulink/Sources/From Workspace', [mdl '/Wind'], ...
          'VariableName', 'wind_data', 'SampleTime', '0.01', 'Interpolate', 'off', ...
          'OutputAfterFinalValue', 'Holding final value', 'Position', [40 155 150 185]);
add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/Airship step'], ...
          'Position', [240 60 400 190]);
add_block('simulink/Discrete/Unit Delay', [mdl '/State x'], ...
          'InitialCondition', 'x0', 'SampleTime', '0.01', 'Position', [470 105 520 145]);
add_block('simulink/Sinks/To Workspace', [mdl '/Log x'], ...
          'VariableName', 'x_log', 'SaveFormat', 'Array', 'SampleTime', '0.01', ...
          'Position', [600 40 690 70]);

% code of the MATLAB Function block: one physics step
ch = find(sfroot, '-isa', 'Stateflow.EMChart', 'Path', [mdl '/Airship step']);
ch.Script = sprintf(['function x_next = airship_step(x, u_cmd, wind, P)\n' ...
                     '%% One 0.01 s step of the sphere airship (see airship_rk4_step.m)\n' ...
                     'x_next = airship_rk4_step(x, u_cmd, wind, P);\n']);
pd = find(ch, '-isa', 'Stateflow.Data', 'Name', 'P');
pd.Scope = 'Parameter';
xo = find(ch, '-isa', 'Stateflow.Data', 'Name', 'x_next');
xo.Props.Array.Size = sprintf('[%d 1]', nx);

add_line(mdl, 'State x/1', 'Airship step/1', 'autorouting', 'on');
add_line(mdl, 'Thrust command/1', 'Airship step/2', 'autorouting', 'on');
add_line(mdl, 'Wind/1', 'Airship step/3', 'autorouting', 'on');
add_line(mdl, 'Airship step/1', 'State x/1', 'autorouting', 'on');
add_line(mdl, 'State x/1', 'Log x/1', 'autorouting', 'on');

% position (first three states) on a scope
try
    add_block('simulink/Signal Routing/Selector', [mdl '/Position'], ...
              'NumberOfDimensions', '1', 'IndexOptions', 'Index vector (dialog)', ...
              'Indices', '[1 2 3]', 'InputPortWidth', num2str(nx), 'Position', [600 160 640 190]);
    add_block('simulink/Sinks/Scope', [mdl '/Position scope'], 'Position', [700 155 740 195]);
    add_line(mdl, 'State x/1', 'Position/1', 'autorouting', 'on');
    add_line(mdl, 'Position/1', 'Position scope/1', 'autorouting', 'on');
catch err
    warning('Position scope not added (%s); the rest of the model is complete.', err.message);
end

save_system(mdl, fullfile(here, [mdl '.slx']));
fprintf('Built %s.slx with layout %s (%d propellers, thrust %g to %g N)\n', ...
        mdl, layout_name, N, min(P.u_lo), max(P.u_hi));
end
