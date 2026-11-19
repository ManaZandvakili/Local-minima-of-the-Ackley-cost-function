function P = airship_load(layout_name)
%AIRSHIP_LOAD  Parameter struct for the simulator: physical model + one layout.
%   P = airship_load('knee14_reversible')
%
%   Layouts available (data/layouts/*.mat):
%     knee14_reversible, knee8_reversible_capped, knee15_oneway,
%     knee9_oneway_capped, paper16_reversible, paper16_oneway,
%     starlab6_reversible
%   The numbers are exported from the Python code by scripts/export_matlab.py.

if nargin < 1
    layout_name = 'knee14_reversible';
end
here = fileparts(mfilename('fullpath'));
S = load(fullfile(here, 'data', 'airship_params.mat'));
L = load(fullfile(here, 'data', 'layouts', [layout_name '.mat']));
P = S.phys;
lay = L.layout;
P.A = lay.A;                       % 6 x N: thrust of each propeller -> forces and torques
P.u_lo = lay.u_lo(:);              % minimum thrust [N] (0 one-way, -12 reversible)
P.u_hi = lay.u_hi(:);              % maximum thrust [N]
P.positions_m = lay.positions_m;
P.directions = lay.directions;
P.n_props = size(lay.A, 2);
P.r_g = P.r_g(:);
P.d_lin = P.d_lin(:);
P.d_quad = P.d_quad(:);
P.tau_dist = P.tau_dist(:);
P.n_sub = round(P.n_sub);          % physics steps per control step (10)
P = rmfield(P, 'notes');           % Simulink parameters must be numeric
end
