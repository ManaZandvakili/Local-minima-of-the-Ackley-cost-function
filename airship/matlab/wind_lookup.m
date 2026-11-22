function wind_block = wind_lookup(t, WIND)
%WIND_LOOKUP  Wind for the 10 physics steps of the control step that starts at time t (30 x 1).
%   WIND.table: one row per control step, [w1 w2 ... w10] (3 numbers each), from closed_loop_setup.
%#codegen

k = round(t / WIND.dt) + 1;
k = min(max(k, 1), size(WIND.table, 1));
wind_block = WIND.table(k, :).';
end
