function animate_airship(X, t, P, video_file, speed)
%ANIMATE_AIRSHIP  3-D animation of a flight, from run_airship or the Simulink log.
%   animate_airship(X, t, P)                     % play in a figure (2x speed)
%   animate_airship(X, t, P, 'flight.mp4')       % also save it as a video
%   animate_airship(X, t, P, '', 1)              % play in real time
%
%   X  states, one row per time step (13 + N columns), t their times [s],
%   P  parameters from airship_load (gives the propeller positions/directions).
%   Shows the hull turning with the airship, the gondola, every propeller with
%   an arrow for its thrust (orange = pushing, blue = pulling), the path so far
%   and the target (+ at the origin).

if nargin < 4, video_file = ''; end
if nargin < 5 || isempty(speed), speed = 2; end
t = t(:);
nx = size(X, 2);
N = nx - 13;
R0 = 1.575;                                      % hull radius [m]
gondola_body = [0; 0; -1.825];                   % gondola position, body frame [m]
arrow_m_per_N = 0.08;                            % arrow length per newton of thrust
stride = max(1, round(0.1 / max(t(2) - t(1), eps)));   % one frame every 0.1 s

[sx, sy, sz] = sphere(16);
S = R0 * [sx(:) sy(:) sz(:)].';

fig = figure('Name', 'Airship flight', 'Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
lim = 4;
xlim(ax, [-lim lim]); ylim(ax, [-lim lim]); zlim(ax, [-lim lim]);
view(ax, 35, 22);
xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]'); zlabel(ax, 'z [m]');
plot3(ax, 0, 0, 0, 'k+', 'MarkerSize', 14, 'LineWidth', 2);
hull = surf(ax, reshape(S(1, :), size(sx)), reshape(S(2, :), size(sx)), reshape(S(3, :), size(sx)), ...
            'FaceColor', [0.75 0.85 1], 'FaceAlpha', 0.15, 'EdgeColor', [0.45 0.55 0.75], 'EdgeAlpha', 0.35);
trail = plot3(ax, NaN, NaN, NaN, '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.2);
gond = plot3(ax, NaN, NaN, NaN, 's', 'MarkerSize', 9, 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'k');
rope = plot3(ax, NaN, NaN, NaN, 'k-');
props = plot3(ax, NaN, NaN, NaN, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
push = quiver3(ax, zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 1), 0, ...
               'Color', [0.90 0.45 0.10], 'LineWidth', 1.6, 'MaxHeadSize', 0.6);
pull = quiver3(ax, zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 1), 0, ...
               'Color', [0.15 0.45 0.85], 'LineWidth', 1.6, 'MaxHeadSize', 0.6);
ttl = title(ax, '');

v = [];
if ~isempty(video_file)
    v = VideoWriter(video_file, 'MPEG-4');
    v.FrameRate = max(1, round(speed / (stride * (t(2) - t(1)))));
    open(v);
end

frames = 1:stride:size(X, 1);
for k = frames
    p = X(k, 1:3).';
    R = quat_to_rot(X(k, 4:7));
    W = R * S + p;
    set(hull, 'XData', reshape(W(1, :), size(sx)), 'YData', reshape(W(2, :), size(sx)), ...
              'ZData', reshape(W(3, :), size(sx)));
    set(trail, 'XData', X(1:k, 1), 'YData', X(1:k, 2), 'ZData', X(1:k, 3));
    g = R * gondola_body + p;
    b = R * [0; 0; -R0] + p;
    set(gond, 'XData', g(1), 'YData', g(2), 'ZData', g(3));
    set(rope, 'XData', [b(1) g(1)], 'YData', [b(2) g(2)], 'ZData', [b(3) g(3)]);
    Pw = (R * P.positions_m.').' + p.';
    Dw = (R * P.directions.').';
    u = X(k, 14:end).';
    A = Dw .* (u * arrow_m_per_N);
    pos = u >= 0;
    set(props, 'XData', Pw(:, 1), 'YData', Pw(:, 2), 'ZData', Pw(:, 3));
    set(push, 'XData', Pw(pos, 1), 'YData', Pw(pos, 2), 'ZData', Pw(pos, 3), ...
              'UData', A(pos, 1), 'VData', A(pos, 2), 'WData', A(pos, 3));
    set(pull, 'XData', Pw(~pos, 1), 'YData', Pw(~pos, 2), 'ZData', Pw(~pos, 3), ...
              'UData', A(~pos, 1), 'VData', A(~pos, 2), 'WData', A(~pos, 3));
    set(ttl, 'String', sprintf('t = %4.1f s    distance from target %.2f m', t(k), norm(p)));
    drawnow;
    if ~isempty(v)
        writeVideo(v, getframe(fig));
    else
        pause(stride * (t(2) - t(1)) / speed);
    end
    if ~ishandle(fig), break; end
end
if ~isempty(v)
    close(v);
    fprintf('Saved %s\n', video_file);
end
end
