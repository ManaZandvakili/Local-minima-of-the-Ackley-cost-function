clc; clear; close all;

% Ackley function setup
d = 2; % Number of dimensions
numParticles = 50; % Number of particles for PSO
numIterations = 100; % Number of iterations for PSO
tStart = 99; tEnd = 98.9; % Simulated annealing temperature
numCoolingLoops = 7500; numEquilbriumLoops = 50;
frac = (tEnd / tStart)^(1.0 / (numCoolingLoops - 1.0));
step_size = 3; % Perturbation step size

% Ackley function definition
ackley = @(x) -20 * exp(-0.2 * sqrt(sum(x.^2) / d)) - exp(sum(cos(2 * pi * x)) / d) + 20 + exp(1);

% PSO Initialization
particle_positions = -32.768 + (32.768 - (-32.768)) * rand(numParticles, d);
particle_velocities = zeros(numParticles, d);
personal_best_positions = particle_positions;
personal_best_costs = arrayfun(@(idx) ackley(particle_positions(idx, :)), 1:numParticles);
[global_best_cost, gbest_idx] = min(personal_best_costs);
global_best_position = personal_best_positions(gbest_idx, :);

% SA Initialization
x_current = global_best_position; % Start SA from PSO's best position
x_best = x_current; cost_best = ackley(x_best);
minima = {x_best};
distinct_minima = minima; % To store distinct local minima
threshold = 0.8; % Threshold for distinct minima

% Visualization grid for Ackley function
[X, Y] = meshgrid(-32.768:0.05:32.768, -32.768:0.05:32.768);
Z = arrayfun(@(x, y) ackley([x, y]), X, Y);

% Real-time visualization setup
figure(1);

% PSO + SA Hybrid Optimization
tCurrent = tStart;
for i = 1:numCoolingLoops
    % PSO Step
    for iter = 1:numIterations
        for p = 1:numParticles
            % Update particle velocities
            inertia = 0.7;
            self = 2 * rand() * (personal_best_positions(p, :) - particle_positions(p, :));
            social = 2 * rand() * (global_best_position - particle_positions(p, :));
            particle_velocities(p, :) = inertia * particle_velocities(p, :) + self + social;

            % Update particle positions
            particle_positions(p, :) = particle_positions(p, :) + particle_velocities(p, :);
            particle_positions(p, :) = max(min(particle_positions(p, :), 32.768), -32.768);

            % Evaluate cost
            cost = ackley(particle_positions(p, :));
            if cost < personal_best_costs(p)
                personal_best_costs(p) = cost;
                personal_best_positions(p, :) = particle_positions(p, :);
            end
        end

        % Update global best
        [current_global_best_cost, gbest_idx] = min(personal_best_costs);
        if current_global_best_cost < global_best_cost
            global_best_cost = current_global_best_cost;
            global_best_position = personal_best_positions(gbest_idx, :);
        end
    end

    % Simulated Annealing Step
    for j = 1:numEquilbriumLoops
        % Perturb solution
        x_new = x_current + step_size * (2 * rand(1, d) - 1);
        x_new = max(min(x_new, 32.768), -32.768);

        % Compute costs
        cost_current = ackley(x_current);
        cost_new = ackley(x_new);

        % Accept new solution
        delta_cost = cost_new - cost_current;
        if (delta_cost < 0) || (exp(-delta_cost / tCurrent) > rand())
            x_current = x_new;
            cost_current = cost_new;
        end

        % Check if new local minimum
        if cost_current < cost_best
            x_best = x_current;
            cost_best = cost_current;
        end
        
        % Check distinct minima using Euclidean distance
        is_distinct = true;
        for m = 1:length(distinct_minima)
            % Calculate Euclidean distance between current point and existing minima
            distance = norm(x_current - distinct_minima{m});
            if distance < threshold
                is_distinct = false;
                break;
            end
        end
        if is_distinct
            distinct_minima{end + 1} = x_current;
        end

        % Real-time visualization
        if mod(j, 10) == 0 % Update plot every 10 iterations
            minima_coords = cell2mat(minima');
            distinct_coords = cell2mat(distinct_minima');
            clf;
            surf(X, Y, Z, 'EdgeColor', 'none'); % 3D surface plot
            colormap jet; colorbar;
            hold on;
            scatter3(x_current(1), x_current(2), ackley(x_current), 50, 'b', 'filled'); % Current point
            scatter3(minima_coords(:, 1), minima_coords(:, 2), ...
                     arrayfun(@(k) ackley(minima_coords(k, :)), 1:size(minima_coords, 1)), ...
                     40, 'r', 'filled'); % Minima found
            scatter3(distinct_coords(:, 1), distinct_coords(:, 2), ...
                     arrayfun(@(k) ackley(distinct_coords(k, :)), 1:size(distinct_coords, 1)), ...
                     40, 'g', 'filled', 'MarkerEdgeColor', 'k','LineWidth', 2); % Distinct minima
            title(['Ackley Function | Cooling Loop: ', num2str(i), ...
                   ' | Temp: ', num2str(tCurrent), ...
                   ' | Distinct Minima: ', num2str(length(distinct_minima))]);
            xlabel('x_1'); ylabel('x_2'); zlabel('Cost');
            view(135, 45); % Adjust viewing angle
            grid on;
            drawnow;
        end
    end

    % Update temperature
    tCurrent = tCurrent * frac;

    % Restart mechanism for SA
    if mod(i, 100) == 0
        tCurrent = tStart; % Reset temperature
        x_current = global_best_position; % Restart from PSO's best solution
    end
end

% Final count
distinct_coords = cell2mat(distinct_minima');
disp(['Final Number of Distinct Minima (green dots): ', num2str(size(distinct_coords, 1))]);
disp(['Final Number of Global Best (red dots): 1']);
