clc; clear; close all;

% Ackley function setup
d = 2; % Number of dimensions
numCoolingLoops = 5000000; % Number of cooling loops for SA
numEquilbriumLoops = 100;  % Number of steps at each temperature
tStart = 99; tEnd = 98.9; % Simulated annealing temperature
frac = (tEnd / tStart)^(1.0 / (numCoolingLoops - 1));
step_size = 3; % Perturbation step size
threshold = 0.8; % Threshold for distinct minima
Offset = 0.1;

% Ackley function definition
ackley = @(x) -20 * exp(-0.2 * sqrt(sum(x.^2) / numel(x))) - exp(sum(cos(2 * pi * x)) / numel(x)) + 20 + exp(1);

% Camparing algorithm with the neighboring points
Camparing_algo = @(point) ...
    ackley(max(min(point + Offset * [1, 0], 32.768), -32.768)) >= ackley(point) && ...
    ackley(max(min(point + Offset * [-1, 0], 32.768), -32.768)) >= ackley(point) && ...
    ackley(max(min(point + Offset * [0, 1], 32.768), -32.768)) >= ackley(point) && ...
    ackley(max(min(point + Offset * [0, -1], 32.768), -32.768)) >= ackley(point) && ...
    ackley(max(min(point + Offset * [1, 1], 32.768), -32.768)) >= ackley(point) && ...
    ackley(max(min(point + Offset * [-1, -1], 32.768), -32.768)) >= ackley(point) && ...
    ackley(max(min(point + Offset * [1, -1], 32.768), -32.768)) >= ackley(point) && ...
    ackley(max(min(point + Offset * [-1, 1], 32.768), -32.768)) >= ackley(point);

% Euclidean Distance Calculator
Euclidean_Cal = @(found_minima, newSolution, threshold) ...
    all(sqrt(sum((found_minima - newSolution).^2, 2)) >= threshold);

% Simulated Annealing Initialization
x_current = -32.768 + (32.768 - (-32.768)) * rand(1, d); % Start SA from a random initial position
x_best = x_current;
cost_best = ackley(x_best);
found_minima = []; % storing
cost_current = ackley(x_current);
numAcceptedSolutions = 1.0;
DeltaE_avg = 0.0;

% Plot the Ackley Function Surface
[X, Y] = meshgrid(-32.768:0.05:32.768, -32.768:0.05:32.768);
Z = arrayfun(@(x, y) ackley([x, y]), X, Y);

figure;
surf(X, Y, Z, 'EdgeColor', 'none');
colormap jet;
hold on;
xlabel('X1'); ylabel('X2'); zlabel('Cost');
title(['All Distinct Local Minima on Ackley Function']);
view(3);
% Storing the local Minima
minimaPlot = scatter3([], [], [], 40, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Simulated Annealing Optimization Loop
tCurrent = tStart;
for i = 1:numCoolingLoops
    % Output the details for each iteration
    disp(['Cooling Loop: ', num2str(i), ...
          ', Equilibrium Loop: ', num2str(j), ...
          ', Temperature: ', num2str(tCurrent), ...
          ', Current Cost: ', num2str(cost_current), ...
          ', Best Cost: ', num2str(cost_best), ...
          ', Number of Local Minima: ', num2str(size(found_minima, 1))]);

    for j = 1:numEquilbriumLoops
        % Perturb solution
        x_new = x_current + step_size * (2 * rand(1, d) - 1);
        x_new = max(min(x_new, 32.768), -32.768);

        % Compute costs
        cost_current = ackley(x_current);
        cost_new = ackley(x_new);
        
        % Accept new solution
        delta_cost = cost_new - cost_current;
        if delta_cost < 0 
            accept = true;
        else 
            if (i == 1 && j == 1)
                DeltaE_avg = abs(delta_cost);
            end
            p = exp(-delta_cost / (DeltaE_avg * tCurrent));
            accept = (p > rand());
        end

        if accept
            % Update current solution
            x_current = x_new;
            cost_current = cost_new;
            numAcceptedSolutions = numAcceptedSolutions + 1.0;
            DeltaE_avg = (DeltaE_avg * (numAcceptedSolutions - 1.0) + abs(delta_cost)) / numAcceptedSolutions;
        end
    end

    % Checking the distinctness
    if Camparing_algo(x_current)
        if isempty(found_minima) || Euclidean_Cal(found_minima, x_current, threshold)
            found_minima = [found_minima; x_current];
            set(minimaPlot, 'SizeData', 20);
            set(minimaPlot, 'XData', found_minima(:, 1), 'YData', found_minima(:, 2), ...
                            'ZData', arrayfun(@(idx) ackley(found_minima(idx, :)), 1:size(found_minima, 1)));
            drawnow;

            % Restart mechanism: Generate a new random starting point
            x_current = -32.768 + (32.768 - (-32.768)) * rand(1, d);
            cost_current = ackley(x_current);
        end
    end

    tCurrent = frac * tCurrent; % Reduce temperature
end

% Display Results
disp(['Number of Local Minima Found: ', num2str(size(found_minima, 1))]);
