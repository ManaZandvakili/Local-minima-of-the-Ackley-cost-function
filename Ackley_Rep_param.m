% Define the Ackley function
function y = ackley(xx, a, b, c)
    d = length(xx);
    if nargin < 4
        c = 2 * pi;
    end
    if nargin < 3
        b = 0.2;
    end
    if nargin < 2
        a = 20;
    end

    sum1 = 0;
    sum2 = 0;
    for ii = 1:d
        xi = xx(ii);
        sum1 = sum1 + xi^2;
        sum2 = sum2 + cos(c * xi);
    end

    term1 = -a * exp(-b * sqrt(sum1 / d));
    term2 = -exp(sum2 / d);
    y = term1 + term2 + a + exp(1);
end

% Visualization
clc; clear; close all;

% Define the range for x1 and x2
x1 = linspace(-32.768, 32.768, 100);
x2 = linspace(-32.768, 32.768, 100);
[X1, X2] = meshgrid(-32.768:0.05:32.768, -32.768:0.05:32.768);
Z = arrayfun(@(x1, x2) ackley([x1, x2]), X1, X2);

% Compute the Ackley function values
Z = zeros(size(X1));
for i = 1:size(X1, 1)
    for j = 1:size(X1, 2)
        Z(i, j) = ackley([X1(i, j), X2(i, j)]);
    end
end

% Plot the Ackley function
figure;
surf(X1, X2, Z, 'EdgeColor', 'none');
colormap jet;
colorbar;
title('Ackley Function');
xlabel('x_1');
ylabel('x_2');
zlabel('f(x_1, x_2)');
view(-45, 45);
