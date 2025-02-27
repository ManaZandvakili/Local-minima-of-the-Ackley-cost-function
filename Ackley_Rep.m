function Ackley_Rep
    % The range for the grid
    s = -32.768:0.882:32.768;
    [x, y] = meshgrid(s, s);
    
    % Flatten the grid into column vectors
    X = x(:);
    Y = y(:);
    
    % Initialize the z vector
    z = zeros(length(X), 1);
    
    % Compute the Ackley function for each point
    for i = 1:length(X)
        z(i) = ackley([X(i), Y(i)]);
    end
    
    % Reshape z into a grid for plotting
    z = reshape(z, size(x));
    
    % Create the surface plot
    surf(x, y, z);
    
end
