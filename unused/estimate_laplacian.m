function laplac = estimate_laplacian(data, chanlocs, stepSize)

% Number of channels (electrodes)
n = length(data(:,1));

% Preallocate spherical coordinate matrix
a = zeros(n, 2);

% Extract spherical coordinates (assuming they are in degrees)
for i = 1:n
    a(i, 1) = chanlocs(i).sph_phi;
    a(i, 2) = chanlocs(i).sph_theta;
end

% Convert to radians
aux = a * pi / 180;

% Initialise matrices for g function and constants
g = zeros(n);
matriz_g = zeros(n + 1);

% Calculate g and g2 matrices for each channel
for i = 1:n
    matriz_g(i, 1:n) = g2(aux(i, 1), aux(i, 2), aux(:, 1), aux(:, 2));
    g(i, :) = gm_1(aux(i, 1), aux(i, 2), aux(:, 1), aux(:, 2));
end

% Boundary condition: set the last row and column to 1
matriz_g(:, n + 1) = 1;
matriz_g(n + 1, :) = 1;

% Precompute the inverse of matriz_g (you could replace this with a linear solver)
inv_m = inv(matriz_g);

% Preallocate for the Laplacian results
largo = length(data(1, :));
laplac = zeros(n, max(1,round((largo - 1) / stepSize)));

% Loop through time points with the specified stepSize
k = 1;
for j = 1:stepSize:largo
    % Get data at time step j
    V = data(:, j);

    % Calculate constants using the inverse matrix
    constant = inv_m * [V; 0];  % Append 0 to V for the boundary condition

    % Compute the Laplacian for each channel
    for i = 1:n
        laplac(i, k) = g(i, :) * constant(1:n);
    end

    % Increment time step
    k = k + 1;
end

end

%% ========================================================================

function valor = g2(r1, r2, r3, r4)
% Compute g2 function based on spherical harmonics, summing to n=25

% Compute x from the spherical coordinates
x = cos(r1) .* cos(r2) .* cos(r3) .* cos(r4) + ...
    cos(r1) .* sin(r2) .* cos(r3) .* sin(r4) + ...
    sin(r1) .* sin(r3);

% Clip values of x to be between -1 and 1
x = min(max(x, -1), 1);

% Preallocate for Legendre polynomials
P = zeros(25, length(x));

% Compute Legendre polynomials up to order n=25
for n = 1:25
    aux = legendre(n, x);
    P(n, :) = aux(1, :);  % Use only the first order term
end

% Calculate the g2 function using a weighted sum
n = 1:25;
aux2 = ((2 * n + 1) ./ (n .* (n + 1)).^3)';
valor = 1 / (4 * pi) * sum(aux2 .* P, 1);
end

function valor = gm_1(r1, r2, r3, r4)
% Compute gm_1 function, similar to g2 but with different weighting

% Both param: Smaller -> smoother
smoothnesParam = 2;
K = 25;

% Compute x from spherical coordinates
x = cos(r1) .* cos(r2) .* cos(r3) .* cos(r4) + ...
    cos(r1) .* sin(r2) .* cos(r3) .* sin(r4) + ...
    sin(r1) .* sin(r3);

% Clip x to ensure it is between -1 and 1
x = min(max(x, -1), 1);

% Preallocate for Legendre polynomials
L = length(x);
P = zeros(K, L);

% Compute Legendre polynomials up to n=25
for n = 1:K
    aux = legendre(n, x);
    P(n, :) = aux(1, :);
end

% Calculate gm_1 using a different weighting for each term
n = 1:K;
aux2 = (2 * n + 1) ./ (n .* (n + 1)).^smoothnesParam;
aux2 = repmat(aux2',1,L);
valor = 1 / (4 * pi) * sum(aux2 .* P, 1);

end