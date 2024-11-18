function low_res = estimate_invlaplacian(data, chanlocs, stepSize)

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

% Initialise matrices for low-resolution spatial smoothing
smoothing_matrix = zeros(n);

% Calculate the smoothing matrix for each channel based on gm_1 (or g2 if preferred)
for i = 1:n
    smoothing_matrix(i, :) = gm_1(aux(i, 1), aux(i, 2), aux(:, 1), aux(:, 2));
end
smoothing_matrix = abs(smoothing_matrix);

% Preallocate for the low-resolution results
largo = length(data(1, :));
low_res = zeros(n, max(1,round((largo - 1) / stepSize)));

% Loop through time points with the specified stepSize
k = 1;
for j = 1:stepSize:largo
    % Get data at time step j
    V = data(:, j);

    % Apply smoothing to obtain low-resolution spatial components
    for i = 1:n
        w = smoothing_matrix(i, :);
        w = w ./ sum(w);
        low_res(i, k) = w * V;

        % low_res(i, k) = sum(smoothing_matrix(i, :) .* V) / sum(smoothing_matrix(i, :));
    end

    % Increment time step
    k = k + 1;
end

end

%% ========================================================================
function valor = gm_1(r1, r2, r3, r4)
% Compute gm_1 function, similar to g2 but with different weighting

% Both param: Smaller -> smoother
smoothnesParam = 0.1;
K = 15;

% Compute x from spherical coordinates
x = cos(r1) .* cos(r2) .* cos(r3) .* cos(r4) + ...
    cos(r1) .* sin(r2) .* cos(r3) .* sin(r4) + ...
    sin(r1) .* sin(r3);

% Clip x to ensure it is between -1 and 1
x = min(max(x, -1), 1);

% Preallocate for Legendre polynomials
L = length(x);
P = zeros(K, L);

% Compute Legendre polynomials up to n=K
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