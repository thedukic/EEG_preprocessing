function mu = estimate_huber_mean(X, delta)
if nargin < 2, delta = 1.345; end % Default standard tuning constant

% Initial guesses across columns
mu = median(X, 1);
tol = 1e-6;
max_iter = 100;

for iter = 1:max_iter
    mu_old = mu;
    residuals = X - mu; % Vectorised broadcasting subtraction

    % Robust scale estimation using MAD along channels
    sigma = median(abs(residuals), 1) / 0.6745;
    sigma(sigma == 0) = 1e-6; % Prevent division by zero vectorially

    % Calculate Huber weights vectorially
    scaled_res = abs(residuals ./ sigma);
    weights = ones(size(X));

    % Use logical indexing to update outliers across the entire matrix at once
    mask = scaled_res > delta;
    weights(mask) = delta ./ scaled_res(mask);

    % Weighted mean update along dimension 1
    mu = sum(weights .* X, 1) ./ sum(weights, 1);

    % Check global convergence across all columns
    if max(abs(mu - mu_old)) < tol
        break;
    end
end
end