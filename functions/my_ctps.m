function [V, pK] = my_ctps(ECGdata, ECGepochs)
% https://sci-hub.st/https://ieeexplore.ieee.org/document/4536072
% For CA detection, a threshold of pK ≥ 20 was applied.
% ICs with a CTPS just above this threshold may be barely detectable by visual inspection.

[NCHN, NTPT] = size(ECGdata);
assert(NTPT>NCHN);

% Channels are columns
ECGdata = ECGdata';

% Use the Hilbert transform to obtain the instantaneous phase:
inst_phase = angle(hilbert(ECGdata)); % [rad]

% Normalise the instantaneous phase into the range [0,1]:
norm_phase = mod(inst_phase / (2 * pi), 1);

V = NaN(NCHN,1);
pK = NaN(NCHN,1);

for i = 1:NCHN
    % For each R-peak, you can extract 1-second windows around them from the phase data:
    phase_windows = NaN(length(ECGepochs(1,1):ECGepochs(1,2)),size(ECGepochs,1));
    for j = 1:size(ECGepochs,1)
        thisChunk = ECGepochs(j,1):ECGepochs(j,2);
        phase_windows(:, j) = norm_phase(thisChunk,i);
    end

    % Sort the phases
    sorted_phases = sort(phase_windows(:));

    % The empirical cumulative distribution function (ECDF), K
    ecdf = sorted_phases;

    % The uniform CDF for comparison, G
    n = length(sorted_phases);
    uniform_cdf = (1:n)' / n;

    % Kuiper's statistic calculation
    D_plus = max(ecdf - uniform_cdf);
    D_minus = max(uniform_cdf - ecdf);

    % Kuiper statistic
    V(i) = D_plus + D_minus;

    % Compute p-value (this is a simplified approximation)
    % p(i) = exp(-2 * n * V(i)^2);

    % Calculate lambda (λ)
    lambda = V(i) * (sqrt(n) + 0.155 + 0.24 / sqrt(n));

    % Compute the p-value approximation PK using the infinite series
    PK = 0;
    for k = 1:100  % Summing a reasonable number of terms in the infinite series
        PK = PK + (4*k^2 * lambda^2 - 1) * exp(-2 * k^2 * lambda^2);
    end
    PK = 2*PK;  % Multiply by 2 according to the formula

    % Calculate the negative logarithmic p-value
    if PK > 0
        pK(i) = -log10(PK);  % pK = -log10(PK)
    else
        pK(i) = Inf;  % If PK becomes very small (practically zero)
    end

    % % Optional: Plot ECDF vs Uniform CDF
    % figure;
    % plot(sorted_phases, ecdf, 'b', sorted_phases, uniform_cdf, 'r--');
    % legend('Uniform CDF','Empirical CDF');
    % title('ECDF vs Uniform CDF');
end

end