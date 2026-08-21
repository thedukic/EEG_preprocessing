function [V, pK] = my_ctps(ECGdata, ECGepochs, event, fs)
% https://ieeexplore.ieee.org/document/4536072/
% For CA detection, a threshold of pK >= 20 was applied.

[NCHN, NTPT] = size(ECGdata);
assert(NTPT > NCHN);

% Filtering
% Paper: 10-20 Hz filter
[bh, ah] = butter(2, 10/(fs/2), 'high');
[bl, al] = butter(2, 20/(fs/2), 'low');
ECGdata = do_filteringcore(bl, al, ECGdata, event, fs);
ECGdata = do_filteringcore(bh, ah, ECGdata, event, fs);

% Transpose: channels -> columns
ECGdata = ECGdata';

% 1. Calculate continuous instantaneous phase
inst_phase = angle(hilbert(ECGdata)); % [rad]

% 2. Normalize to [0, 1] using modulo as defined in the paper
norm_phase = mod(inst_phase / (2 * pi), 1);

V = NaN(NCHN, 1);
pK = NaN(NCHN, 1);

% Loop through channels
for i_channel = 1:NCHN

    num_samples = length(ECGepochs(1, 1):ECGepochs(1, 2));
    num_trials = size(ECGepochs, 1);

    % 3. Split into time windows (trials)
    phase_windows = NaN(num_samples, num_trials);
    for i_epoch = 1:num_trials
        indx = ECGepochs(i_epoch, 1):ECGepochs(i_epoch, 2);
        phase_windows(:, i_epoch) = norm_phase(indx, i_channel);
    end

    V_time = NaN(num_samples, 1);
    pK_time = NaN(num_samples, 1);

    % 4. Calculate distributions for EACH time point relative to the onset
    for t = 1:num_samples
        % Extract the cross-trial phases for this specific time sample
        phases_at_t = phase_windows(t, :)';

        % Sort the phases for the Empirical CDF
        sorted_phases = sort(phases_at_t);

        step_idx = (1:num_trials)';

        % Kuiper's statistic calculation
        D_plus = max(step_idx/num_trials - sorted_phases);
        D_minus = max(sorted_phases - (step_idx - 1)/num_trials);

        % Kuiper statistic V for this time point
        V_time(t) = D_plus + D_minus;

        % Calculate lambda
        lambda = V_time(t) * (sqrt(num_trials) + 0.155 + 0.24 / sqrt(num_trials));

        % Compute the p-value approximation PK
        PK = 0;
        for k = 1:100
            PK = PK + (4*k^2 * lambda^2 - 1) * exp(-2 * k^2 * lambda^2);
        end
        PK = 2 * PK;

        % Calculate the negative logarithmic p-value
        if PK > 0
            pK_time(t) = -log10(PK);
        else
            pK_time(t) = Inf;
        end
    end

    % Deal with inf values
    if ~all(isfinite(pK_time))
        pK_time(~isfinite(pK_time)) = max(pK_time(isfinite(pK_time)));
    end

    % Extract the maximum significance found within the epoch window
    V(i_channel) = mean(V_time);
    pK(i_channel) = mean(pK_time);
end

[max_val, max_indx] = max(pK);
fprintf('Maximum detected pK = %1.1f with V = %1.2f.\n', max_val, V(max_indx));

end

% figure; hold on;
% plot(V_time);
% plot(isfinite(pK_time));