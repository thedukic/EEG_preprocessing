function EEG = do_cca(EEG, cfg)
% DO_CCA Performs CCA on an EEGLAB dataset.
%
% Inputs:
%   EEG - EEGLAB structure containing continuous or epoched data.
%
% Outputs:
%   W   - Unmixing matrix.
%   S   - Canonical source activations (components).
%   rho - Autocorrelation values for each component.
%   fh  - Figure handle for the debugging visualisations.

fprintf('\n================================\n');
fprintf('CCA: Canonical correlation analysis (time-lagged)\n');
fprintf('================================\n');

% Check if data is EEG or EMG
is_eeg = all(lower(string({EEG.chanlocs.type})) == "eeg");

% Extract and reshape data to 2D (Channels x Time)
if is_eeg
    X = double(reshape(EEG.data, EEG.nbchan, []));
else
    [X, clean_idx, bad_idx, quality_table, clean_epochs] = prepare_clean_emg_cca_data(EEG);
end

% Remove the mean of each row
X = X - mean(X, 2);

[n_chan, n_pnts] = size(X);

% =========================================================================
% CCA
% =========================================================================
% Create the predefined temporally delayed function y(t) = x(t-1)
X0 = X(:, 1:end-1);
X1 = X(:, 2:end);

% 3. Calculate within-set and between-sets covariance matrices
Cxx = (X0 * X0') / (n_pnts - 1);
Cyy = (X1 * X1') / (n_pnts - 1);
Cxy = (X0 * X1') / (n_pnts - 1);
% Cyx = Cxy';

% Solve the eigenvalue problem
% The paper defines the problem as: (Cxx^-1 * Cxy * Cyy^-1 * Cyx) * W = rho^2 * W
% Using pinv() ensures numerical stability for rank-deficient EEG data
% =========================================================================
% Dynamic Regularisation & Rank Handling
% =========================================================================
data_rank = rank(Cxx);
cond_num  = cond(Cxx);

if data_rank < n_chan
    % Exact rank deficiency (e.g. average referenced, interpolated channels)
    % Truncate to the true subspace using PCA
    fprintf('Rank deficiency detected (Rank: %d / %d). Applying PCA reduction...\n', data_rank, n_chan);

    [U, S_val, ~] = svd(Cxx, 'econ');
    P = U(:, 1:data_rank); % Subspace projector

    X0_sub = P' * X0;
    X1_sub = P' * X1;

    Cxx_sub = (X0_sub * X0_sub') / (n_pnts - 1);
    Cyy_sub = (X1_sub * X1_sub') / (n_pnts - 1);
    Cxy_sub = (X0_sub * X1_sub') / (n_pnts - 1);

    M_sub = (Cxx_sub \ Cxy_sub) * (Cyy_sub \ Cxy_sub');
    [W_sub, D] = eig(M_sub);

    % Back-project weights to full sensor space
    W = P * W_sub;

elseif cond_num > 1e6
    % Full rank but ill-conditioned: Apply shrinkage
    fprintf('Ill-conditioned covariance (cond: %.2e). Applying shrinkage...\n', cond_num);
    gamma_reg = 0.01;
    I = eye(n_chan);

    Cxx_reg = (1 - gamma_reg) * Cxx + gamma_reg * (trace(Cxx) / n_chan) * I;
    Cyy_reg = (1 - gamma_reg) * Cyy + gamma_reg * (trace(Cyy) / n_chan) * I;

    M = (Cxx_reg \ Cxy) * (Cyy_reg \ Cxy');
    [W, D] = eig(M);

else
    % Well-conditioned: Direct GEVP solution
    fprintf('Covariance matrices are good.\n');
    M = (Cxx \ Cxy) * (Cyy \ Cxy');
    [W, D] = eig(M);
end

% Extract the eigenvalues (rho^2)
rho_sq = diag(D);

% Prevent small negative values due to numerical noise and take the square root
rho = sqrt(max(rho_sq, 0));

% Sort components by autocorrelation
[rho, sort_idx] = sort(rho, 'descend');
W = W(:, sort_idx);

% Compute the K estimates of the sources using z(t) = W^T * x(t)
X = double(reshape(EEG.data, EEG.nbchan, []));
S = W' * X;

% Update the variable because of possible rank reduction
n_chan = length(rho);
n_pnts = size(X, 2);

% =========================================================================
% Remove bad components
% =========================================================================
if is_eeg
    % Calculate the Power Spectral Density (PSD) for all components in S
    fs = EEG.srate;
    win_len = 2 * fs;

    S_unit = S ./ std(S, 0, 2);
    [psd_S, freq] = pwelch(S_unit', hann(win_len), win_len/2, [], fs);
    psd_S = psd_S'; % Transpose to Components x Frequencies

    % Define frequency bands for the ratio
    idx_low  = (freq >= 3 & freq <= 15);
    idx_high = (freq >= 30 & freq <= min(100, fs/2));

    % Calculate the Spectral Power Ratio for each component
    power_low  = mean(psd_S(:, idx_low), 2);
    power_high = mean(psd_S(:, idx_high), 2);
    spectral_ratio = power_high ./ power_low;

    % -------------------------------------------------------------------------
    % Simple
    % -------------------------------------------------------------------------
    % Define automated thresholds
    autocorr_threshold = 0.5; % Components with rho below this are flagged
    ratio_threshold    = 0.5; % Components with a ratio above this are flagged

    % % Find components that fail BOTH criteria (low autocorrelation AND high spectral ratio)
    % bad_components = find((rho < autocorr_threshold) & (spectral_ratio > ratio_threshold));

    % -------------------------------------------------------------------------
    % Automated Threshold Detection via Curve Gap & Knee Point
    % -------------------------------------------------------------------------
    % 1. Min-Max normalise both metrics to [0, 1]
    norm_minmax = @(v) (v - min(v)) ./ (max(v) - min(v) + eps);
    rho_norm    = norm_minmax(rho);
    ratio_norm  = norm_minmax(spectral_ratio);

    % 2. Calculate the Gap Curve: Delta = rho_norm - ratio_norm
    gap_curve = rho_norm - ratio_norm;

    % Smooth the gap curve slightly to prevent local jitter from triggering false knees
    gap_smooth = movmean(gap_curve, 5);

    % -------------------------------------------------------------------------
    % Method A: Zero-Crossing / Level-Crossing of the Gap
    % -------------------------------------------------------------------------
    % Components where the spectral ratio exceeds the autocorrelation
    bad_idx_gap = find(gap_smooth <= 0.5); % Can adjust to 0.1 or 0.0

    if isempty(bad_idx_gap)
        min_bad_idx_gap = NaN;
    else
        min_bad_idx_gap = min(bad_idx_gap);
    end

    % -------------------------------------------------------------------------
    % Method B: Maximum Distance / Knee Point (Vectorised Kneedle approach)
    % -------------------------------------------------------------------------
    % Form a chord between the first and last points of the gap curve
    n_points = length(gap_smooth);
    chord_x = [1, n_points];
    chord_y = [gap_smooth(1), gap_smooth(end)];

    % Compute perpendicular distance of each point to the chord
    p1 = [1, gap_smooth(1)];
    p2 = [n_points, gap_smooth(end)];
    diff_vec = p2 - p1;

    dist_to_chord = zeros(n_points, 1);
    for i = 1:n_points
        p0 = [i, gap_smooth(i)];
        dist_to_chord(i) = abs(det([diff_vec; p0 - p1])) / norm(diff_vec);
    end

    % The maximum distance is the transition knee
    [~, knee_idx] = max(dist_to_chord);
    bad_idx_knee = (knee_idx:n_points)';

    fprintf('Method A (Gap <= 0.5): %d components to remove.\n', length(bad_idx_gap));
    fprintf('Method B (Knee Point): %d components to remove (from comp %d onwards).\n', ...
        length(bad_idx_knee), knee_idx);

    % -------------------------------------------------------------------------
    % Visualisation: Diagnostic Gap & Knee Curve
    % -------------------------------------------------------------------------
    fh = figure('Color', 'w', 'Position', [200, 200, 700, 450]);
    plot(1:n_chan, rho_norm, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Normalised \rho'); hold on;
    plot(1:n_chan, ratio_norm, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Normalised Spectral Ratio');
    plot(1:n_chan, gap_smooth, 'k--', 'LineWidth', 2, 'DisplayName', 'Smoothed Gap (\Delta)');

    % Highlight the knee
    xline(knee_idx, 'g-', sprintf('Knee Point (Comp %d)', knee_idx), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

    if ~isnan(min_bad_idx_gap)
        xline(min_bad_idx_gap, 'g-', sprintf('Gap point (Comp %d)', min_bad_idx_gap), ...
            'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    end

    yline(0, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'HandleVisibility', 'off');

    xlabel('CCA Component Index');
    ylabel('Normalised Value [0, 1]');
    title('CCA Transition & Knee-Point Detection');
    legend('Location', 'southwest');
    grid on; box off;

    % Save
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_cca_1'], [15 15/1.6]);

    % -------------------------------------------------------------------------
    % Power spectra properties
    % -------------------------------------------------------------------------
    % neg_integral = estimate_spectral_metrics(psd_S, freq, [1 70], [2 40]);
    % threshold_power_int = 4; % empirical for raw EEG signal: 10
    % bad_kurt = find(neg_integral > threshold_power_int);

    % -------------------------------------------------------------------------
    % Kurtosis
    % -------------------------------------------------------------------------
    % 1. Spatial Focalisation Metric (Gini Index or Maximum-to-Mean Ratio)
    % True channel pops concentrate almost all energy in 1-2 electrodes
    A = pinv(W');
    A_abs = abs(A);

    % Spatial Kurtosis across electrodes
    spat_kurt = kurtosis(A_abs, 1, 1)';

    % Max-to-Mean absolute ratio (focal spike index)
    max_to_mean = (max(A_abs, [], 1) ./ mean(A_abs, 1))';

    % % Ultra-Low Frequency Power Concentration (< 1 Hz vs 1-30 Hz)
    % % Wobbles/drift concentrate massive power in sub-delta (<1 Hz)
    % idx_drift = (freq >= 0.1 & freq < 1.0);
    % idx_neural = (freq >= 1.0 & freq <= 30.0);
    %
    % p_drift  = mean(psd_S(:, idx_drift), 2);
    % p_neural = mean(psd_S(:, idx_neural), 2);
    % drift_ratio = p_drift ./ (p_neural + eps);

    % Composite Decision Rule:
    % Flag if: (Extremely focal pop) OR (Massive <1Hz drift with high spatial peak)
    is_focal_pop  = (spat_kurt > 25) & (max_to_mean > 6);
    % is_slow_drift = (drift_ratio > 10) & (spat_kurt > 8);

    % bad_kurt = find(is_focal_pop | is_slow_drift);
    bad_kurt = find(is_focal_pop);

    fprintf('Identified %d low-frequency wobble/drift component(s) at the head.\n', length(bad_kurt));

    % -------------------------------------------------------------------------
    % Do it
    % -------------------------------------------------------------------------
    % EMG CCA are not auto-correlated
    bad_components = bad_idx_knee;
    mask_bad_def   = rho(bad_components) < autocorr_threshold;
    bad_components = bad_components(mask_bad_def);

    % Slow drifts that are not removed by lowpass filtering are auto-correlated
    % mask_bad_def   = rho(bad_kurt) > autocorr_threshold;
    % bad_kurt       = bad_kurt(mask_bad_def);
    % bad_kurt = bad_kurt(bad_kurt <= knee_idx);
    bad_kurt = bad_kurt(rho(bad_kurt) > 0.8);

    bad_components = unique([bad_components(:); bad_kurt(:)]);

    fprintf('Automatically identified %d muscle components for removal.\n', length(bad_components));

    % Zero out the identified muscle components
    if ~isempty(bad_components)
        S_clean = S;
        S_clean(bad_components, :) = 0;

        % Reconstruct the EEG data
        EEG_old = EEG;
        A_clean = pinv(W');
        X_clean = A_clean * S_clean;
        EEG.data = reshape(X_clean, size(EEG.data));
    else
        fprintf('CCA removal will not be done.\n');
        EEG_old = EEG;
    end

    % =========================================================================
    % Visualisations
    % =========================================================================
    % -------------------------------------------------------------------------
    % Visualise 1
    % -------------------------------------------------------------------------
    % Check pwoer change
    check_power_cleaning(EEG_old, EEG, [], 'cca', cfg);

    % -------------------------------------------------------------------------
    % Visualise 2
    % -------------------------------------------------------------------------
    num_plot = 5; % Number of tail components to plot in rows 2 & 3
    num_plot = min(num_plot, n_chan);

    fh = figure('Color', 'w', 'Position', [100, 100, max(900, num_plot * 260), 900], 'Name', 'CCA Metrics & Artefact QA');
    t = tiledlayout(3, num_plot, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Top-Left Panel: Dual-Axis Component Plot ---
    % Spans the left half of the top row
    nexttile(1, [1 floor(num_plot/2) + (mod(num_plot,2) ~= 0)]);
    yyaxis left
    plot(1:n_chan, rho, '-o', 'Color', [0.2 0.4 0.7], 'LineWidth', 1.5, 'MarkerSize', 4);
    ylabel('Autocorrelation (\rho)');
    ylim([0 1.05]);

    yyaxis right
    plot(1:n_chan, spectral_ratio, '-s', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5, 'MarkerSize', 4);
    ylabel('Spectral Ratio (20-100Hz / 1-15Hz)');

    title('Autocorrelation & Spectral Ratio per Component', 'FontWeight', 'bold');
    xlabel('CCA Component (Sorted by \rho)');
    xlim([1 n_chan]);
    grid on;

    % --- Top-Right Panel: Scatter Plot (rho vs Spectral Ratio) ---
    % Spans the remaining columns of the top row
    nexttile(floor(num_plot/2) + 1 + (mod(num_plot,2) ~= 0), [1 ceil(num_plot/2) - (mod(num_plot,2) ~= 0)]);
    scatter(rho, spectral_ratio, 35, 1:n_chan, 'filled');
    colormap(gca, 'parula');
    cb = colorbar;
    cb.Label.String = 'Component Index';
    xlabel('Autocorrelation (\rho)');
    ylabel('Spectral Ratio');
    title('Covariation: \rho vs Spectral Ratio', 'FontWeight', 'bold');
    grid on;

    % Highlight the lowest components in the scatter plot
    hold on;
    scatter(rho(end-num_plot+1:end), spectral_ratio(end-num_plot+1:end), 70, 'r', 'LineWidth', 1.5);

    % --- Row 2: Time Courses of the Tail Components ---
    time_vec = (0:(min(n_pnts, 5000)-1)) / EEG.srate;

    for i = 0:(num_plot - 1)
        comp_idx = n_chan - i;
        nexttile(num_plot + i + 1, [1 1]);
        plot(time_vec, S(comp_idx, 1:length(time_vec)), 'k');
        title({sprintf('Comp %d', comp_idx), ...
            sprintf('\\rho=%.2f, Ratio=%.2f', rho(comp_idx), spectral_ratio(comp_idx))}, ...
            'FontSize',8);
        xlabel('Time (s)');
        pbaspect([1.618 1 1]); box off;
    end

    % --- Row 3: Field Distributions (Topographies) ---
    for i = 0:(num_plot - 1)
        comp_idx = n_chan - i;
        % nexttile(2*num_plot + i + 1, [1 1]);
        % topoplot(A(:, comp_idx), EEG.chanlocs, 'electrodes', 'on');
        % title(sprintf('Map %d', comp_idx));
        mytopoplot(A(:, comp_idx), [], sprintf('Map %d', comp_idx), nexttile(2*num_plot + i + 1, [1 1]));
    end

    % Save
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_cca_2'], [30 16]);

    % -------------------------------------------------------------------------
    % Visualise 3
    % -------------------------------------------------------------------------
    num_remove = length(bad_components);
    bad_components = flip(bad_components);

    if num_remove == 0
        fprintf('No bad CCA components identified to plot.\n');
    else
        % Fixed 8-column layout (or fewer if total components < 8)
        n_cols = min(8, num_remove);
        n_rows = ceil(num_remove / n_cols);

        % Create figure
        fh = figure('Color', 'w');
        tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Plot all topographies dynamically
        for i = 1:num_remove
            comp_idx = bad_components(i);
            mytopoplot(A(:, comp_idx), [], sprintf('Map %d', comp_idx), nexttile);
        end

        % Dynamically scale saved figure dimensions (width x height) based on row/column count
        fig_width  = max(12, n_cols * 3.5);
        fig_height = max(6,  n_rows * 3.5);

        % Save
        save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_cca_overview'], [fig_width, fig_height]);
    end

elseif ~is_eeg
    % Check
    comp_idx = 1;
    fh = diagnose_emg_component(EEG, W, S, comp_idx);

    % Save
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_cca_overview'], [fig_width, fig_height]);

end


% =========================================================================
% Diagnose specific components
% =========================================================================
% for i = 55
%     diagnose_eeg_component(EEG, W, S, rho, i, [20, 30]);
% end

end

% % function [neg_integral, slope_low, slope_discrepancy] =
% estimate_spectral_metrics(powspctrm, freqs, interest_range, censor_range)
% - already a fucntion file
% % % ESTIMATE_SPECTRAL_METRICS Computes slope discrepancy and negative integral metrics.
% % %
% % % Inputs:
% % %   powspctrm      : Matrix of power spectra (num_channels x freqs)
% % %   freqs          : Vector of frequencies (1 x freqs or freqs x 1)
% % %   interest_range : 1x2 vector of frequencies to include (e.g., [1 70])
% % %   censor_range   : 1x2 vector of frequencies to exclude (e.g., [3 30])
% % %
% % % Outputs:
% % %   slope_discrepancy : Vector (num_channels x 1) of internal consistency differences
% % %   neg_integral      : Vector (num_channels x 1) of the total negative residual area
% % 
% % % Ensure column vector for frequencies
% % freqs = freqs(:);
% % num_channels = size(powspctrm, 1);
% % 
% % log_freqs = log10(freqs);
% % log_power_spectra = log10(powspctrm);
% % 
% % % --- FREQUENCY MASKING ---
% % interest_idx1 = freqs >= interest_range(1) & freqs <= interest_range(2);
% % interest_idx2 = freqs < censor_range(1)    | freqs > censor_range(2);
% % interest_idx3 = freqs < 49 | freqs > 51;
% % interest_idx4 = freqs < 99 | freqs > 101;
% % interest_idx  = interest_idx1 & interest_idx2 & interest_idx3 & interest_idx4;
% % 
% % freqs_interest = freqs(interest_idx);
% % log_freqs_select = log_freqs(interest_idx);
% % log_power_select = log_power_spectra(:, interest_idx);
% % 
% % % Identify the low/high frequency 'islands' within the filtered data
% % idx_low  = find(freqs_interest <= censor_range(1));
% % idx_high = find(freqs_interest >= censor_range(2));
% % 
% % % Fallback if data is entirely outside the censored range
% % if isempty(idx_low) && ~isempty(idx_high)
% %     idx_low  = find(freqs_interest <= freqs_interest(15));
% %     idx_high = find(freqs_interest >= freqs_interest(end-15));
% % end
% % 
% % % --- ALLOCATION & MODEL FITTING ---
% % slope_low = NaN(num_channels, 1);
% % slope_discrepancy = NaN(num_channels, 1);
% % ap_fit = NaN(num_channels, length(freqs));
% % 
% % for i_channel = 1:num_channels
% %     power_channel = log_power_select(i_channel, :);
% % 
% %     % Global fit to calculate the aperiodic baseline
% %     mdl = fitlm(log_freqs_select, power_channel(:));
% %     intercept = mdl.Coefficients.Estimate(1);
% %     slope     = mdl.Coefficients.Estimate(2);
% % 
% %     ap_fit_log = intercept + (slope * log_freqs');
% %     ap_fit(i_channel, :) = 10.^ap_fit_log;
% % 
% %     % Sub-segment fits for internal consistency
% %     mdl_low  = fitlm(log_freqs_select(idx_low), power_channel(idx_low));
% %     mdl_high = fitlm(log_freqs_select(idx_high), power_channel(idx_high));
% % 
% %     slope_low(i_channel, 1) = mdl_low.Coefficients.Estimate(2);
% %     slope_discrepancy(i_channel, 1) = abs(mdl_low.Coefficients.Estimate(2) - mdl_high.Coefficients.Estimate(2));
% % end
% % 
% % % --- NEGATIVE RESIDUAL ESTIMATION ---
% % periodic_estimate = log_power_spectra - log10(ap_fit);
% % periodic_estimate = periodic_estimate(:, interest_idx1);
% % neg_mask          = periodic_estimate < 0;
% % 
% % % Area of 'impossible' power per channel
% % neg_integral = sum(abs(periodic_estimate .* neg_mask), 2);
% % end

function fh = diagnose_eeg_component(EEG, W, S, rho, comp_idx, time_window)
% DIAGNOSE_CCA_COMPONENT Generates a comprehensive QA diagnostic panel
% for a specific CCA component.
%
% Inputs:
%   EEG         - EEGLAB data structure
%   W           - Unmixing matrix (n_chan x n_chan)
%   S           - Component activation matrix (n_chan x n_samples)
%   rho         - Vector of autocorrelation values (sorted, n_chan x 1)
%   comp_idx    - Index of the component to diagnose (1 to n_chan)
%   time_window - [Optional] [start_sec, end_sec] for time plot (default: [0, 5])

if nargin < 6 || isempty(time_window)
    time_window = [0, min(5, EEG.pnts / EEG.srate)];
end

n_chan = size(S, 1);
if comp_idx < 1 || comp_idx > n_chan
    error('comp_idx must be between 1 and %d.', n_chan);
end

fs = EEG.srate;
comp_signal = double(S(comp_idx, :));

% ---------------------------------------------------------------------
% 1. Compute Metrics (PSD, Spectral Ratio, Peak Frequency)
% ---------------------------------------------------------------------
win_len = min(4 * fs, length(comp_signal));
[psd_val, freqs] = pwelch(comp_signal, hann(win_len), win_len/2, 1024, fs);

idx_low  = (freqs >= 1 & freqs <= 15);
idx_high = (freqs >= 20 & freqs <= min(100, fs/2));

p_low  = mean(psd_val(idx_low));
p_high = mean(psd_val(idx_high));
spec_ratio = p_high / (p_low + eps);

[~, max_idx] = max(psd_val(freqs >= 1 & freqs <= 45));
valid_freqs = freqs(freqs >= 1 & freqs <= 45);
peak_freq = valid_freqs(max_idx);

% Forward model for scalp topography (inverse of W')
A = pinv(W');
topo_weights = A(:, comp_idx);

% ---------------------------------------------------------------------
% 2. Create Diagnostic Visualisation Layout
% ---------------------------------------------------------------------
fh = figure('Color', 'w', 'Position', [150, 150, 1100, 750], ...
    'Name', sprintf('CCA Component %d Diagnostics', comp_idx));
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel 1: Scalp Topography ---
nexttile;
if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
    topoplot(topo_weights, EEG.chanlocs, 'electrodes', 'on', 'style', 'both');
    title(sprintf('Topography (Comp %d)', comp_idx), 'FontWeight', 'bold');
    colorbar;
else
    axis off;
    text(0.5, 0.5, 'No chanlocs found', 'HorizontalAlignment', 'center');
end

% --- Panel 2: Power Spectral Density (PSD) ---
nexttile;
plot(freqs, 10 * log10(psd_val), 'k-', 'LineWidth', 1.5);
hold on;
xline(peak_freq, 'b--', sprintf('Peak: %.1f Hz', peak_freq), ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');
xlim([1, min(100, fs/2)]);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power / Frequency (dB/Hz)');
title('Power Spectrum (1-100 Hz)', 'FontWeight', 'bold');
box off;

% --- Panel 3: Quantitative Summary Card ---
nexttile;
axis off;

% Classification heuristic for quick visual guidance
if rho(comp_idx) < 0.5 && spec_ratio > 0.6
    classification_str = 'Likely Muscle / High-Freq Noise';
    status_color = [0.8 0.1 0.1];
elseif rho(comp_idx) >= 0.7 && spec_ratio < 0.4
    classification_str = 'Likely Cortical / Neural Signal';
    status_color = [0.1 0.6 0.2];
else
    classification_str = 'Mixed / Review Required';
    status_color = [0.85 0.5 0.1];
end

summary_text = {
    sprintf('\\bfComponent Index:\\rm %d of %d', comp_idx, n_chan);
    '';
    sprintf('\\bfAutocorrelation (\\rho):\\rm %.3f', rho(comp_idx));
    sprintf('\\bfSpectral Ratio (20-100Hz / 1-15Hz):\\rm %.3f', spec_ratio);
    sprintf('\\bfDominant Low Peak:\\rm %.1f Hz', peak_freq);
    '';
    sprintf('\\bfDiagnostic Suggestion:\\rm');
    sprintf('\\color[rgb]{%.2f,%.2f,%.2f}\\bf%s', ...
    status_color(1), status_color(2), status_color(3), classification_str)
    };

text(0.1, 0.5, summary_text, 'FontSize', 11, 'Interpreter', 'tex');

% --- Panel 4: Time Series Activation Snippet ---
nexttile([1, 2]);
p_start = max(1, round(time_window(1) * fs));
p_end   = min(length(comp_signal), round(time_window(2) * fs));
t_axis  = (p_start:p_end) / fs;

plot(t_axis, comp_signal(p_start:p_end), 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Activation (a.u.)');
title(sprintf('Time Series Snapshot (%.1f s - %.1f s)', time_window(1), time_window(2)), ...
    'FontWeight', 'bold');
xlim([time_window(1), time_window(2)]);
grid on;
box off;

% --- Panel 5: Spectrogram (Time-Frequency) ---
nexttile;
nfft = 256;
noverlap = round(nfft * 0.8);
spectrogram(comp_signal, kaiser(nfft, 5), noverlap, nfft, fs, 'yaxis');
ylim([0, min(80, fs/2)]);
colormap(gca, 'parula');
title('Time-Frequency Spectrogram', 'FontWeight', 'bold');
ylabel('Frequency (Hz)');
xlabel('Time (s)');
colorbar off;
end

function fh = diagnose_emg_component(EEG, W, S, comp_idx)
% =========================================================================
% Trigger-Locked BSS-CCA EMG Diagnostic (Multi-Condition: 20s, 30s, 50s)
% Time Window: -3 s before Start Press to +3 s after Stop Press
% Includes raw time-series trace (first 10 s) at the bottom.
% =========================================================================
fs = EEG.srate;

% Channel labels
if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs)
    chan_labels = {EEG.chanlocs.labels};
else
    chan_labels = arrayfun(@(x) sprintf('Ch %d', x), 1:size(W, 1), 'UniformOutput', false);
end

% -------------------------------------------------------------------------
% 1. Define Multi-Condition Markers
% -------------------------------------------------------------------------
markers_onset = {'condition 20', 'condition 30', 'condition 50'};
markers_start = {'condition 21', 'condition 31', 'condition 51'};
markers_stop  = {'condition 22', 'condition 32', 'condition 52'};

% Normalise event types to lower-case string array for robust matching
all_types = lower(strtrim(string({EEG.event.type})));
all_lat   = [EEG.event.latency];

% Find all start-pressing events across all conditions
is_start_event = ismember(all_types, markers_start);
start_indices  = find(is_start_event);

n_trials = length(start_indices);
if n_trials == 0
    error('No events found matching start markers (condition 21, 31, 51).');
end

comp_signal = double(S(comp_idx, :));

% Pre-allocate trial metadata
trial_cond     = strings(n_trials, 1);
lat_start      = NaN(n_trials, 1);
stop_durations = NaN(n_trials, 1);
start_base_dur = NaN(n_trials, 1);

for tr = 1:n_trials
    idx_tr = start_indices(tr);
    t_start = all_lat(idx_tr);
    lat_start(tr) = t_start;
    trial_cond(tr) = all_types(idx_tr);

    % Match corresponding Stop marker (must occur after current Start)
    stop_cand_idx = find(ismember(all_types, markers_stop) & (all_lat >= t_start));
    if ~isempty(stop_cand_idx)
        stop_durations(tr) = (all_lat(stop_cand_idx(1)) - t_start) / fs;
    end

    % Match corresponding Onset marker (must occur before current Start)
    onset_cand_idx = find(ismember(all_types, markers_onset) & (all_lat <= t_start));
    if ~isempty(onset_cand_idx)
        start_base_dur(tr) = (all_lat(onset_cand_idx(end)) - t_start) / fs;
    end
end

% -------------------------------------------------------------------------
% 2. Epoch Construction (-3 s before Start to +3 s after Max Stop)
% -------------------------------------------------------------------------
t_pre = 3.0;
max_stop_dur = nanmax(stop_durations);
if isnan(max_stop_dur)
    t_post = 6.0;
else
    t_post = max_stop_dur + 3.0;
end

win_samples = round(-t_pre * fs) : round(t_post * fs);
t_axis = win_samples / fs;

trial_matrix = NaN(n_trials, length(win_samples));
for tr = 1:n_trials
    t_start = lat_start(tr);
    lat_range = round(t_start + win_samples);

    valid_mask = (lat_range >= 1) & (lat_range <= length(comp_signal));
    if any(valid_mask)
        trial_matrix(tr, valid_mask) = comp_signal(lat_range(valid_mask));
    end
end

% -------------------------------------------------------------------------
% 3. Calculate Rectified Envelope & Condition Metrics
% -------------------------------------------------------------------------
trial_rect = abs(trial_matrix);
trial_env  = movmean(trial_rect, round(0.05 * fs), 2); % 50 ms smoothing

% Grand average
mean_env = nanmean(trial_env, 1);
sem_env  = nanstd(trial_env, [], 1) ./ sqrt(sum(~isnan(trial_env(:, 1))));

avg_stop_time  = nanmean(stop_durations);
avg_start_base = nanmean(start_base_dur);

% -------------------------------------------------------------------------
% 4. Multi-Panel Diagnostic Visualisation Layout (4 Rows x 2 Columns)
% -------------------------------------------------------------------------
fh = figure('Color', 'w', 'Position', [100, 50, 1200, 950], ...
    'Name', sprintf('EMG BSS-CCA Comp %d [Conds: 20/30 Series]', comp_idx));
t = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel 1: Spatial Projection Vector (A) ---
nexttile(1, [1 1]);
A = pinv(W');
spatial_weights = A(:, comp_idx);
bar(spatial_weights, 'FaceColor', [0.2 0.45 0.7]);
set(gca, 'XTick', 1:length(chan_labels), 'XTickLabel', chan_labels);
xtickangle(45);
ylabel('Weight (a.u.)');
title(sprintf('Spatial Projection (Comp %d)', comp_idx), 'FontWeight', 'bold');
grid on; box off;

% --- Panel 2: Power Spectral Density ---
nexttile(2, [1 1]);
win_len = min(2 * fs, length(comp_signal));
[psd_val, freqs] = pwelch(comp_signal, hann(win_len), win_len/2, 2048, fs);
plot(freqs, 10 * log10(psd_val), 'r', 'LineWidth', 1.5);
xlim([0, min(500, fs/2)]);
xlabel('Frequency (Hz)');
ylabel('Power (dB/Hz)');
title('Power Spectral Density', 'FontWeight', 'bold');
grid on; box off;

% --- Panel 3: Trial-by-Trial Raster Heatmap ---
nexttile(3, [1 2]);
imagesc(t_axis, 1:n_trials, trial_env);
colormap(gca, 'hot');
cb = colorbar;
cb.Label.String = 'Envelope (a.u.)';

hold on;
xline(0, 'c-', 'Start Press (21/31)', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
if ~isnan(avg_start_base)
    xline(avg_start_base, 'c--', 'Avg Onset (20/30)', 'LineWidth', 1.2, 'LabelOrientation', 'horizontal');
end

% Overlay exact stop markers per trial
valid_stops = ~isnan(stop_durations);
plot(stop_durations(valid_stops), find(valid_stops), 'c.', 'MarkerSize', 9);

xlim([t_axis(1), t_axis(end)]);
xlabel('Time Relative to Start Press (s)');
ylabel('Trial #');
title('Trial-by-Trial EMG Envelope (-3 s before Start to +3 s after Stop)', 'FontWeight', 'bold');
box off;

% --- Panel 4: Mean Contraction Profile (Overall + Per Condition) ---
nexttile(5, [1 2]);
fill([t_axis, fliplr(t_axis)], [mean_env + sem_env, fliplr(mean_env - sem_env)], ...
    [0.2 0.4 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Overall SEM');
hold on;
plot(t_axis, mean_env, 'Color', [0.1 0.2 0.7], 'LineWidth', 2.2, 'DisplayName', 'Overall Mean');

% Plot condition-specific averages if multiple conditions exist
cond_palette = [0.85 0.33 0.1; 0.47 0.67 0.19; 0.49 0.18 0.56];
unique_conds = unique(trial_cond);

for c = 1:length(unique_conds)
    c_mask = (trial_cond == unique_conds(c));
    if sum(c_mask) > 1
        c_mean = nanmean(trial_env(c_mask, :), 1);
        col = cond_palette(mod(c-1, size(cond_palette, 1)) + 1, :);
        plot(t_axis, c_mean, 'Color', col, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Mean %s (n=%d)', unique_conds(c), sum(c_mask)));
    end
end

xline(0, 'g-', 'Start Press', 'LineWidth', 1.5, 'HandleVisibility', 'off');
if ~isnan(avg_stop_time)
    xline(avg_stop_time, 'r--', sprintf('Avg Stop: %.2f s', avg_stop_time), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    % Highlight active contraction window
    patch([0 avg_stop_time avg_stop_time 0], ...
        [min(ylim) min(ylim) max(ylim) max(ylim)], ...
        [0.3 0.8 0.3], 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

xlim([t_axis(1), t_axis(end)]);
xlabel('Time Relative to Start Press (s)');
ylabel('Activation (a.u.)');
title('Averaged Contraction Profile', 'FontWeight', 'bold');
grid on; box off;
legend('Location', 'northeast');

% --- Panel 5: Raw Continuous Component Trace (First 10 s) ---
nexttile(7, [1 2]);
raw_dur_s = 15;
n_raw_samples = min(round(raw_dur_s * fs), length(comp_signal));
t_raw = (0:(n_raw_samples - 1)) / fs;

plot(t_raw, comp_signal(1:n_raw_samples), 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
hold on;

% Mark any events that fall within this first 10-second window
events_in_win = find(all_lat >= 1 & all_lat <= n_raw_samples);
for ev_idx = events_in_win
    ev_time = (all_lat(ev_idx) - 1) / fs;
    ev_label = string(EEG.event(ev_idx).type);
    xline(ev_time, 'm--', ev_label, 'LabelVerticalAlignment', 'top', ...
        'LabelOrientation', 'horizontal', 'LineWidth', 1, 'HandleVisibility', 'off');
end

xlim([0, t_raw(end)]);
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');
title(sprintf('Raw Component Time Series (First %.1f s)', t_raw(end)), 'FontWeight', 'bold');
grid on; box off;

end

function [X, clean_idx, bad_idx, quality_table, clean_epochs] = prepare_clean_emg_cca_data(EEG, varargin)
% =========================================================================
% PREPARE_CLEAN_EMG_CCA_DATA (Adaptive Threshold Version)
% Uses data-driven relative outlier detection (Tukey IQR) tailored for
% monopolar surface EMG with baseline ECG and submaximal contractions.
% =========================================================================

% --- 1. Parse Input Options with Calibrated Defaults ---
p = inputParser;
addParameter(p, 'PreTime', 3.0, @isnumeric);
addParameter(p, 'PostBufferTime', 3.0, @isnumeric);
addParameter(p, 'BaselineWindow', [-2.5, -0.5], @(x) isnumeric(x) && length(x)==2);

% Tolerant activation ratio for submaximal / monopolar tasks
addParameter(p, 'MinActivationRatio', 1.2, @isnumeric);

% Multipliers for relative IQR outlier rejection (Tukey bounds)
addParameter(p, 'IQRMultiplierBaseNoise', 2.5, @isnumeric);
addParameter(p, 'IQRMultiplierBaseKurt',  3.0, @isnumeric);
addParameter(p, 'IQRMultiplierPeakAmp',   3.0, @isnumeric);

% Optional hard lower bound on whole-epoch kurtosis (catches flatlines / saturation)
addParameter(p, 'MinEpochKurtosis', 1.5, @isnumeric);
addParameter(p, 'PlotDiagnostics', true, @islogical);
parse(p, varargin{:});

t_pre         = p.Results.PreTime;
t_post_buf    = p.Results.PostBufferTime;
base_win      = p.Results.BaselineWindow;
min_act_ratio = p.Results.MinActivationRatio;
iqr_k_base    = p.Results.IQRMultiplierBaseNoise;
iqr_k_kurt    = p.Results.IQRMultiplierBaseKurt;
iqr_k_peak    = p.Results.IQRMultiplierPeakAmp;
min_ep_kurt   = p.Results.MinEpochKurtosis;
do_plot       = p.Results.PlotDiagnostics;

fs     = EEG.srate;
n_chan = EEG.nbchan;

% --- 2. Identify Event Markers ---
markers_onset = {'condition 20', 'condition 30', 'condition 50'};
markers_start = {'condition 21', 'condition 31', 'condition 51'};
markers_stop  = {'condition 22', 'condition 32', 'condition 52'};

all_types = lower(strtrim(string({EEG.event.type})));
all_lat   = [EEG.event.latency];

start_indices = find(ismember(all_types, markers_start));
n_trials      = length(start_indices);

if n_trials == 0
    error('No events found matching Start markers (condition 21, 31, 51).');
end

% --- 3. Match Latencies & Durations ---
lat_start      = NaN(n_trials, 1);
stop_durations = NaN(n_trials, 1);
trial_cond     = strings(n_trials, 1);

for tr = 1:n_trials
    idx_tr = start_indices(tr);
    t_s = all_lat(idx_tr);
    lat_start(tr) = t_s;
    trial_cond(tr) = all_types(idx_tr);

    stop_cand = find(ismember(all_types, markers_stop) & (all_lat >= t_s));
    if ~isempty(stop_cand)
        stop_durations(tr) = (all_lat(stop_cand(1)) - t_s) / fs;
    end
end

max_stop = nanmax(stop_durations);
if isnan(max_stop)
    max_stop = 3.0;
end
t_post = max_stop + t_post_buf;

win_samples = round(-t_pre * fs) : round(t_post * fs);
t_axis      = win_samples / fs;
n_samples   = length(win_samples);

% --- 4. Extract Multi-Channel Epochs ---
raw_epochs = NaN(n_chan, n_samples, n_trials);
data_cont  = double(EEG.data);

for tr = 1:n_trials
    t_s = lat_start(tr);
    lat_range = round(t_s + win_samples);

    valid_mask = (lat_range >= 1) & (lat_range <= size(data_cont, 2));
    if any(valid_mask)
        raw_epochs(:, valid_mask, tr) = data_cont(:, lat_range(valid_mask));
    end
end

% --- 5. Extract Per-Trial Metrics ---
base_mask = (t_axis >= base_win(1)) & (t_axis <= base_win(2));

rms_base_global   = NaN(n_trials, 1);
rms_act_global    = NaN(n_trials, 1);
max_base_kurtosis = NaN(n_trials, 1);
min_epoch_kurt_tr = NaN(n_trials, 1);
peak_amp_global   = NaN(n_trials, 1);

for tr = 1:n_trials
    trial_data = raw_epochs(:, :, tr);
    if all(isnan(trial_data(:)))
        continue;
    end

    t_stop = stop_durations(tr);
    if isnan(t_stop) || t_stop < 0.5
        t_stop = max_stop;
    end
    act_mask = (t_axis >= 0.2) & (t_axis <= (t_stop - 0.1));

    % 1. Global RMS Energy
    env_tr = sqrt(nanmean(trial_data.^2, 1));
    rms_base_global(tr) = sqrt(nanmean(env_tr(base_mask).^2));
    rms_act_global(tr)  = sqrt(nanmean(env_tr(act_mask).^2));

    % 2. Baseline Kurtosis across channels
    base_data = trial_data(:, base_mask);
    kurt_base_ch = NaN(n_chan, 1);
    for ch = 1:n_chan
        sig_b = base_data(ch, :);
        valid_b = sig_b(~isnan(sig_b));
        if length(valid_b) > 10
            kurt_base_ch(ch) = kurtosis(valid_b);
        end
    end
    max_base_kurtosis(tr) = nanmax(kurt_base_ch);

    % 3. Whole-Epoch Kurtosis (check for flatlines/saturation)
    kurt_epoch_ch = NaN(n_chan, 1);
    for ch = 1:n_chan
        sig_e = trial_data(ch, :);
        valid_e = sig_e(~isnan(sig_e));
        if length(valid_e) > 10
            kurt_epoch_ch(ch) = kurtosis(valid_e);
        end
    end
    min_epoch_kurt_tr(tr) = nanmin(kurt_epoch_ch);

    % 4. Absolute Peak across channels
    peak_amp_global(tr) = nanmax(abs(trial_data(:)));
end

act_ratio = rms_act_global ./ (rms_base_global + eps);

% --- 6. Adaptive Relative Outlier Detection ---

% Relative Outlier Threshold 1: Baseline RMS (Tonic noise)
q75_b = prctile(rms_base_global, 75);
iqr_b = iqr(rms_base_global);
thresh_base_rms = q75_b + iqr_k_base * iqr_b;
flag_noisy_baseline = rms_base_global > thresh_base_rms;

% Relative Outlier Threshold 2: Baseline Kurtosis (Outlier pops only)
q75_k = prctile(max_base_kurtosis, 75);
iqr_k = iqr(max_base_kurtosis);
thresh_base_kurt = q75_k + iqr_k_kurt * iqr_k;
flag_baseline_pop = max_base_kurtosis > thresh_base_kurt;

% Relative Outlier Threshold 3: Peak Amplitude (Severe movement/jumps)
q75_p = prctile(peak_amp_global, 75);
iqr_p = iqr(peak_amp_global);
thresh_peak_amp = q75_p + iqr_k_peak * iqr_p;
flag_extreme_peak = peak_amp_global > thresh_peak_amp;

% Dynamic Activation Contrast
flag_poor_contrast = act_ratio < min_act_ratio;

% Flatlined / Saturated Signals
flag_low_kurtosis = min_epoch_kurt_tr < min_ep_kurt;

% Combine Flags
is_bad = flag_noisy_baseline | flag_baseline_pop | flag_poor_contrast | ...
    flag_extreme_peak | flag_low_kurtosis;

clean_idx = find(~is_bad);
bad_idx   = find(is_bad);

% Summary table
quality_table = table((1:n_trials)', trial_cond, act_ratio, rms_base_global, ...
    max_base_kurtosis, min_epoch_kurt_tr, peak_amp_global, ...
    flag_noisy_baseline, flag_baseline_pop, flag_poor_contrast, ...
    flag_extreme_peak, flag_low_kurtosis, ~is_bad, ...
    'VariableNames', {'Trial', 'Condition', 'ActivationRatio', 'BaselineRMS', ...
    'MaxBaseKurtosis', 'MinEpochKurtosis', 'PeakAmp', ...
    'Flag_NoisyBaseline', 'Flag_BasePop', 'Flag_PoorContrast', ...
    'Flag_ExtremePeak', 'Flag_LowKurtosis', 'KeepTrial'});

fprintf('\n================== EMG Quality Assessment Summary ==================\n');
fprintf('Total Trials Processed:  %d\n', n_trials);
fprintf('Retained Clean Trials:   %d (%.1f%%)\n', length(clean_idx), 100 * length(clean_idx) / n_trials);
fprintf('Rejected Bad Trials:     %d (%.1f%%)\n', length(bad_idx), 100 * length(bad_idx) / n_trials);
fprintf(' - Tonic Baseline Noise: %d (Thresh: %.2f)\n', sum(flag_noisy_baseline), thresh_base_rms);
fprintf(' - Baseline Spikes/Pops: %d (Thresh: %.2f)\n', sum(flag_baseline_pop), thresh_base_kurt);
fprintf(' - Poor Dynamic Contrast:%d (Thresh: %.2f)\n', sum(flag_poor_contrast), min_act_ratio);
fprintf(' - Extreme Peak Outliers:%d (Thresh: %.2f)\n', sum(flag_extreme_peak), thresh_peak_amp);
fprintf(' - Low Kurtosis/Clipping:%d (Thresh: %.2f)\n', sum(flag_low_kurtosis), min_ep_kurt);
fprintf('====================================================================\n\n');

% --- 7. Construct Concatenated Clean Matrix X ---
clean_epochs = raw_epochs(:, :, clean_idx);
X = double(reshape(clean_epochs, n_chan, []));
valid_cols = ~any(isnan(X), 1);
X = X(:, valid_cols);

% --- 8. Diagnostic Visualisation ---
if do_plot
    fh = figure('Color', 'w', 'Position', [80, 80, 1200, 750], 'Name', 'EMG Epoch Screening Diagnostics (Adaptive)');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Panel 1: Baseline RMS
    nexttile(1);
    bar(1:n_trials, rms_base_global, 'FaceColor', [0.3 0.65 0.4]); hold on;
    yline(thresh_base_rms, 'r--', sprintf('IQR Limit (%.2f)', thresh_base_rms), 'LineWidth', 1.5);
    plot(find(flag_noisy_baseline), rms_base_global(flag_noisy_baseline), 'rx', 'MarkerSize', 8, 'LineWidth', 1.8);
    xlabel('Trial #'); ylabel('RMS (a.u.)');
    title('Baseline Noise (Tonic)', 'FontWeight', 'bold'); grid on; box off;

    % Panel 2: Baseline Kurtosis
    nexttile(2);
    bar(1:n_trials, max_base_kurtosis, 'FaceColor', [0.8 0.4 0.3]); hold on;
    yline(thresh_base_kurt, 'r--', sprintf('Relative IQR Limit (%.1f)', thresh_base_kurt), 'LineWidth', 1.5);
    plot(find(flag_baseline_pop), max_base_kurtosis(flag_baseline_pop), 'rx', 'MarkerSize', 8, 'LineWidth', 1.8);
    xlabel('Trial #'); ylabel('Kurtosis');
    title('Baseline Kurtosis (Outliers Only)', 'FontWeight', 'bold'); grid on; box off;

    % Panel 3: Activation Ratio
    nexttile(3);
    bar(1:n_trials, act_ratio, 'FaceColor', [0.25 0.5 0.75]); hold on;
    yline(min_act_ratio, 'r--', sprintf('Min (%.2f)', min_act_ratio), 'LineWidth', 1.5);
    plot(find(flag_poor_contrast), act_ratio(flag_poor_contrast), 'rx', 'MarkerSize', 8, 'LineWidth', 1.8);
    xlabel('Trial #'); ylabel('Active / Base RMS Ratio');
    title('Activation Ratio (Contrast)', 'FontWeight', 'bold'); grid on; box off;

    % Panel 4: Extreme Peaks
    nexttile(4);
    bar(1:n_trials, peak_amp_global, 'FaceColor', [0.85 0.55 0.15]); hold on;
    yline(thresh_peak_amp, 'r--', sprintf('Spike Limit (%.1f)', thresh_peak_amp), 'LineWidth', 1.5);
    plot(find(flag_extreme_peak), peak_amp_global(flag_extreme_peak), 'rx', 'MarkerSize', 8, 'LineWidth', 1.8);
    xlabel('Trial #'); ylabel('Max Amplitude (a.u.)');
    title('Extreme Sensor Outliers', 'FontWeight', 'bold'); grid on; box off;

    % Panel 5: Whole-Epoch Min Kurtosis
    nexttile(5);
    bar(1:n_trials, min_epoch_kurt_tr, 'FaceColor', [0.55 0.35 0.75]); hold on;
    yline(min_ep_kurt, 'r--', sprintf('Min (%.1f)', min_ep_kurt), 'LineWidth', 1.5);
    plot(find(flag_low_kurtosis), min_epoch_kurt_tr(flag_low_kurtosis), 'rx', 'MarkerSize', 8, 'LineWidth', 1.8);
    xlabel('Trial #'); ylabel('Min Kurtosis across Channels');
    title('Channel Flatline / Saturation', 'FontWeight', 'bold'); grid on; box off;

    % Panel 6: Clean vs Rejected Waveforms
    nexttile(6);
    global_env_all = squeeze(sqrt(nanmean(raw_epochs.^2, 1)))';
    if ~isempty(clean_idx)
        plot(t_axis, nanmean(global_env_all(clean_idx, :), 1), 'Color', [0.1 0.35 0.75], ...
            'LineWidth', 2.0, 'DisplayName', sprintf('Clean (n=%d)', length(clean_idx)));
        hold on;
    end
    if ~isempty(bad_idx)
        plot(t_axis, nanmean(global_env_all(bad_idx, :), 1), 'Color', [0.85 0.2 0.2], ...
            'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', sprintf('Rejected (n=%d)', length(bad_idx)));
    end
    xline(0, 'k:', 'Start Press');
    xlabel('Time Relative to Start (s)'); ylabel('Global RMS (a.u.)');
    title('Clean vs Rejected Waveforms', 'FontWeight', 'bold');
    grid on; box off; legend('Location', 'northeast');
end

end