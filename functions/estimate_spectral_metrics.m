function [neg_integral, slope_low, slope_discrepancy] = estimate_spectral_metrics(powspctrm, freqs, interest_range, censor_range, plot_indx)
% ESTIMATE_SPECTRAL_METRICS Computes slope discrepancy and negative integral metrics.
%
% Inputs:
%   powspctrm      : Matrix of power spectra (num_channels x freqs)
%   freqs          : Vector of frequencies (1 x freqs or freqs x 1)
%   interest_range : 1x2 vector of frequencies to include (e.g., [1 70])
%   censor_range   : 1x2 vector of frequencies to exclude (e.g., [3 30])
%   plot_indx      : (Optional) Indices of channels/components to visualise (e.g., [1, 5, 10])
%
% Outputs:
%   neg_integral      : Vector (num_channels x 1) of the total negative residual area
%   slope_low         : Vector (num_channels x 1) of low-frequency band slopes
%   slope_discrepancy : Vector (num_channels x 1) of internal consistency differences

if nargin < 5
    plot_indx = [];
end

% Ensure column vector for frequencies
freqs = freqs(:);
num_channels = size(powspctrm, 1);

log_freqs = log10(freqs);
log_power_spectra = log10(powspctrm);

% --- FREQUENCY MASKING ---
interest_idx1 = freqs >= interest_range(1) & freqs <= interest_range(2);
interest_idx2 = freqs < censor_range(1)    | freqs > censor_range(2);
interest_idx3 = freqs < 49 | freqs > 51;
interest_idx4 = freqs < 99 | freqs > 101;
interest_idx  = interest_idx1 & interest_idx2 & interest_idx3 & interest_idx4;

freqs_interest   = freqs(interest_idx);
log_freqs_select = log_freqs(interest_idx);
log_power_select = log_power_spectra(:, interest_idx);

% Identify the low/high frequency 'islands' within the filtered data
idx_low  = find(freqs_interest <= censor_range(1));
idx_high = find(freqs_interest >= censor_range(2));

% Fallback if data is entirely outside the censored range
if isempty(idx_low) && ~isempty(idx_high)
    idx_low  = find(freqs_interest <= freqs_interest(15));
    idx_high = find(freqs_interest >= freqs_interest(end-15));
end

% --- ALLOCATION & MODEL FITTING ---
slope_low         = NaN(num_channels, 1);
slope_discrepancy = NaN(num_channels, 1);
ap_fit            = NaN(num_channels, length(freqs));

% Storage for optional plotting lines
if ~isempty(plot_indx)
    fit_params = struct('int_glob', [], 'slope_glob', [], ...
        'int_low', [], 'slope_low', [], ...
        'int_high', [], 'slope_high', []);
    fit_params = repmat(fit_params, num_channels, 1);
end

for i_channel = 1:num_channels
    power_channel = log_power_select(i_channel, :);

    % Global fit to calculate the aperiodic baseline
    mdl = fitlm(log_freqs_select, power_channel(:));
    intercept = mdl.Coefficients.Estimate(1);
    slope     = mdl.Coefficients.Estimate(2);

    ap_fit_log = intercept + (slope * log_freqs');
    ap_fit(i_channel, :) = 10.^ap_fit_log;

    % Sub-segment fits for internal consistency
    mdl_low  = fitlm(log_freqs_select(idx_low), power_channel(idx_low));
    mdl_high = fitlm(log_freqs_select(idx_high), power_channel(idx_high));

    slope_low(i_channel, 1) = mdl_low.Coefficients.Estimate(2);
    slope_discrepancy(i_channel, 1) = abs(mdl_low.Coefficients.Estimate(2) - mdl_high.Coefficients.Estimate(2));

    if ~isempty(plot_indx) && ismember(i_channel, plot_indx)
        fit_params(i_channel).int_glob   = intercept;
        fit_params(i_channel).slope_glob = slope;
        fit_params(i_channel).int_low    = mdl_low.Coefficients.Estimate(1);
        fit_params(i_channel).slope_low  = mdl_low.Coefficients.Estimate(2);
        fit_params(i_channel).int_high   = mdl_high.Coefficients.Estimate(1);
        fit_params(i_channel).slope_high = mdl_high.Coefficients.Estimate(2);
    end
end

% --- NEGATIVE RESIDUAL ESTIMATION ---
periodic_estimate = log_power_spectra - log10(ap_fit);
periodic_estimate = periodic_estimate(:, interest_idx1);
neg_mask          = periodic_estimate < 0;

% Area of 'impossible' power per channel
neg_integral = sum(abs(periodic_estimate .* neg_mask), 2);

% --- VISUALISATION ---
if ~isempty(plot_indx)
    valid_plots = plot_indx(plot_indx >= 1 & plot_indx <= num_channels);
    n_plots     = length(valid_plots);

    if n_plots > 0
        cols = ceil(sqrt(n_plots));
        rows = ceil(n_plots / cols);

        fh = figure('Color', 'w', 'Name', 'Log-Log Spectral Fits', ...
            'Units', 'normalized', 'Position', [0.1, 0.1, min(0.85, 0.22*cols), min(0.85, 0.25*rows)]);
        th = tiledlayout(fh, rows, cols, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Frequency evaluation vectors for plotting fit lines
        x_eval_all  = log_freqs(interest_idx1);
        x_eval_low  = log_freqs_select(idx_low);
        x_eval_high = log_freqs_select(idx_high);

        for p = 1:n_plots
            ch = valid_plots(p);
            ax = nexttile(th);
            set(ax, 'Color', 'w');
            hold(ax, 'on');

            % 1. Plot raw log-log spectrum across the interest range
            plot(ax, log_freqs(interest_idx1), log_power_spectra(ch, interest_idx1), ...
                'k-', 'LineWidth', 1.3, 'DisplayName', 'Data');

            % 2. Global aperiodic linear fit
            y_glob = fit_params(ch).int_glob + fit_params(ch).slope_glob * x_eval_all;
            plot(ax, x_eval_all, y_glob, 'k--', 'LineWidth', 1.0, 'DisplayName', 'Global 1/f');

            % 3. Low-band fit
            y_low = fit_params(ch).int_low + fit_params(ch).slope_low * x_eval_low;
            plot(ax, x_eval_low, y_low, 'b-', 'LineWidth', 2.0, 'DisplayName', sprintf('Low (s=%.2f)', fit_params(ch).slope_low));

            % 4. High-band fit
            y_high = fit_params(ch).int_high + fit_params(ch).slope_high * x_eval_high;
            plot(ax, x_eval_high, y_high, 'r-', 'LineWidth', 2.0, 'DisplayName', sprintf('High (s=%.2f)', fit_params(ch).slope_high));

            % 5. Highlight censored band
            x_censor = log10(censor_range);
            xline(ax, x_censor(1), ':k', 'Alpha', 0.4, 'HandleVisibility', 'off');
            xline(ax, x_censor(2), ':k', 'Alpha', 0.4, 'HandleVisibility', 'off');

            grid(ax, 'on'); box(ax, 'on');
            xlim(ax, [log10(interest_range(1)), log10(interest_range(2))]);

            title(ax, sprintf('IC/Ch %d | \\Delta Slope: %.2f | NegInt: %.2f', ...
                ch, slope_discrepancy(ch), neg_integral(ch)), 'FontSize', 9, 'FontWeight', 'bold');

            if mod(p - 1, cols) == 0
                ylabel(ax, 'log_{10} Power', 'FontSize', 8);
            end
            if p > (rows - 1) * cols
                xlabel(ax, 'log_{10} Frequency (Hz)', 'FontSize', 8);
            end

            if p == 1
                legend(ax, 'Location', 'southwest', 'FontSize', 7);
            end
        end
    end
end

end