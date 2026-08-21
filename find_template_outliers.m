function [valid_idx, bad_mask, stats] = find_template_outliers(epoch_data, time_vec, cfg)
% FIND_TEMPLATE_OUTLIERS Rejects outlier event epochs (Blinks, Saccades, R-peaks)
%
% Inputs:
%   epoch_data : Matrix of size [n_samples x n_epochs] (single channel or component)
%   time_vec   : Vector of time points relative to event [1 x n_samples] (in seconds)
%   cfg        : Configuration structure with fields:
%                  .corr_thresh - Minimum correlation to median template (default: 0.70)
%                  .mad_thresh  - Max robust Z-score for amplitude/variance (default: 3.5)
%                  .peak_win    - [t_min t_max] search window for peak in s (default: [-0.08 0.08], or [] to skip)
%                  .max_abs_amp - Hard maximum amplitude cutoff in uV (default: Inf)
%                  .plot_diag   - true/false to generate diagnostic plot (default: false)
%
% Outputs:
%   valid_idx  : Vector of retained epoch indices
%   bad_mask   : Logical vector [1 x n_epochs] where true indicates an outlier
%   stats      : Structure containing individual metric values per epoch

if nargin < 3; cfg = struct(); end
if ~isfield(cfg, 'corr_thresh'); cfg.corr_thresh = 0.70; end
if ~isfield(cfg, 'mad_thresh');  cfg.mad_thresh  = 3.5;  end
if ~isfield(cfg, 'peak_win');    cfg.peak_win    = [-0.08 0.08]; end % Widened for blinks
if ~isfield(cfg, 'max_abs_amp'); cfg.max_abs_amp = Inf; end
if ~isfield(cfg, 'plot_diag');   cfg.plot_diag   = false; end

[n_samples, n_epochs] = size(epoch_data);
time_vec = time_vec(:)';
assert(length(time_vec) == n_samples, 'Length of time_vec must match number of rows in epoch_data.');

% -------------------------------------------------------------------------
% 1. Robust Median Template & Morphological Correlation
% -------------------------------------------------------------------------
if isempty(cfg.corr_thresh)
    bad_corr = false(1, n_epochs);
    r_vals  = nan(1, n_epochs);
else
    median_template = median(epoch_data, 2);

    % Zero-mean correlation for each epoch against median template
    r_vals = zeros(1, n_epochs);
    template_norm = median_template - mean(median_template);
    norm_t = norm(template_norm);

    for e = 1:n_epochs
        ep = epoch_data(:, e) - mean(epoch_data(:, e));
        denom = norm_t * norm(ep);
        if denom > 0
            r_vals(e) = dot(template_norm, ep) / denom;
        else
            r_vals(e) = 0;
        end
    end
    bad_corr = r_vals < cfg.corr_thresh;
end

% -------------------------------------------------------------------------
% 2. Peak Centring Check (Polarity-aware & Drift-robust)
% -------------------------------------------------------------------------
if isempty(cfg.peak_win)
    bad_latency = false(1, n_epochs);
    peak_times  = nan(1, n_epochs);
else
    % Identify dominant polarity of the median template near t = 0
    [~, idx_zero] = min(abs(time_vec));
    template_polarity = sign(median_template(idx_zero));
    if template_polarity == 0
        template_polarity = 1;
    end

    % Local search zone around t = 0 (prevents epoch edge drift from hijacking peak detection)
    search_margin = max(0.20, max(abs(cfg.peak_win)) * 1.5);
    search_mask   = (time_vec >= -search_margin) & (time_vec <= search_margin);
    search_idx    = find(search_mask);

    % Find maximum aligned with template polarity within the central zone
    signed_epochs = template_polarity * epoch_data(search_idx, :);
    [~, local_max_idx] = max(signed_epochs, [], 1);
    peak_idx = search_idx(local_max_idx);
    peak_times = time_vec(peak_idx);

    bad_latency = (peak_times < cfg.peak_win(1)) | (peak_times > cfg.peak_win(2));
end

% -------------------------------------------------------------------------
% 3. Amplitude & Variance Outlier Detection via MAD
% -------------------------------------------------------------------------
p2p_vals = max(epoch_data, [], 1) - min(epoch_data, [], 1);
std_vals = std(epoch_data, 0, 1);

% Robust Z-scores via Median Absolute Deviation
p2p_mad_z = abs(p2p_vals - median(p2p_vals)) / (1.4826 * mad(p2p_vals, 1));
std_mad_z = abs(std_vals - median(std_vals)) / (1.4826 * mad(std_vals, 1));

bad_amp = (p2p_mad_z > cfg.mad_thresh) | (std_mad_z > cfg.mad_thresh) | ...
    (max(abs(epoch_data), [], 1) > cfg.max_abs_amp);

% -------------------------------------------------------------------------
% 4. Combine Outlier Rejection Criteria
% -------------------------------------------------------------------------
bad_mask  = bad_corr | bad_latency | bad_amp;
valid_idx = find(~bad_mask);

% Output diagnostic metrics
stats.correlation   = r_vals;
stats.peak_times    = peak_times;
stats.p2p_amplitude = p2p_vals;
stats.p2p_mad_z     = p2p_mad_z;
stats.std_mad_z     = std_mad_z;
stats.bad_corr      = bad_corr;
stats.bad_latency   = bad_latency;
stats.bad_amp       = bad_amp;

fprintf('[Epoch Cleaning] Evaluated %d epochs: %d retained, %d rejected (%.1f%%)\n', ...
    n_epochs, length(valid_idx), sum(bad_mask), (sum(bad_mask)/n_epochs)*100);
fprintf('  - Rejected for Low Correlation (< %.2f) : %d\n', cfg.corr_thresh, sum(bad_corr));
if ~isempty(cfg.peak_win)
    fprintf('  - Rejected for Peak Misalignment       : %d\n', sum(bad_latency));
else
    fprintf('  - Peak Misalignment Check              : Skipped\n');
end
fprintf('  - Rejected for Amplitude / MAD Outlier : %d\n', sum(bad_amp));

% -------------------------------------------------------------------------
% 5. Optional Diagnostic Visualisation
% -------------------------------------------------------------------------
if cfg.plot_diag
    figure('Color', 'w', 'Name', 'Epoch Cleaning Diagnostics');

    subplot(1, 2, 1);
    if any(bad_mask)
        plot(time_vec, epoch_data(:, bad_mask), 'Color', [0.85 0.3 0.3 0.4], 'LineWidth', 0.8); hold on;
    end
    plot(time_vec, epoch_data(:, ~bad_mask), 'Color', [0.2 0.7 0.2 0.3], 'LineWidth', 0.8); hold on;
    plot(time_vec, median_template, 'k', 'LineWidth', 2.0);
    xlim([time_vec(1), time_vec(end)]);
    xlabel('Time (s)'); ylabel('Amplitude (\muV)');
    title(sprintf('Epochs (Kept: %d, Discarded: %d)', length(valid_idx), sum(bad_mask)));
    grid on;

    subplot(1, 2, 2);
    scatter(r_vals(~bad_mask), p2p_mad_z(~bad_mask), 25, [0 0.5 0], 'filled', 'MarkerFaceAlpha', 0.6); hold on;
    if any(bad_mask)
        scatter(r_vals(bad_mask), p2p_mad_z(bad_mask), 25, [0.8 0 0], 'filled', 'MarkerFaceAlpha', 0.6);
    end
    xline(cfg.corr_thresh, '--r', 'Correlation Cutoff');
    yline(cfg.mad_thresh, '--r', 'MAD Cutoff');
    xlabel('Template Correlation (r)'); ylabel('P2P Amplitude MAD Z-Score');
    title('Decision Space');
    grid on;
end
end