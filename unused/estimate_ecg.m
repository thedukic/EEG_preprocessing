function [opt_ecg, w_opt, chan_info, fh] = estimate_ecg(DATA, cfg)
% ESTIMATE_ECG
% Identifies cranial leads (LM, RM mandatory; LEL, REL optional) from an
% EEGLAB dataset, derives an optimal pseudo-ECG time-series using Generalized
% Eigendecomposition (GED), and visualises diagnostics.
%
% Inputs:
%   DATA      : EEGLAB dataset structure (2D continuous or 3D epoched)
%   cfg       : (Optional) Configuration structure:
%                 cfg.qrs_band    : QRS passband in Hz (Default: [10 25])
%                 cfg.gamma_reg   : Regularisation parameter (Default: 0.01)
%                 cfg.project_raw : If true, project broadband data; if false,
%                                   project QRS bandpass filtered data (Default: false)
%                 cfg.do_plot     : Render diagnostic figure (Default: true)
%                 cfg.plot_sec    : Duration in seconds to plot (Default: 10)
%
% Outputs:
%   opt_ecg   : Optimally combined pseudo-ECG signal [1 x n_pnts (x n_trials)]
%   w_opt     : [N x 1] Spatial filter weights
%   chan_info : Struct containing resolved channel indices, labels, and eigenvalues
%   fh        : Figure handle to the diagnostic plot

% -------------------------------------------------------------------------
% 1. Parameter Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'qrs_band'),    cfg.qrs_band    = [10, 25]; end % Hz
if ~isfield(cfg, 'gamma_reg'),   cfg.gamma_reg   = 0.01;     end
if ~isfield(cfg, 'project_raw'), cfg.project_raw = false;    end
if ~isfield(cfg, 'do_plot'),     cfg.do_plot     = true;     end
if ~isfield(cfg, 'plot_sec'),    cfg.plot_sec    = 10;       end % Seconds

fs = DATA.srate;
all_labels = {DATA.chanlocs.labels};

% -------------------------------------------------------------------------
% 2. Identify Cranial Channel Indices (LM/RM mandatory, LEL/REL optional)
% -------------------------------------------------------------------------
target_leads = {'LEL', 'REL', 'LM', 'RM'};
chan_idx = [];
matched_labels = {};

for k = 1:length(target_leads)
    lbl = target_leads{k};
    idx = find(strcmpi(all_labels, lbl), 1);

    if ~isempty(idx)
        chan_idx(end+1) = idx; %#ok<AGROW>
        matched_labels{end+1} = DATA.chanlocs(idx).labels; %#ok<AGROW>
    else
        if ismember(lbl, {'LM', 'RM'})
            error('estimate_ecg:ChannelNotFound', ...
                'Mandatory cranial lead %s not found in DATA.chanlocs.', lbl);
        else
            fprintf('Optional cranial lead %s not found; omitting from GED.\n', lbl);
        end
    end
end

n_cranial = length(chan_idx);
fprintf('Extracting pseudo-ECG using %d cranial leads: %s (Indices: %s)\n', ...
    n_cranial, strjoin(matched_labels, ', '), mat2str(chan_idx));

% -------------------------------------------------------------------------
% 3. Extract and Reshape Data Matrix
% -------------------------------------------------------------------------
is_epoched = (ndims(DATA.data) == 3);
if is_epoched
    [~, n_pnts, n_trials] = size(DATA.data);
    cranial_raw = double(DATA.data(chan_idx, :, :));
    cranial_2d  = reshape(cranial_raw, n_cranial, n_pnts * n_trials);
else
    n_trials   = 1;
    n_pnts     = size(DATA.data, 2);
    cranial_2d = double(DATA.data(chan_idx, :));
end

total_pnts = size(cranial_2d, 2);

% Demean each channel
cranial_2d = cranial_2d - mean(cranial_2d, 2);

% -------------------------------------------------------------------------
% 4. Create QRS Target Signal (10-25 Hz Bandpass)
% -------------------------------------------------------------------------
[b_qrs, a_qrs] = butter(2, cfg.qrs_band / (fs / 2), 'bandpass');
data_qrs_2d = filtfilt(b_qrs, a_qrs, cranial_2d')';

% -------------------------------------------------------------------------
% 5. Compute Covariance Matrices & Regularise
% -------------------------------------------------------------------------
% S: Target covariance (QRS power)
% R: Reference covariance (Broadband power)
S = (data_qrs_2d * data_qrs_2d') / (total_pnts - 1);
R = (cranial_2d  * cranial_2d')  / (total_pnts - 1);

gamma = cfg.gamma_reg * (trace(R) / size(R, 1));
R_reg = R + gamma * eye(size(R));

% -------------------------------------------------------------------------
% 6. Generalized Eigendecomposition (GED)
% -------------------------------------------------------------------------
[W, D] = eig(S, R_reg);
[evals, sort_idx] = sort(diag(D), 'descend');
W = W(:, sort_idx);

w_opt = W(:, 1);
w_opt = w_opt / norm(w_opt); % Normalise to unit vector

% -------------------------------------------------------------------------
% 7. Project Signal & Normalise Polarity
% -------------------------------------------------------------------------
if cfg.project_raw
    opt_ecg_2d = w_opt' * cranial_2d;
else
    opt_ecg_2d = w_opt' * data_qrs_2d;
end

% Ensure positive R-peak polarity (positive skewness)
if skewness(opt_ecg_2d) < 0
    opt_ecg_2d = -opt_ecg_2d;
    w_opt      = -w_opt;
end

% -------------------------------------------------------------------------
% 8. Format Output Dimensions & Diagnostic Struct
% -------------------------------------------------------------------------
if is_epoched
    opt_ecg = reshape(opt_ecg_2d, 1, n_pnts, n_trials);
else
    opt_ecg = opt_ecg_2d;
end

chan_info.chan_indices   = chan_idx;
chan_info.chan_labels    = matched_labels;
chan_info.eigenvalues    = evals;
chan_info.max_eigenvalue = evals(1);
chan_info.w_opt          = w_opt;

fprintf('GED pseudo-ECG extracted successfully (Max Eigenvalue SNR: %.2f).\n', evals(1));

% -------------------------------------------------------------------------
% 9. Diagnostic Visualisation
% -------------------------------------------------------------------------
fh = [];
if ~cfg.do_plot, return; end

n_plot_pnts = min(total_pnts, round(cfg.plot_sec * fs));
t_axis = (0:(n_plot_pnts - 1)) / fs;

fh = figure('Color', 'w', 'Position', [100, 100, 1100, 700], 'Name', 'Pseudo-ECG GED Diagnostics');

% --- Panel 1: Time-Series (Raw Cranial Leads vs GED ECG) ---
subplot(2, 2, [1, 3]);
hold on;

spacing = 3.5;
raw_snippet = cranial_2d(:, 1:n_plot_pnts);
for k = 1:n_cranial
    sig_norm = (raw_snippet(k, :) - mean(raw_snippet(k, :))) / std(raw_snippet(k, :));
    plot(t_axis, sig_norm + (n_cranial + 1 - k) * spacing, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
    text(-0.02 * cfg.plot_sec, (n_cranial + 1 - k) * spacing, matched_labels{k}, ...
        'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'FontSize', 9);
end

% Synthesised pseudo-ECG at bottom
ecg_snippet = opt_ecg_2d(1:n_plot_pnts);
ecg_norm = (ecg_snippet - mean(ecg_snippet)) / std(ecg_snippet);
plot(t_axis, ecg_norm, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.8);
text(-0.02 * cfg.plot_sec, 0, 'GED ECG', ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', 'Color', [0.85 0.1 0.1], 'FontSize', 9);

% Detect and mark candidate R-peaks
min_dist = round(0.50 * fs); % Min 500 ms RR-interval (max 120 BPM)
thresh_peak = 1.5;           % 1.5 SD threshold
[pks, locs] = findpeaks(ecg_norm, 'MinPeakDistance', min_dist, 'MinPeakHeight', thresh_peak);
if ~isempty(locs)
    plot(t_axis(locs), pks, 'v', 'MarkerEdgeColor', [0.1 0.5 0.1], ...
        'MarkerFaceColor', [0.2 0.8 0.2], 'MarkerSize', 6);
end

xlim([0, t_axis(end)]);
ylim([-3, (n_cranial + 1) * spacing + 1]);
xlabel('Time (s)', 'FontSize', 10);
ylabel('Normalised Amplitude (z-scored & stacked)', 'FontSize', 10);
title(sprintf('Cranial Inputs vs. Synthesised Pseudo-ECG (First %.1fs)', t_axis(end)), ...
    'FontWeight', 'bold', 'FontSize', 11);
grid on;
set(gca, 'YTick', []);

% --- Panel 2: Spatial Filter Weights ---
subplot(2, 2, 2);
bar(w_opt, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'k');
set(gca, 'XTick', 1:n_cranial, 'XTickLabel', matched_labels, 'FontSize', 9);
ylabel('Filter Weight (a.u.)', 'FontSize', 10);
xlabel('Cranial Lead', 'FontSize', 10);
title('Spatial Filter Weights (w_{opt})', 'FontWeight', 'bold', 'FontSize', 11);
grid on;
ylim([-1.1, 1.1]);

for k = 1:n_cranial
    y_pos = w_opt(k) + sign(w_opt(k)) * 0.12;
    text(k, y_pos, sprintf('%.2f', w_opt(k)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

% --- Panel 3: GED Eigenspectrum ---
subplot(2, 2, 4);
stem(1:n_cranial, evals, 'filled', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'MarkerSize', 7);
comp_labels = cell(1, n_cranial);
for k = 1:n_cranial
    comp_labels{k} = sprintf('Comp %d', k);
end
set(gca, 'XTick', 1:n_cranial, 'XTickLabel', comp_labels, 'FontSize', 9);
ylabel('Eigenvalue \lambda (QRS / Broadband SNR)', 'FontSize', 10);

if n_cranial > 1
    title(sprintf('GED Eigenspectrum (\\lambda_1 / \\lambda_2 Ratio: %.1f)', evals(1) / max(evals(2), eps)), ...
        'FontWeight', 'bold', 'FontSize', 11);
else
    title(sprintf('GED Eigenspectrum (\\lambda_1: %.1f)', evals(1)), ...
        'FontWeight', 'bold', 'FontSize', 11);
end
grid on;
xlim([0.5, n_cranial + 0.5]);

drawnow;

% Save figure if pipeline paths exist
if isfield(DATA, 'ALSUTRECHT') && isfield(DATA.ALSUTRECHT, 'subject')
    if isfield(DATA.ALSUTRECHT.subject, 'figures')
        save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_ecg_ged_diagnostic'], [20 15]);
    elseif isfield(DATA.ALSUTRECHT, 'figures') && isfield(DATA.ALSUTRECHT.figures, 'output_dir')
        save_figure(fh, DATA.ALSUTRECHT.figures.output_dir, [DATA.ALSUTRECHT.subject.id '_ecg_ged_diagnostic'], [20 15]);
    end
end

end