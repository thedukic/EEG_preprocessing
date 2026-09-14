function [deviant_indices, fh] = check_channelcov(CovMatrices, subjects, options)
% CHECK_CHANNELCOV Identifies technical electrode defects (bad/disconnected/bridged channels)
% using Spatial Profile Laplacian error across BioSemi 128 channels.
%
% Inputs:
%   CovMatrices         : [N_chans x N_chans x N_sub] array (Covariance matrices recommended)
%   subjects            : Cell array of subject IDs (length M)
%   options             : (Optional) Struct with fields:
%       .z_threshold    : Outlier Z-score threshold (Default: 20)
%       .num_plots      : Maximum number of worst subjects to plot (Default: 5)
%
% Outputs:
%   deviant_indices     : Indices of flagged outlier subjects
%   fh                  : Figure handle

% -------------------------------------------------------------------------
% 1. Setup and Channel Location Loading
% -------------------------------------------------------------------------
if nargin < 3 || isempty(options), options = struct(); end
if ~isfield(options, 'z_threshold'), options.z_threshold = 20; end
if ~isfield(options, 'num_plots'),   options.num_plots = 5;     end

[n_channels, ~, n_subjects] = size(CovMatrices);
assert(n_subjects == numel(subjects), ...
    'Number of matrices (%d) must match number of subjects (%d).', n_subjects, numel(subjects));

% Load true BioSemi 128 channel labels and coordinates
mat_data = load('biosemi128_eeglab.mat', 'chanlocs');
chanlocs = mat_data.chanlocs(1:n_channels);
chan_labels = {chanlocs.labels};

% -------------------------------------------------------------------------
% 2. Build 3D Spatial Neighbourhood Matrix (Topological)
% -------------------------------------------------------------------------
% find_neighbours requires a struct with .X, .Y, .Z fields, which chanlocs has.
% It returns a zero-padded matrix where each row 'c' contains the neighbours
% of channel 'c', including channel 'c' itself.
raw_neighbours = find_neighbours(chanlocs);

% -------------------------------------------------------------------------
% 3. Compute Spatial Profile Laplacian Error (Profile Inconsistency)
% -------------------------------------------------------------------------
% chan_lap_err: [n_subjects x n_channels]
% Measures how much each channel's covariance profile violates local scalp smoothness
chan_lap_err = zeros(n_subjects, n_channels);

for i = 1:n_subjects
    R = CovMatrices(:, :, i);
    for c = 1:n_channels
        % Extract neighbours: remove zero-padding and exclude the channel itself
        c_neigh = raw_neighbours(c, :);
        valid_neigh = c_neigh(c_neigh > 0 & c_neigh ~= c);

        neigh_mean_profile = mean(R(valid_neigh, :), 1);
        diff_profile = R(c, :) - neigh_mean_profile;
        chan_lap_err(i, c) = norm(diff_profile, 2);
    end
end

% -------------------------------------------------------------------------
% 4. Cohort-Standardised Outlier Scoring
% -------------------------------------------------------------------------
cohort_med = median(chan_lap_err, 1);
cohort_mad = median(abs(chan_lap_err - cohort_med), 1) * 1.4826;

% Standardised anomaly matrix
chan_z_matrix = (chan_lap_err - cohort_med) ./ (cohort_mad + eps);

% Subject-level peak channel anomaly
[subj_peak_z, worst_ch_idx] = max(chan_z_matrix, [], 2);

% Flag subjects crossing threshold
deviant_indices = find(subj_peak_z > options.z_threshold);
[~, sort_order] = sort(subj_peak_z(deviant_indices), 'descend');
deviant_indices = deviant_indices(sort_order); % Sorted worst to least worst

% Group template for visualisation
CovMatricesMean = trimmean(CovMatrices, 20, 'round', 3);

% -------------------------------------------------------------------------
% 5. Console Diagnostic Summary & Full Ranking
% -------------------------------------------------------------------------
fprintf('\n===================================================================\n');
fprintf('  Spatial Profile Laplacian Diagnostics (BioSemi 128)\n');
fprintf('===================================================================\n');
fprintf('  Total Subjects:           %d\n', n_subjects);
fprintf('  Outlier Threshold:        Cohort Z > %.2f\n', options.z_threshold);
fprintf('  Flagged Outlier Subjects: %d\n', length(deviant_indices));
fprintf('  -----------------------------------------------------------------\n');

% Show top 5 highest-scoring subjects across the cohort regardless of cutoff
[all_sorted_z, all_order] = sort(subj_peak_z, 'descend');
fprintf('  Top Ranked Anomaly Scores in Cohort:\n');
for r = 1:min(5, n_subjects)
    s_i = all_order(r);
    ch_i = worst_ch_idx(s_i);
    flag_str = '';
    if subj_peak_z(s_i) > options.z_threshold, flag_str = ' [FLAGGED]'; end
    fprintf('    %d) %-15s | Peak Z = %5.2f | Lead: %s%s\n', ...
        r, subjects{s_i}, all_sorted_z(r), chan_labels{ch_i}, flag_str);
end
fprintf('===================================================================\n\n');

% -------------------------------------------------------------------------
% 6. Visualisation
% -------------------------------------------------------------------------
% Only plot the top 'num_plots' highest scoring subjects
top_plot_indices = all_order(1:min(options.num_plots, n_subjects));
n_plots = length(top_plot_indices);

if n_plots == 0
    fh = [];
    return;
end

Ncol = 4;
Nrow = n_plots;
fig_height = min(1200, 300 + (Nrow * 240));
fh = figure('Color', 'w', 'Name', 'Spatial Laplacian Channel Diagnostics', ...
    'Position', [100, 100, 1250, fig_height]);

th = tiledlayout(1 + Nrow, Ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
title(th, sprintf('Top %d Worst Participants (Threshold Z > %.1f)', ...
    n_plots, options.z_threshold), 'FontWeight', 'bold', 'FontSize', 12);

max_val = max(abs(CovMatricesMean(:)));
if max_val == 0, max_val = 1; end
c_lim = [-max_val, max_val];

try
    div_cmap = brewermap(256, '*RdBu');
catch
    div_cmap = jet(256);
end

% Top Row, Panel 1: Group Template Matrix
ax1 = nexttile(1);
imagesc(CovMatricesMean, c_lim);
colormap(ax1, div_cmap);
axis square; colorbar;
title('Group Average Template', 'FontWeight', 'bold', 'FontSize', 9);
xlabel('Channels'); ylabel('Channels');

% Top Row, Panels 2-4: Cohort Score Histogram
nexttile([1 3]);
histogram(subj_peak_z, 25, 'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'k');
hold on;
xline(options.z_threshold, '--r', sprintf('Threshold (Z = %.1f)', options.z_threshold), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
xlabel('Peak Spatial Laplacian Z-Score', 'FontSize', 9);
ylabel('Subject Count', 'FontSize', 9);
title('Distribution of Peak Channel Anomaly Scores', 'FontWeight', 'bold', 'FontSize', 10);
grid on;

% Outlier Rows
for i = 1:n_plots
    s_idx = top_plot_indices(i);
    diff_mat = CovMatrices(:, :, s_idx) - CovMatricesMean;
    ch_z     = chan_z_matrix(s_idx, :);
    [sorted_ch_z, ch_order] = sort(ch_z, 'descend');

    % Flag text for the title
    if subj_peak_z(s_idx) > options.z_threshold
        status = 'FLAGGED';
    else
        status = 'OK';
    end

    % Col 1: Subject Matrix
    ax_mat = nexttile;
    imagesc(CovMatrices(:, :, s_idx), c_lim);
    colormap(ax_mat, div_cmap);
    axis square; colorbar;
    title(sprintf('%s [%s]\n(Peak Z = %.2f)', subjects{s_idx}, status, subj_peak_z(s_idx)), ...
        'Interpreter', 'none', 'FontWeight', 'bold', 'FontSize', 9);

    % Col 2: Difference Matrix
    ax_diff = nexttile;
    max_d = max(abs(diff_mat(:)));
    if max_d == 0, max_d = 1e-4; end
    imagesc(diff_mat, [-max_d, max_d]);
    colormap(ax_diff, div_cmap);
    axis square; colorbar;
    title('Residual (Subj - Mean)', 'FontWeight', 'bold', 'FontSize', 9);

    % Col 3: Direct Call to mytopoplot (Laplacian Error Topomap)
    ax_topo = nexttile;
    top_k = min(10, n_channels);
    mytopoplot(ch_z(:), ch_order <= top_k, sprintf('Worst: %s (Z=%.1f)', chan_labels{ch_order(1)}, sorted_ch_z(1)), ax_topo);

    % Col 4: Top 5 Worst Leads Bar Plot
    nexttile;
    bar(sorted_ch_z(1:top_k), 'FaceColor', [0.85 0.25 0.25], 'EdgeColor', 'k');
    set(gca, 'XTick', 1:top_k, 'XTickLabel', chan_labels(ch_order(1:top_k)), 'FontSize', 8);
    ylabel('Cohort Z-Score', 'FontSize', 8);
    title(sprintf('Top Bad Lead: %s', chan_labels{ch_order(1)}), 'FontWeight', 'bold', 'FontSize', 9);
    grid on;
    yline(options.z_threshold, '--k', 'LineWidth', 1.0);
end

drawnow;

end