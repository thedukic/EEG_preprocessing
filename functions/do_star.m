function DATA = do_star(DATA, type_data, cfg)
% DO_STAR Sparse Time Artefact Removal (STAR) with Dual-Pass Mask Transfer
%
% Syntax:
%   DATA = do_star(DATA, type_data, cfg)
%
% Description:
%   Implements an improved version of the Sparse Time Artefact Removal (STAR)
%   algorithm (de Cheveigné, 2016, J Neurosci Methods) to suppress channel-specific
%   artefacts (electrode pops, contact glitches, and localised myogenic bursts)
%   without discarding channels, dropping trials, or reducing data rank.
%   https://www.sciencedirect.com/science/article/pii/S0165027016000066
%
% Key Improvements in this Version:
%
%   1. Dual-Pass Mask Transfer (0.5 Hz Detection -> 0.1 Hz Target Repair):
%      Artefact detection is performed on a temporary 0.5 Hz high-pass filtered
%      stream (DATA_TRAIN). This prevents sub-delta fluctuations and slow skin
%      potentials from corrupting the excentricity metric. The resulting sample-
%      and channel-specific intervention mask (ww) is then transferred directly
%      to guide repairs on the target 0.1 Hz data.
%
%   2. Boundary Offset Stitching (Elimination of Edge Jumps):
%      Resolves edge step discontinuities at repair transitions. Because residual
%      low-frequency power (0.1 to 0.5 Hz) creates local baseline offsets between
%      the target electrode and its neighbours' spatial projection, hard sample
%      replacements can inject artificial high-frequency pops. This implementation
%      measures boundary mismatches at the entry and exit of each contiguous bad
%      segment and applies a linear baseline-matching ramp across the patch for
%      seamless docking.
%
%   3. Spatially Constrained Projections (KNN Geometry):
%      Replaces unconstrained global projection (all N-1 channels) with a K-nearest
%      neighbourhood (e.g. 8 to 10 electrodes for 128-channel montages). This
%      prevents spatial overfitting and ensures projections do not pull in distant,
%      uncorrelated myogenic activity from opposite cranial margins.
%
%   4. Calibrated Subspace Regularisation (threshold_pca):
%      Tunes PCA eigenvalue thresholding (typically threshold_pca = 0.01) to eliminate
%      collinear or degenerate dimensions while preserving the multi-dimensional
%      topographical structure of local sensor clusters.
%
%   5. Vectorised Performance & Precomputed Projections:
%      Replaces iterative filter calls with vectorised moving averages (movmean)
%      and computes projection matrices once per iteration from the clean covariance
%      matrix c0, avoiding redundant matrix inversions in deep repair loops.
%
% Inputs:
%   DATA       - EEGLAB data structure containing .data, .chanlocs, and .srate
%   type_data  - String specifying data modality ('eeg' supported; 'emg' placeholder)
%   cfg        - Configuration structure containing filtering and figure display fields
%
% Outputs:
%   DATA       - Processed EEGLAB structure with repaired continuous EEG data
%
% TODO:
% 1. If done on trials, correct the whole trial
% 2.

fprintf('\n================================\n');
fprintf('STAR: Sparse time artifact removal (%s)\n', type_data);
fprintf('================================\n');

% Define params
if strcmpi(type_data, 'eeg')
    num_iter      = 2;    % Number of STAR iterations
    num_neigh     = 10;   % Number of neighouring channels
    num_deep      = 3;    % Maximum number of channels to fix at each sample
    num_smooth    = 32;   % Samples for smoothing applied on excentricity, too short -> too sensitive
    threshold_exc = 3;    % Threshold for excentricity, higher -> looser
    threshold_pca = 0.01; % Threshold for discarding weak PCs (percent of the max{PCs} of C of that neigh group of channels)

elseif strcmpi(type_data, 'emg')
    error('Not supported.');

    num_iter      = 2;
    num_deep      = 3;
    num_smooth    = 32;
    threshold_exc = 5;
    threshold_pca = 0.01;
end

% Demean
DATA.data = DATA.data - mean(DATA.data, 2);

% Filter for training
% cut_off = 0.5;
DATA_TRAIN = do_filtering_method(DATA, cfg.star.cut_off, cfg);

% Define neighbours (only for EEG data)
if strcmpi(type_data, 'eeg')
    % EEG channels
    mask_channel = strcmp({DATA.chanlocs.type}, 'EEG');

    % Define neighbours
    fprintf('Extracting channel neighbours...\n');

    % Set number of neighbours (via distance)
    channel_neighbours = find_neighbours_knn(DATA.chanlocs(mask_channel), num_neigh);

    % % Variable number of neighbours (via triangulation)
    % channel_neighbours = find_neighbours(DATA.chanlocs(mask_channel));
    % num_channel = sum(mask_channel);
    % for i_channel = 1:num_channel
    %     mask = channel_neighbours(i_channel, :) == i_channel;
    %     channel_neighbours(i_channel, mask) = 0;
    % end

elseif strcmpi(type_data, 'emg')
    % mask_channel = strcmp({DATA.chanlocs.type}, 'EMG');
    % channel_neighbours = [];
end

% Extract
x_old = double(DATA.data(mask_channel, :))';
x_train_old = double(DATA_TRAIN.data(mask_channel, :))';

x_new = x_old;
x_train_new = x_train_old;

% Must be continuous
assert(ismatrix(x_new));
assert(ismatrix(x_train_new));

% Run
fprintf('Applying %d iterations.\n', num_iter);
for i_iter = 1:num_iter
    fprintf('Iteration: %d\n', i_iter);

    % 1. Detect artifacts on the 0.5 Hz filtered training data
    [x_train_new, w, ww] = nt_star_improved(x_train_new, threshold_exc, channel_neighbours, num_deep, threshold_pca, num_smooth);

    % 2. Extract target covariance restricted to clean samples (using w directly)
    c0_target = nt_cov(x_new, [], w);

    % % 3. Apply repairs to the target 0.1 Hz data
    % for i_channel = 1:size(x_new, 2)
    %     bad_samples = find(ww(:, i_channel) == 0);
    %     if isempty(bad_samples), continue; end
    %
    %     if ~isempty(channel_neighbours)
    %         oC = channel_neighbours(i_channel, :);
    %     else
    %         oC = setdiff(1:size(x_new, 2), i_channel);
    %     end
    %     oC(oC > size(x_new, 2) | oC <= 0) = [];
    %     if isempty(oC), continue; end
    %
    %     c_sub = c0_target(oC, oC);
    %     [topcs, eigenvalues] = nt_pcarot(c_sub);
    %     keep = (eigenvalues / max(eigenvalues)) > Tpca;
    %     topcs = topcs(:, keep);
    %
    %     b = (c0_target(i_channel, oC) * topcs) / (topcs' * c_sub * topcs);
    %     projWeights = topcs * b';
    %
    %     x_new(bad_samples, i_channel) = x_new(bad_samples, oC) * projWeights;
    % end

    % 3. Apply repairs to the target data with boundary offset matching
    for i_channel = 1:size(x_new, 2)
        bad_idx = find(ww(:, i_channel) == 0, 1);
        if isempty(bad_idx), continue; end

        if ~isempty(channel_neighbours)
            oC = channel_neighbours(i_channel, :);
        else
            oC = setdiff(1:size(x_new, 2), i_channel);
        end
        oC(oC > size(x_new, 2) | oC <= 0) = [];
        if isempty(oC), continue; end

        c_sub = c0_target(oC, oC);
        [topcs, eigenvalues] = nt_pcarot(c_sub);
        keep = (eigenvalues / max(eigenvalues)) > threshold_pca;
        topcs = topcs(:, keep);

        b = (c0_target(i_channel, oC) * topcs) / (topcs' * c_sub * topcs);
        projWeights = topcs * b';

        % Project entire channel candidate
        y_full = x_new(:, oC) * projWeights;

        % Identify contiguous segments of bad samples
        diff_bad = diff([0; ww(:, i_channel) == 0; 0]);
        seg_starts = find(diff_bad == 1);
        seg_ends   = find(diff_bad == -1) - 1;

        for i_seg = 1:numel(seg_starts)
            idx_start = seg_starts(i_seg);
            idx_end   = seg_ends(i_seg);
            L = idx_end - idx_start + 1;

            % Measure offsets at boundary edges
            if idx_start > 1
                d_start = x_new(idx_start - 1, i_channel) - y_full(idx_start);
            else
                d_start = 0;
            end

            if idx_end < size(x_new, 1)
                d_end = x_new(idx_end + 1, i_channel) - y_full(idx_end);
            else
                d_end = 0;
            end

            % Create linear baseline correction ramp
            offset_ramp = linspace(d_start, d_end, L)';

            % Splice seamlessly
            x_new(idx_start:idx_end, i_channel) = y_full(idx_start:idx_end) + offset_ramp;
        end
    end
end

% % Check
% DATA_NEW = DATA;
% DATA_NEW.data(mask_channel, :) = x_new';
%
% cut_off = 1;
% DATA_PLOT     = do_filtering_method(DATA, cut_off, cfg);
% DATA_NEW_PLOT = do_filtering_method(DATA_NEW, cut_off, cfg);
% vis_artifacts(DATA_NEW_PLOT, DATA_PLOT);

% Store
DATA.data(mask_channel, :) = x_new';

% Plot
fh = visualise_sparse_repairs(x_old, x_new, ww, DATA.chanlocs(mask_channel), DATA.srate, cfg.figure.visible);
save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_star'], [30 15]);

end


function neighbours = find_neighbours_knn(locs, K)
% Defines spatial neighbourhoods for TFCE using a K-Nearest Neighbours approach.
%
% Inputs:
%   locs : Structure containing .X, .Y, .Z coordinates
%   K    : The exact number of neighbours to assign to each location

% Extract coordinates
x = [locs(:).X]';
y = [locs(:).Y]';
z = [locs(:).Z]';
points = [x, y, z];
num_channel = size(points, 1);

% Calculate the complete 3D distance matrix instantly without a loop
distances = pdist2(points, points);

% Preallocate the exact size required for the TFCE matrix
neighbours = zeros(num_channel, K);

for i_channel = 1:num_channel
    % Sort distances for the i-th region
    [~, indx_sorted] = sort(distances(i_channel, :));

    % Extract the closest K channels.
    % We index from 2 to (K+1) because index 1 is a distance of 0 (the channel itself).
    indx_neighbour = indx_sorted(2:(K + 1));

    % Store in the padded matrix
    neighbours(i_channel, :) = indx_neighbour;
end

end

function fh = visualise_sparse_repairs(x_raw, x_clean, ww, chanlocs, srate, figure_visible)
% VISUALISE_SPARSE_REPAIRS Generates a diagnostic summary figure for nt_star.
%
% Inputs:
%   x_raw    : The original data matrix before nt_star (samples x channels)
%   x_clean  : The output data matrix from nt_star (samples x channels)
%   ww       : The intervention mask from nt_star (samples x channels)
%   chanlocs : EEGLAB chanlocs structure (optional; pass [] if not available)
%   srate    : Sampling rate in Hz for the time axis (optional)

if nargin < 5 || isempty(srate)
    srate = 1;
    time_axis = 1:size(x_raw, 1);
    time_label = 'Samples';
else
    time_axis = (0:size(x_raw, 1)-1) / srate;
    time_label = 'Time (s)';
end

num_channels = size(ww, 2);

% Calculate the percentage of samples repaired per channel
% repairs_per_chan = sum(ww == 0, 1);
% percent_repaired = (repairs_per_chan / size(ww, 1)) * 100;
percent_repaired = 1 - mean(ww, 1);

% Initialise the figure
fh = figure('Name', 'Sparse Artifact Repairs Summary', 'Color', 'w', 'Position', [100, 100, 1400, 800], 'Visible', figure_visible);
% tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');
% =========================================================================
% PANEL 1: Intervention Heatmap
% =========================================================================
ax1 = nexttile(tl);

% Transpose ww so time is on the X axis to match the butterfly plot
imagesc(ax1, time_axis, 1:num_channels, ww');

% Custom colour map: Red for 0 (fixed), Light grey for 1 (untouched)
colormap(ax1, [0.85 0.33 0.10; 0.95 0.95 0.95]);

xlabel(ax1, time_label);
ylabel(ax1, 'Channel Index');
title(ax1, 'Raster Map of Algorithmic Interventions');
set(ax1, 'YDir', 'normal'); % Set channel 1 at the bottom

% =========================================================================
% PANEL 2: Spatial Distribution of Repairs
% =========================================================================
ax2 = nexttile(tl);

if nargin >= 4 && ~isempty(chanlocs)
    % If channel locations are provided, project the repair density onto the scalp
    mytopoplot(percent_repaired, false(1, 128), 'Spatial Density of Repairs (%)', ax2);
    % title(ax2, 'Spatial Density of Repairs (%)');
    colorbar(ax2);
else
    % Fallback to a standard bar plot if locations are missing
    bar(ax2, percent_repaired, 'FaceColor', [0.15 0.55 0.82], 'EdgeColor', 'none');
    xlabel(ax2, 'Channel Index');
    ylabel(ax2, 'Data Repaired (%)');
    title(ax2, 'Proportion of Repaired Data per Channel');
    grid(ax2, 'on');
    ax2.GridAlpha = 0.1;
    ax2.Box = 'off';
    xlim(ax2, [0 num_channels+1]);
end

% % =========================================================================
% % PANEL 3: Global Butterfly Overlay
% % =========================================================================
% ax3 = nexttile(tl, [1 2]); % Span the entire bottom row
% hold(ax3, 'on');
%
% % Plot the raw data in faint red
% plot(ax3, time_axis, x_raw, 'Color', [0.85 0.33 0.10 0.4], 'LineWidth', 1);
%
% % Plot the cleaned data in solid blue on top
% plot(ax3, time_axis, x_clean, 'Color', [0.0 0.4 0.7 0.8], 'LineWidth', 1);
%
% xlabel(ax3, time_label);
% ylabel(ax3, 'Amplitude (\muV)');
% title(ax3, 'Global Butterfly Overlay (Red: Raw Wobbles, Blue: Cleaned Signal)');
% grid(ax3, 'on');
% ax3.GridAlpha = 0.1;
% ax3.Box = 'off';
% xlim(ax3, [min(time_axis) max(time_axis)]);
%
% % Set dynamic Y-limits based on the cleaned data to prevent massive raw spikes
% % from zooming the plot out too far
% y_limit = max(abs(x_clean(:))) * 1.5;
% if y_limit > 0
%     ylim(ax3, [-y_limit y_limit]);
% end

end

function [x, w, ww] = nt_star_improved(x, thresh, closest, depth, PCA_THRESH, NSMOOTH)
% [y, w, ww] = nt_star_improved(x, thresh, closest, depth, PCA_THRESH, NSMOOTH)
%
% Sensor noise suppression via Sparse Time Artifact Removal (STAR).
%
% y:          Denoised data
% w:          0 for parts needing fixing, 1 elsewhere (time*1)
% ww:         0 for parts needing fixing, 1 elsewhere (time*chans)
%
% x:          Data to denoise (time*chans or time*chans*trials)
% thresh:     Threshold for excentricity measure (default: 1)
% closest:    Indices of closest channels per row (default: all)
% depth:      Maximum number of channels to fix per sample (default: 1)
% PCA_THRESH: Threshold for discarding weak PCs (default: 1e-15)
% NSMOOTH:    Smoothing window in samples for excentricity (default: 10)

MINPROP = 0.3;  % Minimum proportion of artifact-free data at first iteration
NITER   = 3;    % Iterations to refine c0
VERBOSE = 1;

if nargin < 1, error('More inputs needed.'); end
if nargin < 2 || isempty(thresh), thresh = 1; end
if nargin < 3, closest = []; end
if nargin < 4 || isempty(depth), depth = 1; end
if nargin < 5 || isempty(PCA_THRESH), PCA_THRESH = 1e-15; end
if nargin < 6 || isempty(NSMOOTH), NSMOOTH = 10; end

if ~isempty(closest) && size(closest, 1) ~= size(x, 2)
    error('closest array should have as many rows as channels of x');
end
if any(closest < 0, 'all'), error('Invalid channel indices.'); end

if nargout == 0
    [y, w, ww] = nt_star_improved(x, thresh, closest, depth, PCA_THRESH, NSMOOTH);
    disp([mean(w(:)), mean(ww(:))]);
    figure(1); clf;
    subplot(3, 1, 1); plot(x); title('Original');
    subplot(3, 1, 2); plot(y); title('Cleaned');
    subplot(3, 1, 3); plot(w, '.-'); ylim([0 1.1]); title('Clean Mask w');
    clear x w ww;
    return;
end

[num_time, num_channel, ~] = size(x);
x = nt_unfold(x);

p0 = nt_wpwr(x);
mn = mean(x, 1);
x  = nt_demean(x);
nn = sqrt(mean(x.^2, 1));
nn(nn == 0) = 1; % Guard against division by zero on dead channels
x  = x ./ nn;
p00 = nt_wpwr(x);

% Flag and isolate flat / NaN channels
iNan  = find(all(isnan(x), 1));
iZero = find(all(x == 0, 1));
x(:, iNan)  = randn(size(x, 1), numel(iNan));
x(:, iZero) = randn(size(x, 1), numel(iZero));

x  = nt_demean(x);
c0 = nt_cov(x);

% Parse neighbour list for each channel
oChanList = cell(num_channel, 1);
for i_channel = 1:num_channel
    if ~isempty(closest)
        oC = closest(i_channel, :);
    else
        oC = setdiff(1:num_channel, i_channel);
    end
    oC(oC > num_channel | oC <= 0) = [];
    oChanList{i_channel} = oC;
end

% -------------------------------------------------------------------------
% Phase 1: Iterative clean mask and covariance estimation
% -------------------------------------------------------------------------
i_iter = NITER;
while i_iter > 0
    w = ones(size(x, 1), 1);

    for i_channel = 1:num_channel
        oC = oChanList{i_channel};
        if isempty(oC), continue; end

        c_sub = c0(oC, oC);
        [topcs, eigenvalues] = nt_pcarot(c_sub);
        keep = (eigenvalues / max(eigenvalues)) > PCA_THRESH;
        topcs = topcs(:, keep);

        % Projection vector: x(:, oC) * P_chan predicts x(:, iChan)
        b = (c0(i_channel, oC) * topcs) / (topcs' * c_sub * topcs);
        P_chan = topcs * b';

        y  = x(:, oC) * P_chan;
        dx = abs(y - x(:, i_channel)) + eps;

        % Excentricity relative to current clean samples
        d = dx / mean(dx(w > 0));
        if NSMOOTH > 0
            d = movmean(d, NSMOOTH, 1);
        end

        d = d / thresh;
        w = min(w, double(d < 1));
    end

    prop = mean(w);
    if VERBOSE > 0, disp(['Proportion artifact free: ', num2str(prop)]); end

    if i_iter == NITER && prop < MINPROP
        thresh = 1.1 * thresh;
        if VERBOSE > 0, disp(['Warning: nt_star increasing threshold to ', num2str(thresh)]); end
        w = ones(size(w));
    else
        i_iter = i_iter - 1;
    end

    x  = nt_demean(x, w);
    c0 = nt_cov(x, [], w);
end

% -------------------------------------------------------------------------
% Phase 2: Excentricity ranking and local channel repair
% -------------------------------------------------------------------------
std_clean = sqrt(mean(x(w > 0, :).^2, 1));
std_clean(std_clean == 0) = 1;
xx = abs(x) ./ std_clean;

if NSMOOTH > 0
    xx = movmean(xx, NSMOOTH, 1);
end

% Precompute projection vectors for all channels (computed once from stable c0)
projWeights = cell(num_channel, 1);
for i_channel = 1:num_channel
    oC = oChanList{i_channel};
    if isempty(oC), continue; end

    c_sub = c0(oC, oC);
    [topcs, eigenvalues] = nt_pcarot(c_sub);
    keep = (eigenvalues / max(eigenvalues)) > PCA_THRESH;
    topcs = topcs(:, keep);

    b = (c0(i_channel, oC) * topcs) / (topcs' * c_sub * topcs);
    projWeights{i_channel} = topcs * b';
end

% Identify worst channels at contaminated points
depth = min(depth, num_channel - 1);
[~, rank] = sort(xx, 2, 'descend');
rank(w > 0, :) = 0; % Do not touch samples confirmed artifact-free

ww = ones(size(x));
for i_depth = 1:depth
    iFixed = num_channel;
    nFixed = 0;

    for i_channel = 1:num_channel
        bad_samples = find(rank(:, i_depth) == i_channel);

        if i_depth > 1 && ~isempty(bad_samples)
            bad_samples(xx(bad_samples, i_channel) < thresh) = [];
        end

        if isempty(bad_samples)
            iFixed = iFixed - 1;
            continue;
        end

        nFixed = nFixed + numel(bad_samples);
        ww(bad_samples, i_channel) = 0;

        oC = oChanList{i_channel};
        if ~isempty(oC) && ~isempty(projWeights{i_channel})
            % Fix bad samples using precomputed spatial projection
            x(bad_samples, i_channel) = x(bad_samples, oC) * projWeights{i_channel};
        end
    end

    if VERBOSE > 0
        disp(['Depth: ', num2str(i_depth), ...
            ', channels active: ', num2str(iFixed), ...
            ', samples repaired: ', num2str(nFixed), ...
            ', power ratio: ', num2str(nt_wpwr(x) / p00)]);
    end
end

% Restore channel scaling, offsets, and original structure
x = nt_demean(x);
x = (x .* nn) + mn;

x  = nt_fold(x, num_time);
w  = nt_fold(w, num_time);
ww = nt_fold(ww, num_time);

x(:, iNan, :)  = nan;
x(:, iZero, :) = 0;

if VERBOSE > 0, disp(['Final power ratio: ', num2str(nt_wpwr(x) / p0)]); end
end