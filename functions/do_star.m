function DATA = do_star(DATA, type_data, cfg)
% Reduce sparse artifacts
% https://www.sciencedirect.com/science/article/pii/S0165027016000066
%
% -> channel-specific artifats
% -> slow electrode drifts and pops
% -> does not fix EMG, EOG, ECG
%
% TODO:
% 1. If done on trials, correct the whole trial?
% 2.
%

fprintf('\n================================\n');
fprintf('STAR: Sparse time artifact removal (%s)\n', type_data);
fprintf('================================\n');

% Define params
if strcmpi(type_data, 'eeg')
    Niter   = 2;    % Number of STAR iterations
    Nneigh  = 8;    % Number of neighouring channels
    Texc    = 3;    % Threshold for excentricity, higher -> looser
    Ndeep   = 2;    % Maximum number of channels to fix at each sample
    Tpca    = 0.15; % Threshold for discarding weak PCs (percent of the max{PCs} of C of that neigh group of channels)
    Nsmooth = 64;   % Samples for smoothing applied on excentricity, too short -> too sensitive

elseif strcmpi(type_data, 'emg')
    Niter   = 2;
    Texc    = 5;
    Ndeep   = 3;
    Tpca    = 0.15;
    Nsmooth = 4*32;
end

% Define neighbours (only for EEG data)
if strcmpi(type_data, 'eeg')
    % EEG channels
    mask_channel = strcmp({DATA.chanlocs.type}, 'EEG');

    % Define neighbours
    fprintf('Extracting channel neighbours...\n');

    % Set number of neighbours (via distance)
    channel_neighbours = find_neighbours_knn(DATA.chanlocs(mask_channel), Nneigh);

    % % Variable number of neighbours (via triangulation)
    % channel_neighbours = find_neighbours(DATA.chanlocs(mask_channel));
    % num_channel = sum(mask_channel);
    % for i_channel = 1:num_channel
    %     mask = channel_neighbours(i_channel, :) == i_channel;
    %     channel_neighbours(i_channel, mask) = 0;
    % end

elseif strcmpi(type_data, 'emg')
    error('Not supported.');
    % mask_channel = strcmp({DATA.chanlocs.type}, 'EMG');
    % channel_neighbours = [];
end

% Extract
x_old = double(DATA.data(mask_channel, :))';
x_clean = x_old;
assert(ismatrix(x_old));

% Run
fprintf('Applying %d iterations.\n', Niter);
for i_iter = 1:Niter
    fprintf('Iteration: %d\n', i_iter);
    [x_clean, w, ww] = nt_star(x_clean, Texc, channel_neighbours, Ndeep, Tpca, Nsmooth);
end

% % Check
% EEGNEW = DATA;
% EEGNEW.data(mask_channel,:) = x_clean';
% vis_artifacts(EEGNEW, DATA);

% Store
DATA.data(mask_channel, :) = x_clean';

% Plot
fh = visualise_sparse_repairs(x_old, x_clean, ww, DATA.chanlocs(mask_channel), DATA.srate, cfg.figure.visible);
save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_star'], [30 15]);

end


function roiNeighbours = find_neighbours_knn(roilocs, K)
% Defines spatial neighbourhoods for TFCE using a K-Nearest Neighbours approach.
%
% Inputs:
%   roilocs : Structure containing .X, .Y, .Z coordinates
%   K       : The exact number of neighbours to assign to each location

% Extract coordinates
x = [roilocs(:).X]';
y = [roilocs(:).Y]';
z = [roilocs(:).Z]';
regionCentres = [x, y, z];

numRegions = size(regionCentres, 1);

% Calculate the complete 3D distance matrix instantly without a loop
distMatrix = pdist2(regionCentres, regionCentres);

% Preallocate the exact size required for the TFCE matrix
roiNeighbours = zeros(numRegions, K);

for i = 1:numRegions
    % Sort distances for the i-th region
    [~, indxSorted] = sort(distMatrix(i, :));

    % Extract the closest K channels.
    % We index from 2 to (K+1) because index 1 is a distance of 0 (the channel itself).
    neighbourIndices = indxSorted(2:(K + 1));

    % Store in the padded matrix
    roiNeighbours(i, :) = neighbourIndices;
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
percent_repaired = 1- mean(ww, 1);

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