function plot_toposnapshots(data, srate, window_ms)
% PLOT_DATA_SNAPSHOTS Plots snapshots of channel activity averaged over time windows.
%
% This function takes a 2D data matrix (channels x time), a sampling rate,
% and a window length in milliseconds, and plots the average channel activity
% (the 'snapshot') for each defined time window.
%
% Syntax:
%   plot_data_snapshots(data, srate, window_ms)
%
% Inputs:
%   data        - 2D matrix (M x N), where M is the number of channels
%                 and N is the number of time points.
%   srate       - The sampling rate of the data in Hz (samples per second).
%   window_ms   - The desired window length for averaging in milliseconds (ms).
%
% Example Usage (run this section for a demonstration):
%   srate = 1000; % 1000 Hz (1 sample per ms)
%   num_channels = 32;
%   total_time_s = 5;
%   time_points = srate * total_time_s;
%
%   % Create synthetic data with some evolving spatial pattern
%   t = (0:time_points-1) / srate;
%   time_mod = repmat(cos(2*pi*1*t), num_channels, 1); % 1 Hz oscillation
%   channel_mod = repmat((1:num_channels)', 1, time_points); % Channel gradient
%   noise = 0.5 * randn(num_channels, time_points);
%   data_sim = (channel_mod .* time_mod) + noise;
%
%   window_ms = 250; % Average over 250 ms windows
%   plot_data_snapshots(data_sim, srate, window_ms);

% --- Input Validation ---
if nargin < 3
    error('Not enough input arguments. Requires data, srate, and window_ms.');
end
if isempty(data)
    disp('Data matrix is empty. Nothing to plot.');
    return;
end

% --- Parameter Calculation ---
[num_channels, num_timepoints] = size(data);

% Convert window_ms to samples
window_samples = round(window_ms / 1000 * srate);

if window_samples < 1
    error('Window length is too short. It must correspond to at least one sample.');
end

% Determine the number of non-overlapping windows
num_windows = floor(num_timepoints / window_samples);

if num_windows < 1
    disp('Not enough data to form a single window of the specified size.');
    return;
end

% --- Windowing and Averaging ---

% Pre-allocate matrix to store the averaged snapshots (channels x windows)
averaged_snapshots = zeros(num_channels, num_windows);

for w = 1:num_windows
    % Calculate start and end indices for the current window
    start_idx = (w - 1) * window_samples + 1;
    end_idx = w * window_samples;

    % Extract the data segment
    segment = data(:, start_idx:end_idx);

    % Calculate the average across the time dimension (2)
    averaged_snapshots(:, w) = mean(segment, 2);
end

% --- Visualization ---
figure('Name', sprintf('Averaged Snapshots (Window: %d ms)', window_ms), 'Color', 'w');

% % Determine the maximum number of snapshots to display (e.g., up to 12)
% max_plots = min(num_windows, 12);
%
% % Calculate subplot layout
% rows = ceil(sqrt(max_plots));
% cols = ceil(max_plots / rows);

% Calculate subplot layout
max_plots = num_windows;
rows = 2;
cols = ceil(max_plots / rows);

% Plot each snapshot
tiledlayout(rows, cols, "TileSpacing", "compact", "Padding", "compact");
for w = 1:max_plots
    % Calculate time boundaries for the title
    time_start_ms = (w - 1) * window_ms;
    time_end_ms = w * window_ms;

    % Plot
    mytopoplot(averaged_snapshots(:, w), [], [], nexttile); axis tight; drawnow;
    title(sprintf('%d - %d ms', time_start_ms, time_end_ms), 'FontSize', 10);
end

sgtitle(sprintf('Spatial Snapshots Averaged over %d ms Windows (Total windows: %d)', window_ms, num_windows), 'FontWeight', 'bold', 'FontSize', 14);

end