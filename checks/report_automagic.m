function metrics = report_automagic(EEG, num_bad_channels, varargin)
% CHECK_AUTOMAGICMETRICS Estimates quality measures from Pedroni et al. (2019)
%
% Inputs:
%   EEG              - EEGLAB data structure (epoched or continuous)
%   num_bad_channels - Integer, number of channels identified as bad/interpolated
%   varargin         - Optional name-value pairs for calculation thresholds (uV):
%                      'ThresholdOHA' (default 30 uV)
%                      'ThresholdStd' (default 15 uV)
%
% Outputs:
%   metrics          - Struct containing RBC, OHA, THV, CHV and Status.

% Parse inputs and set defaults based on Section 2.9.1
p = inputParser;
addRequired(p, 'EEG');
addRequired(p, 'num_bad_channels');
addParameter(p, 'ThresholdOHA', 30); % Voltage magnitude threshold (uV)
addParameter(p, 'ThresholdStd', 15); % Standard deviation threshold (uV)
parse(p, EEG, num_bad_channels, varargin{:});

thresh_x = p.Results.ThresholdOHA;
thresh_std = p.Results.ThresholdStd;

% Extract
data = EEG.data(1:128,:);
[n_chans, total_points] = size(data);
total_data_points = numel(data);

%% 1. Ratio of Bad Channels (RBC)
% Defined as the ratio of identified bad and consequently interpolated channels.
% Note: Ensure num_bad_channels includes interpolated ones.
metrics.RBC = num_bad_channels / n_chans;

%% 2. Ratio of Data with Overall High Amplitude (OHA)
% The ratio of data points (channel x time) where absolute voltage > x.
% Formula: Sum(|v| > x) / N [cite: 391]
high_amp_points = sum(abs(data(:)) > thresh_x);
metrics.OHA = high_amp_points / total_data_points;

%% 3. Ratio of Timepoints of High Variance (THV)
% The ratio of timepoints where the SD across all channels exceeds x.
% Formula: Sum(std_channels(v) > x) / T
std_across_channels = std(data, 0, 1);
metrics.THV = mean(std_across_channels > thresh_std);

%% 4. Ratio of Channels of High Variance (CHV)
% The ratio of channels where the SD across all timepoints exceeds x.
% Formula: Sum(std_time(v) > x) / C [cite: 403]
std_across_time = std(data, 0, 2);
metrics.CHV = mean(std_across_time > thresh_std);

%% 5. Classification (Adjusted to "Good", "OK", "Bad" from Section 2.9.2)
% The paper suggests determining objective exclusion criteria by categorising
% into three categories[cite: 418].

% Define Cutoffs (Percentages) - These are examples. The paper states users
% can modify default values and category cutoffs.
% "Bad" Thresholds (Strict rejection)
bad_RBC = 0.20;
bad_OHA = 0.20; % See Fig 2B
bad_THV = 0.25; % See Fig 2B
bad_CHV = 0.45; % See Fig 2B

% "OK" Thresholds (Lenient acceptance) - Adjust based on your population
ok_RBC = 0.15;
ok_OHA = 0.10;
ok_THV = 0.15;
ok_CHV = 0.15;

% Logic: Check if it exceeds "Bad" limits first
if metrics.RBC > bad_RBC || metrics.OHA > bad_OHA || ...
        metrics.THV > bad_THV || metrics.CHV > bad_CHV
    metrics.status = 'Bad';

    % If not Bad, check if it exceeds "OK" limits (meaning it is between Good and Bad)
elseif metrics.RBC > ok_RBC || metrics.OHA > ok_OHA || ...
        metrics.THV > ok_THV || metrics.CHV > ok_CHV
    metrics.status = 'OK';

else
    % If below all thresholds, it is Good
    metrics.status = 'Good';
end

end