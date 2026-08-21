function [eyeBlinksMask, eyeBlinksEpochs, BlinkMaxLatency, eyeBlinkData, brainData, threshold] = detect_veog(DATA, winBlink, trIQR, optVisible)

% -------------------------------------------------------------------------
% 1. Extract and Filter VEOG
% -------------------------------------------------------------------------
chaneog = strcmp({DATA.chanlocs.labels}, 'VEOG');
assert(sum(chaneog) == 1, 'Could not uniquely identify VEOG channel.');
eyeBlinkData = DATA.data(chaneog, :);

% Temporarily filter VEOG, 0.2-10 Hz
[bl, al] = butter(2, 10/(DATA.srate/2), 'low');
[bh, ah] = butter(2, 0.2/(DATA.srate/2), 'high');
eyeBlinkData = do_filteringcore(bl, al, eyeBlinkData, DATA.event, DATA.srate);
eyeBlinkData = do_filteringcore(bh, ah, eyeBlinkData, DATA.event, DATA.srate);

brainData = []; % Placeholder to match function outputs

% -------------------------------------------------------------------------
% 2. Thresholding (RELAX Approach)
% -------------------------------------------------------------------------
EOGIQR = iqr(eyeBlinkData);
EOG75P = prctile(eyeBlinkData, 75);

threshold = EOG75P + (trIQR * EOGIQR);
threshold = max(threshold, 100); % Hard minimum of 100uV

fprintf('Eye blinks (L = +-%d ms) were detected using a threshold of 75PRC + %dIQR = %1.0fuV.\n', winBlink, trIQR, threshold);

BlinkIndexMetric = double(eyeBlinkData > threshold);

% -------------------------------------------------------------------------
% 3. Detect Blinks (Duration > 50ms)
% -------------------------------------------------------------------------
% Initialize output variable to prevent crash if no blinks are found
BlinkMaxLatency = [];

% Pad with zeros to guarantee start/end pairs line up perfectly
paddedMetric = [0, BlinkIndexMetric, 0];
crossings = diff(paddedMetric);

ix_blinkstart = find(crossings == 1);
ix_blinkend   = find(crossings == -1) - 1; % -1 corrects for the padding shift

mspersmpl = 1000 / DATA.srate;
min_duration_smpl = round(50 / mspersmpl);

if ~isempty(ix_blinkstart)
    BlinkLengths = ix_blinkend - ix_blinkstart;
    valid_runs = find(BlinkLengths > min_duration_smpl);

    if ~isempty(valid_runs)
        BlinkMaxLatency = zeros(1, length(valid_runs));
        for x = 1:length(valid_runs)
            start_idx = ix_blinkstart(valid_runs(x));
            end_idx   = ix_blinkend(valid_runs(x));

            % Find local maximum within the thresholded window
            [~, max_offset] = max(eyeBlinkData(start_idx:end_idx));
            BlinkMaxLatency(x) = start_idx + max_offset - 1;
        end
    end
end

% -------------------------------------------------------------------------
% 4. Epoching and Boundary Checks
% -------------------------------------------------------------------------
winBlinksmpl = round(winBlink / mspersmpl);

if isempty(BlinkMaxLatency)
    eyeBlinksEpochs = [];
    eyeBlinksMask = false(size(eyeBlinkData));
    fprintf('\nNo valid eye blinks detected.\n');
    return; % Exit early, nothing more to do
else
    fprintf('\nNumber of blinks: %d\n', length(BlinkMaxLatency));
end

% Create epochs
eyeBlinksEpochs = [BlinkMaxLatency' - winBlinksmpl, BlinkMaxLatency' + winBlinksmpl];

% Vectorized boundary check: keep only epochs that fit entirely in the data
valid_epochs = (eyeBlinksEpochs(:,1) >= 1) & (eyeBlinksEpochs(:,2) <= length(eyeBlinkData));

if sum(~valid_epochs) > 0
    warning('Excluded %d blinks because their window fell outside recorded data boundaries.', sum(~valid_epochs));
end

eyeBlinksEpochs = eyeBlinksEpochs(valid_epochs, :);
BlinkMaxLatency = BlinkMaxLatency(valid_epochs);

% -------------------------------------------------------------------------
% 5. Find outliers
% -------------------------------------------------------------------------
% Create time axis centered around the peak (0 ms)
T = ((-winBlinksmpl:winBlinksmpl) * mspersmpl) / 1000; % Time in seconds

% Extract the blinks
num_epoch = size(eyeBlinksEpochs, 1);
N_samples = (winBlinksmpl * 2) + 1;
blinks_data = zeros(num_epoch, N_samples);
for i_epoch = 1:num_epoch
    raw_segment = eyeBlinkData(eyeBlinksEpochs(i_epoch, 1) : eyeBlinksEpochs(i_epoch, 2));
    blinks_data(i_epoch, :) = raw_segment - mean(raw_segment);
end

% Find outliers
[valid_idx, bad_mask, stats] = find_template_outliers(blinks_data', T);

eyeBlinksEpochs = eyeBlinksEpochs(valid_idx, :);
BlinkMaxLatency = BlinkMaxLatency(valid_idx);

% -------------------------------------------------------------------------
% 6. Create the Mask
% -------------------------------------------------------------------------
num_epoch = size(eyeBlinksEpochs, 1);
eyeBlinksMask = false(size(eyeBlinkData));

for i_epoch = 1:num_epoch
    eyeBlinksMask(eyeBlinksEpochs(i_epoch, 1) : eyeBlinksEpochs(i_epoch, 2)) = true;
end

% -------------------------------------------------------------------------
% 7. Visualisation
% -------------------------------------------------------------------------
if num_epoch > 0
    % Extract and baseline-correct epochs for visualization
    N_samples = (winBlinksmpl * 2) + 1;
    blinks_data = zeros(num_epoch, N_samples);

    for i_epoch = 1:num_epoch
        raw_segment = eyeBlinkData(eyeBlinksEpochs(i_epoch, 1) : eyeBlinksEpochs(i_epoch, 2));
        blinks_data(i_epoch, :) = raw_segment - mean(raw_segment);
    end

    blinks_avg = mean(blinks_data, 1);

    fh = figure('Name', 'VEOG Blinks QA', 'Color', 'w', 'Position', [100, 100, 800, 450], 'Visible', optVisible);

    % Setup compact single-tile layout to match processing standard assets
    t = tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'tight');
    ax = nexttile(t);
    hold(ax, 'on');

    % % Create time axis centered around the peak (0 ms)
    % T = ((-winBlinksmpl:winBlinksmpl) * mspersmpl) / 1000; % Time in seconds

    % Generate smooth, high-contrast density cloud palette (YlGn/Greens variation)
    colors = brewermap(num_epoch, 'YlGn');

    % Plot individual blink traces with heavy alpha transparency
    for i_epoch = 1:num_epoch
        h_line = plot(ax, T, blinks_data(i_epoch, :), 'LineWidth', 1, 'Color', colors(i_epoch, :));
        h_line.Color(4) = 0.10; % 10% opacity
    end

    % Overlay the bold Grand Average trend
    h_avg = plot(ax, T, blinks_avg, 'Color', [0.10, 0.45, 0.30], 'LineWidth', 2.5);

    % Axis Styling & Grid Typography
    grid(ax, 'on');
    set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'Layer', 'top');
    set(ax, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11);

    axis(ax, 'tight');
    xlim(ax, [T(1), T(end)]);

    % Absolute, non-overlapping formatting labels
    xlabel(ax, 'Time Relative to Blink Peak (s)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ax, 'Amplitude (\muV)', 'FontSize', 12, 'FontWeight', 'bold');

    title_str = sprintf('Detected VEOG Eyeblinks (N = %d)', num_epoch);
    title(ax, title_str, 'FontSize', 13, 'FontWeight', 'bold');

    legend(h_avg, 'Grand Average Blink', 'Location', 'NorthEast', 'Box', 'off');
    hold(ax, 'off');

    % Save
    save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_detected_veog'], [20 11]);
end




% function [eyeBlinksMask, eyeBlinksEpochs, BlinkMaxLatency, eyeBlinkData, brainData, threshold] = detect_veog(DATA, winBlink, trIQR)
%
% % winBlink = 150;
% % winSaccade = 200;
% % trIQR = 3;
%
% % -------------------------------------------------------------------------
% % Apporach 1 (RELAX)
% chaneog = strcmp({DATA.chanlocs.labels}, 'VEOG');
% assert(sum(chaneog) == 1);
%
% eyeBlinkData = DATA.data(chaneog,:);
% % eyeBlinkData = eyeBlinkData - trimmean(eyeBlinkData,10);
%
% % Temporarily filter VEOG, 0-10 Hz
% [bl, al] = butter(2, 10/(DATA.srate/2),'low');
% [bh, ah] = butter(2, 0.2/(DATA.srate/2),'high');
%
% eyeBlinkData = do_filteringcore(bl,al,eyeBlinkData,DATA.event,DATA.srate);
% eyeBlinkData = do_filteringcore(bh,ah,eyeBlinkData,DATA.event,DATA.srate);
%
% % chaneeg   = strcmp({DATA.chanlocs.type}, 'EEG');
% % brainData = DATA.data(chaneeg,:);
% % brainData = do_filteringcore(bl,al,brainData,DATA.event,DATA.srate);
% % brainData = do_filteringcore(bh,ah,brainData,DATA.event,DATA.srate);
% brainData = [];
%
% % eyeBlinkData = eyeBlinkData - trimmean(eyeBlinkData, 10);
% % eyeBlinkData = abs(eyeBlinkData);
% % figure; histogram(eyeBlinkData);
%
% % -------------------------------------------------------------------------
% % Tresholding
% EOGIQR = iqr(eyeBlinkData);
% EOG75P = prctile(eyeBlinkData, 75);
%
% % if EOGIQR > 100
% %     % If eyeblinks are too frequent and small?
% %     threshold = EOG75P + EOGIQR;
% %     fprintf('Eye blinks were detected using a threshold of 75PRC+IQR = %1.0fuV.\n',threshold);
% % else
% %     threshold = EOG75P + 1*EOGIQR;
% %     fprintf('Eye blinks were detected using a threshold of 75PRC+2IQR = %1.0fuV.\n',threshold);
% % end
% threshold = EOG75P + (trIQR * EOGIQR);
% threshold = max(threshold, 100);
% fprintf('Eye blinks (L = +-%d ms) were detected using a threshold of 75PRC + %dIQR = %1.0fuV.\n', winBlink, trIQR, threshold);
%
% % Treshold the EOG signal
% BlinkIndexMetric = double(eyeBlinkData > threshold);
%
% % figure; hold on;
% % plot(EEG.times(1:40*256),dataeog(1:40*256));
% % plot(EEG.times(1:40*256),threshold*BlinkIndexMetric(1:40*256));
%
% % -------------------------------------------------------------------------
% % Check that blinks exceed the IQR threshold for more than 50ms,
% % and to detect the blink peak within the period that exceeds the threshold:
% ix_blinkstart = find(diff(BlinkIndexMetric) == 1) + 1;  % indices where BlinkIndexMetric goes from 0 to 1
% ix_blinkend   = find(diff(BlinkIndexMetric) == -1);     % indices where BlinkIndexMetric goes from 1 to 0
%
% % # [ms] per 1 [sample]
% mspersmpl = 1000 / DATA.srate;
%
% if ~isempty(ix_blinkstart)
%     if ix_blinkend(1,1) < ix_blinkstart(1,1); ix_blinkend(:,1) = []; end % if the first downshift occurs before the upshift, remove the first value in end
%     if ix_blinkend(1,size(ix_blinkend,2)) < ix_blinkstart(1,size(ix_blinkstart,2)); ix_blinkstart(:,size(ix_blinkstart,2)) = [];end % if the last upshift occurs after the last downshift, remove the last value in start
%
%     BlinkThresholdExceededLength = ix_blinkend - ix_blinkstart; % length of consecutive samples where blink threshold was exceeded
%     BlinkRunIndex = find(BlinkThresholdExceededLength>round(50/mspersmpl)); % find locations where blink threshold was exceeded by more than 50ms
%     % find latency of the max voltage within each period where the blink
%     % threshold was exceeded:
%     if size(BlinkRunIndex,2) > 0
%         % continuousEEG.RELAX.IQRmethodDetectedBlinks = 1;
%         % epochedEEG.RELAX.IQRmethodDetectedBlinks    = 1;
%         for x = 1:size(BlinkRunIndex,2)
%             o = ix_blinkstart(BlinkRunIndex(x));
%             c = ix_blinkend(BlinkRunIndex(x));
%             [~,I] = max(eyeBlinkData(1,o:c),[],2);
%             BlinkMaxLatency(1,x) = o + I;
%         end
%     end
% end
%
% winBlinksmpl = round(winBlink/mspersmpl);
% eyeBlinksEpochs = [BlinkMaxLatency-winBlinksmpl; BlinkMaxLatency+winBlinksmpl]';
%
% % Blinks cannot be outside of the recorded data
% eyeBlinksEpochs(eyeBlinksEpochs < 1) = 1;
% eyeBlinksEpochs(eyeBlinksEpochs > DATA.pnts) = DATA.pnts;
%
% % -------------------------------------------------------------------------
% % Check if all blinks fit in the data
% L = 2 * winBlinksmpl + 1;
% NTRL = length(BlinkMaxLatency);
% EyeBlinksGood = false(NTRL,1);
%
% for i = 1:NTRL
%     if length(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2)) == L
%         EyeBlinksGood(i) = true;
%     end
% end
%
% eyeBlinksEpochs = eyeBlinksEpochs(EyeBlinksGood,:);
%
% % Report
% if ~all(EyeBlinksGood)
%     warning('There are %d blinks excluded because their window was outside of the recorded data.',sum(~EyeBlinksGood));
% end
%
% % -------------------------------------------------------------------------
% % mask = zeros(size(dataeog));
% % for i = 1:size(eyeBlinksEpochs,1)
% %     mask(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2)) = mask(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2))+1;
% % end
% % figure; hold on;
% % plot(EEG.times(274803:end),dataeog(274803:end));
% % plot(EEG.times(274803:end),threshold*mask(274803:end));
%
% % Create the mask
% eyeBlinksMask = false(size(eyeBlinkData));
% for i = 1:size(eyeBlinksEpochs,1)
%     eyeBlinksMask(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2)) = true;
% end
%
% % figure; hold on;
% % % times = EEG.times/1000; % EEGLAB time is in [ms] - > [s]
% % % plot(times(1:20*256),dataeog(1:20*256));
% % timeaxis = (0:length(dataeog))/EEG.srate;
% % plot(timeaxis(1:20*256),dataeog(1:20*256));
% % plot(locs(1:11),qrspeaks(1:11),'ro');
% %
% % figure; hold on;
% % plot(dataeog(1:30*256));
% % plot(threshold*mask(1:30*256));
%
% end