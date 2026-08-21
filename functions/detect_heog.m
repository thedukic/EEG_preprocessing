function [saccadesMask, saccadesEpochs, saccadesMaxLatency, saccadesData, brainData, threshold] = detect_heog(DATA, winSaccade, trIQRsaccade, optVisible)

% -------------------------------------------------------------------------
% Initialize Outputs (Prevents crashes if no saccades are found)
% -------------------------------------------------------------------------
saccadesMask       = [];
saccadesEpochs     = [];
saccadesMaxLatency = [];
brainData          = [];
threshold          = NaN;

winBlink    = 200;
trIQRblink  = 1;
do_combined = false;

% -------------------------------------------------------------------------
% Extract and Filter HEOG
% -------------------------------------------------------------------------
chanheog = strcmp({DATA.chanlocs.labels}, 'HEOG');
assert(sum(chanheog) == 1, 'Could not uniquely identify HEOG channel.');
saccadesData = DATA.data(chanheog,:);

% Temporarily filter HEOG, 0.2-30 Hz
[bl, al] = butter(2, 30/(DATA.srate/2), 'low');
[bh, ah] = butter(2, 0.2/(DATA.srate/2), 'high');
saccadesData = do_filteringcore(bl, al, saccadesData, DATA.event, DATA.srate);
saccadesData = do_filteringcore(bh, ah, saccadesData, DATA.event, DATA.srate);

% -------------------------------------------------------------------------
% HEOG Thresholding
% -------------------------------------------------------------------------
EOGIQR = iqr(saccadesData);
EOG75P = prctile(saccadesData, 75);
threshold = EOG75P + (trIQRsaccade * EOGIQR);
fprintf('Eye saccades (L = +-%d ms) were detected using a threshold of 75PRC + %dIQR = %1.0fuV.\n', winSaccade, trIQRsaccade, threshold);

saccadesMaskTmp = saccadesData > threshold;

if ~any(saccadesMaskTmp)
    warning('No activity exceeded the HEOG threshold.');
    return; % Early exit
end

% -------------------------------------------------------------------------
% 1. Exclude False Positives using VEOG
% -------------------------------------------------------------------------
fprintf('\nEye blinks (L = +-%d ms) from VEOG will be used to remove false saccades in HEOG.\n', winBlink);
[noiseMaskBlink, ~] = detect_veog(DATA, winBlink, trIQRblink, optVisible);

jump = find(diff([false, saccadesMaskTmp, false]) ~= 0);
durall = jump(2:2:end) - jump(1:2:end);

% Minimum of 20 ms of HEOG duration
mspersamp = 1000/DATA.srate;
mindrftdur = round(20/mspersamp);
saccadesEpochsTmp = find(durall > mindrftdur);
num_heog = length(saccadesEpochsTmp);

if num_heog == 0
    warning('No saccades were detected longer than 20ms...');
    return;
end

jumpStart = jump(1:2:end);
jumpStop  = jump(2:2:end);
EOGfocussamples = round(winSaccade/mspersamp);

jumpStart = jumpStart - EOGfocussamples;
jumpStop = jumpStop + EOGfocussamples;

% Boundary checks
jumpStart(jumpStart < 1) = 1;
jumpStop(jumpStop > DATA.pnts) = DATA.pnts;

% Remove those that are actually blinks captured by HEOG
actuallyBlinks = false(num_heog, 1);
for i = 1:num_heog
    idx = saccadesEpochsTmp(i);
    actuallyBlinks(i) = any(noiseMaskBlink(jumpStart(idx):jumpStop(idx)));
end
saccadesEpochsTmp(actuallyBlinks) = [];
num_heog = length(saccadesEpochsTmp);

if num_heog == 0
    warning('All detected saccades were flagged as blinks...');
    return;
end

% -------------------------------------------------------------------------
% 2. Template Matching via Cross-Correlation
% -------------------------------------------------------------------------
load('modelSaccade.mat', 'modelSaccade');
[saccadesData_norm, L_current] = interp_saccades(saccadesData, jumpStart, jumpStop, saccadesEpochsTmp);

saccadesData_norm = zscore(saccadesData_norm);

% Correlate with the models/priors
[r1, ~] = corr(saccadesData_norm, modelSaccade(:,1));
[r2, ~] = corr(saccadesData_norm, modelSaccade(:,2));

% Detect the true saccades
Rmin = 0.7;
indx1 = r1 > Rmin;
indx2 = r2 > Rmin;
indx = indx1 & indx2;

% Remove ambiguous mixed-direction saccades
if any(indx)
    indx1(indx) = false;
    indx2(indx) = false;
end
assert(~any(indx1 & indx2));

indx = indx1 | indx2;

if ~any(indx)
    warning('\nNo HEOG events matched the saccade template (R > 0.7)...');
    return;
else

    fprintf('\nNumber of saccades: %d\n', num_heog);
end

% -------------------------------------------------------------------------
% Save Grand Average Contribution
% -------------------------------------------------------------------------
% avg1 = mean(saccadesData_norm(:,indx1), 2);
% avg2 = mean(saccadesData_norm(:,indx2), 2);
%
% save_path = fullfile(EEG.ALSUTRECHT.subject.mycodes, 'files', 'weights', 'modelsaccadestmp');
% if exist(save_path, 'dir')
%     save(fullfile(save_path, [DATA.ALSUTRECHT.subject.id '_saccade_avg.mat']), "avg1", "avg2");
% end

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------
fh = figure('Name', 'HEOG Saccades QA', 'Color', 'w', 'Position', [100, 100, 800, 450], 'Visible', optVisible);

% Setup a compact single-tile layout to eliminate thick white space borders
t = tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'tight');
ax = nexttile(t);
hold(ax, 'on');

T = linspace(-0.5, 0.5, size(saccadesData_norm, 1));

% Refined, earthy color palette for left/right saccades
% Primary Red/Burgundy for class 1, Deep Slate Blue/Teal for class 2
c_heog = [0.65, 0.15, 0.15; ...  % Class 1 (e.g., Left)
    0.15, 0.40, 0.55];     % Class 2 (e.g., Right)

% Step 1: Plot individual traces with high transparency
% This allows hundreds of step functions to cluster cleanly into a density map
for i = 1:num_heog
    y = saccadesData_norm(:, i);
    if indx1(i)
        h_line = plot(ax, T, y, 'LineWidth', 1, 'Color', c_heog(1, :));
        h_line.Color(4) = 0.08; % 8% opacity
    elseif indx2(i)
        h_line = plot(ax, T, y, 'LineWidth', 1, 'Color', c_heog(2, :));
        h_line.Color(4) = 0.08; % 8% opacity
    end
end

% Step 2: Plot Grand Averages safely on top
% Enforcing a clean 3.0 LineWidth keeps the trend lines bold without cluttering the canvas
h_plots = [];
legend_labels = {};

if any(indx1)
    avg1 = mean(saccadesData_norm(:, indx1), 2);
    h1 = plot(ax, T, avg1, 'LineWidth', 3.0, 'Color', c_heog(1, :));
    h_plots = [h_plots, h1];
    legend_labels = [legend_labels, sprintf('Direction 1 (N = %d)', sum(indx1))];
end

if any(indx2)
    avg2 = mean(saccadesData_norm(:, indx2), 2);
    h2 = plot(ax, T, avg2, 'LineWidth', 3.0, 'Color', c_heog(2, :));
    h_plots = [h_plots, h2];
    legend_labels = [legend_labels, sprintf('Direction 2 (N = %d)', sum(indx2))];
end

% Step 3: Refined Axis Styling & Grid Typography
grid(ax, 'on');
set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'Layer', 'top');
set(ax, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 11);

axis(ax, 'tight');
xlim(ax, [T(1), T(end)]);

% Clear, non-overlapping labels
xlabel(ax, 'Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax, 'HEOG Amplitude (Z-Score)', 'FontSize', 12, 'FontWeight', 'bold');

% Clean Title incorporating metadata dynamically
title_str = sprintf('Detected HEOG Saccades (N = %d)', sum(indx1) + sum(indx2));
title(ax, title_str, 'FontSize', 13, 'FontWeight', 'bold');

% Add legend pointing only to the clean average traces
if ~isempty(h_plots)
    legend(h_plots, legend_labels, 'Location', 'NorthEast', 'Box', 'off');
end

hold(ax, 'off');

% Save configuration using your standard resolution settings
% plotX = 20; plotY = 11;
% set(fh, 'InvertHardCopy', 'Off', 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
% print(fh, fullfile(DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_detected_heog']), '-dtiff', '-r150'); close(fh);
save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_detected_heog'], [20 11]);

% -------------------------------------------------------------------------
% Finalize Outputs (Calculate Max Latencies)
% -------------------------------------------------------------------------
TEOG = mean(L_current(indx)) / DATA.srate * 1000;
fprintf('Average saccade duration: %2.0f ms.\n', TEOG);

valid_saccades = saccadesEpochsTmp(indx);
N_FINAL = length(valid_saccades);
saccadesMaxLatency = zeros(1, N_FINAL);

for i = 1:N_FINAL
    idx = valid_saccades(i);
    [~, peak_offset] = max(abs(saccadesData(jumpStart(idx):jumpStop(idx))));
    saccadesMaxLatency(i) = jumpStart(idx) + peak_offset - 1;
end

if do_combined
    fprintf('L/R saccades will be combined in the output.\n');
    saccadesEpochs = [jumpStart(valid_saccades); jumpStop(valid_saccades)]';
    saccadesMask = false(size(saccadesData));
    for i = 1:N_FINAL
        saccadesMask(saccadesEpochs(i,1):saccadesEpochs(i,2)) = true;
    end
else
    fprintf('L/R saccades will be separated in the output.\n');
    saccadesEpochsTmp_lr = {saccadesEpochsTmp(indx1), saccadesEpochsTmp(indx2)};
    saccadesEpochs = cell(1,2);
    saccadesMask = cell(1,2);
    for i_lr = 1:2
        N_LR = length(saccadesEpochsTmp_lr{i_lr});
        saccadesMask{i_lr} = false(size(saccadesData));
        saccadesEpochs{i_lr} = [jumpStart(saccadesEpochsTmp_lr{i_lr}); jumpStop(saccadesEpochsTmp_lr{i_lr})]';
        for i = 1:N_LR
            saccadesMask{i_lr}(saccadesEpochs{i_lr}(i,1) : saccadesEpochs{i_lr}(i,2)) = true;
        end
    end
end

end

% =========================================================================
% Helper function
% =========================================================================
function [saccadesData_norm, L_current] = interp_saccades(saccadesData, jumpStart, jumpStop, saccadesEpochs)
L_target = 512;
X_target = linspace(1, L_target, L_target);
NHEOG = length(saccadesEpochs);
saccadesData_norm = NaN(L_target, NHEOG);
L_current = NaN(NHEOG,1);

for i = 1:NHEOG
    current_signal = saccadesData(jumpStart(saccadesEpochs(i)):jumpStop(saccadesEpochs(i)));
    L_current(i) = length(current_signal);
    X_current = linspace(1, L_target, L_current(i));
    % Changed to 'pchip' to prevent artificial splining artifacts
    Y_interpolated = interp1(X_current, current_signal, X_target, 'pchip');
    saccadesData_norm(:, i) = Y_interpolated(:);
end
end