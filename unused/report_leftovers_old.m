function [DATA, flag_redo] = report_leftovers_old(DATA, tag_figure, cfg)

fprintf('\n================================\n');
fprintf('Detecting leftovers\n');
fprintf('================================\n');

%% ========================================================================
fprintf('\n--------------------------------\n');
fprintf('Muscle activity leftovers\n');
fprintf('--------------------------------\n');

emgSlopeThreshold = cfg.bch.emgSlopeThreshold;
emgSlopeDuration  = cfg.bch.emgSlopeTime;

% Estimate log-log power spectra
slopesChannelsxEpochs = detect_emg(DATA, cfg.bch);
[NCHANEEG, NTRL] = size(slopesChannelsxEpochs);

% Strong slow drifts are reflected as very steep negative slopes of the power spectrum
badchn = sum(slopesChannelsxEpochs > emgSlopeThreshold, 2);
badchn = badchn ./ NTRL;
badElectrodes = {DATA.chanlocs(find(badchn > emgSlopeDuration)).labels};

% Threshold and rank the muscle artefacts
slopesChannelsxEpochs(slopesChannelsxEpochs < emgSlopeThreshold) = NaN;
slopesChannelsxEpochs = slopesChannelsxEpochs - emgSlopeThreshold;

% Sum muscle slopes across all channels
slopesEpochs = sum(slopesChannelsxEpochs, 1, 'omitnan');
proportionOfDataShowingMuscleActivityTotal = mean(slopesEpochs > 0);

% Log Muscle Leftovers
fprintf('Total amount of leftover muscle artefact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);
fprintf(DATA.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(DATA.ALSUTRECHT.subject.fid,'Leftovers: muscle artefacts\n');
fprintf(DATA.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(DATA.ALSUTRECHT.subject.fid,'Muscle log-log slope threshold: %1.2f\n', emgSlopeThreshold);
fprintf(DATA.ALSUTRECHT.subject.fid,'Total amount of leftover muscle artefact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);

DATA.ALSUTRECHT.leftovers.muscle1 = proportionOfDataShowingMuscleActivityTotal;

% =========================================================================
% Muscle
% =========================================================================
fprintf('\n--------------------------------\n');
fprintf('Eye blink leftovers\n');
fprintf('--------------------------------\n');

% Minimum number of blinks for stats below
num_min = 3;

chanlocs = readlocs('biosemi128_eeglab.ced');
myCmap1 = brewermap(128, '*RdBu');
myCmap2 = brewermap(128, 'Reds'); % Changed to clearly show absolute amplitude

% Select only EEG
chaneeg = strcmp({DATA.chanlocs.type}, 'EEG');
dataeeg = DATA.data(chaneeg, :);

% =========================================================================
% Eye
% =========================================================================
% Detect eye blinks (Short Window)
blink_duration = 150;
blink_iqr = 3;
[~, eyeBlinksEpochs, BlinkMaxLatency, dataeog, ~, threshold] = detect_veog(DATA, blink_duration, blink_iqr, cfg.figure.visible);

% Find and remove multi-blinks
if ~isempty(eyeBlinksEpochs)
    multiBlink = detect_multiblinks(eyeBlinksEpochs, 0);
    eyeBlinksEpochs(multiBlink, :) = [];
    NTRL_1 = size(eyeBlinksEpochs, 1);
else
    NTRL_1 = 0;
end

if NTRL_1 > num_min
    L = mode(diff(eyeBlinksEpochs')) + 1;
    timeBlink0 = (0:L-1) ./ DATA.srate * 1000;
    timeBlink1_1 = timeBlink0 - blink_duration;

    dataeegepoched_1 = NaN(NCHANEEG, L, NTRL_1);
    dataeogepoched_1 = NaN(L, NTRL_1);

    for i = 1:NTRL_1
        dataeegepoched_1(:, :, i) = dataeeg(:, eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
        dataeogepoched_1(:, i)    = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    end

    baselineTime = timeBlink0(end) * [0.05 0.95];
    timesel = timeBlink0 < baselineTime(1) | timeBlink0 > baselineTime(2);

    dataeegepoched_1 = dataeegepoched_1 - mean(dataeegepoched_1(:, timesel, :), 2);
    dataeegepoched_1 = dataeegepoched_1 - mean(dataeegepoched_1, 1);

    dataeegepoched_stat = squeeze(mean(dataeegepoched_1, 2));
    [~, ~, ~, stats] = ttest(dataeegepoched_stat');
    tstatMean = mean(stats.tstat);

else
    warning('No data to make an estimate of blink leftovers (Short Window).');
    stats.tstat = NaN;
    tstatMean = NaN;
end

% =========================================================================
% Continuous Time-Series Cross-Correlation
% Correlate the continuous EEG directly with the continuous VEOG
corr_continuous = abs(corr(dataeeg', dataeog'));
meanFrontalCorr = mean(corr_continuous(ismember({DATA.chanlocs(:).labels}, cfg.ica.blinkchans)));
r_max = max(meanFrontalCorr, 0.15);

% =========================================================================
% Detect eye blinks (Long Window) - Shift to Peak-to-Peak Amplitude
blink_duration = 2000;
[~, eyeBlinksEpochs, ~, dataeog, ~, threshold] = detect_veog(DATA, blink_duration, blink_iqr, cfg.figure.visible);

if ~isempty(eyeBlinksEpochs)
    multiBlink = detect_multiblinks(eyeBlinksEpochs, 0);
    fprintf('Number of detected blinks is %d.\n', size(eyeBlinksEpochs, 1));
    fprintf('Number of detected multiple blinks within each evaluation window is %d.\n', sum(multiBlink));

    eyeBlinksEpochs(multiBlink, :) = [];
    NTRL_2 = size(eyeBlinksEpochs, 1);
else
    NTRL_2 = 0;
end

if NTRL_2 > num_min
    L = mode(diff(eyeBlinksEpochs')) + 1;
    timeBlink0 = (0:L-1) ./ DATA.srate * 1000;
    timeBlink1_2 = timeBlink0 - blink_duration;

    dataeegepoched_2 = NaN(NCHANEEG, L, NTRL_2);
    dataeogepoched_2 = NaN(L, NTRL_2);
    for i = 1:NTRL_2
        dataeegepoched_2(:, :, i) = dataeeg(:, eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
        dataeogepoched_2(:, i)   = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    end

    % Define baseline and active windows
    [~, col_500ms]  = min(abs(timeBlink0 - 500));
    [~, col_1500ms] = min(abs(timeBlink0 - 1500));
    [~, col_2500ms] = min(abs(timeBlink0 - 2500));
    [~, col_3500ms] = min(abs(timeBlink0 - 3500));
    col_4000ms = L;

    % Baseline correct data using standard windows
    dataeegepoched_2 = dataeegepoched_2 - mean(dataeegepoched_2(:, [1:col_500ms, col_3500ms:col_4000ms], :), 2);

    % Average across trials to get the ERP
    mean_blink_ERP = mean(dataeegepoched_2, 3);

    % Calculate Peak-to-Peak (PtP) in the active window
    active_window = mean_blink_ERP(:, col_1500ms:col_2500ms);
    BlinkPtP = max(active_window, [], 2) - min(active_window, [], 2);
    meanFrontalBlinkLeftoverPtP = mean(BlinkPtP(ismember({DATA.chanlocs(:).labels}, cfg.ica.blinkchans)));

else
    warning('No data to make an estimate of blink leftovers (Long Window).');
    BlinkPtP = NaN;
    meanFrontalBlinkLeftoverPtP = NaN;
end

% Log Eye Blink Leftovers
fprintf('\nAverage Frontal Peak-to-Peak Leftover: %1.2f uV\n', meanFrontalBlinkLeftoverPtP);
fprintf('Average Frontal VEOG-EEG Correlation: r = %1.2f\n', meanFrontalCorr);

fprintf(DATA.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(DATA.ALSUTRECHT.subject.fid,'Leftovers: eye blink artefacts\n');
fprintf(DATA.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(DATA.ALSUTRECHT.subject.fid,'Average Frontal Peak-to-Peak Leftover: %1.2f uV\n', meanFrontalBlinkLeftoverPtP);
fprintf(DATA.ALSUTRECHT.subject.fid,'Average Frontal VEOG-EEG Correlation: r = %1.2f\n', meanFrontalCorr);

DATA.ALSUTRECHT.leftovers.blinksPtP = BlinkPtP;
DATA.ALSUTRECHT.leftovers.blinksTstat = stats.tstat;
DATA.ALSUTRECHT.leftovers.corrContinuous = corr_continuous;

% =========================================================================
% Plot
% =========================================================================
close all hidden;
fh = figure('Visible', cfg.figure.visible);
tiledlayout(3, 2, "TileSpacing", "compact", "Padding", "compact");

if NTRL_1 > num_min
    % Plot 1: EOG Traces
    nexttile(1); hold on;
    plot(timeBlink1_1, 0 * ones(size(timeBlink1_1)), 'Color', 0 * [1 1 1], 'LineStyle', '-');
    plot(timeBlink1_1, dataeogepoched_1, 'LineWidth', 1.2);
    plot(timeBlink1_1, threshold * ones(size(timeBlink1_1)), 'Color', 0.5 * [1 1 1], 'LineStyle', '--');
    set(gca, 'ColorOrder', [0 0 0; brewermap(NTRL_1, 'Spectral'); 0.5 * [1 1 1]]);
    title(['Detected blinks, N = ' num2str(NTRL_1)]);
    pbaspect([1.618 1 1]); ylabel('EOG amplitude (\muV)');

    % Plot 2: Topoplot of T-statistics
    nexttile(2); hold on;
    topoplot(stats.tstat, chanlocs, 'maplimits', max(abs(stats.tstat))*[-1 1], 'headrad', 0.5, 'colormap', myCmap1, 'whitebk', 'on', 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    title('EEG timelocked to the eye blinks');
    hcb = colorbar;
    hcb.Title.String = "T-value";
end

if NTRL_2 > num_min
    % Plot 3: Long Window EOG Traces
    nexttile(3); hold on;
    plot(timeBlink1_2, 0 * ones(size(timeBlink1_2)), 'Color', 0 * [1 1 1], 'LineStyle', '-');
    plot(timeBlink1_2, dataeogepoched_2, 'LineWidth', 1.2);
    plot(timeBlink1_2, threshold * ones(size(timeBlink1_2)), 'Color', 0.5 * [1 1 1], 'LineStyle', '--');
    set(gca, 'ColorOrder', [0 0 0; brewermap(NTRL_2, 'Spectral'); 0.5 * [1 1 1]]);
    title(['Detected blinks, N = ' num2str(NTRL_2)]);
    pbaspect([1.618 1 1]); ylabel('EOG (\muV)');

    % Plot 4: Topoplot of Peak-to-Peak Amplitude
    nexttile(4);
    maxBlinkPtP = prctile(BlinkPtP, 95);
    maxBlinkPtP = max(maxBlinkPtP, 2.0); % Ensure colour bar is scaled nicely
    topoplot(BlinkPtP, chanlocs, 'maplimits', [0 maxBlinkPtP], 'headrad', 0.5, 'colormap', myCmap2, 'whitebk', 'on', 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    title({'Peak-to-Peak Amplitude', [num2str(round(meanFrontalBlinkLeftoverPtP, 1)) ' \muV Frontal']});
    hcb = colorbar;
    hcb.Title.String = "\muV";

    % Plot 6: Frontal ERP Trace (Visual validation of the PtP amplitude)
    nexttile(6); hold on;
    maskChanBlink = ismember({DATA.chanlocs(:).labels}, cfg.ica.blinkchans);
    plot(timeBlink1_2, mean_blink_ERP(maskChanBlink, :), 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1);
    plot(timeBlink1_2, mean(mean_blink_ERP(maskChanBlink, :), 1), 'k', 'LineWidth', 2);
    xline(timeBlink1_2(col_1500ms), '--k'); xline(timeBlink1_2(col_2500ms), '--k');
    title('Frontal Channels ERP');
    ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
    pbaspect([1.618 1 1]);
end

nexttile(5);
topoplot(corr_continuous, chanlocs, 'maplimits', [0 r_max], 'headrad', 0.5, 'colormap', myCmap2, 'whitebk', 'on', 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
title({'Continuous VEOG-EEG Correlation', ['Frontal Mean r = ' num2str(round(meanFrontalCorr, 2))]});
hcb = colorbar;
hcb.Title.String = "Pearson r";

% Save Figure
% plotX = 25; plotY = 35;
% set(fh, 'InvertHardCopy', 'Off', 'Color', [1 1 1]);
% set(fh, 'PaperPositionMode', 'Manual', 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 plotX plotY], 'PaperSize', [plotX plotY]);
% print(fh, fullfile(DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_leftovers_' num2str(tag_figure)]), '-dtiff', '-r200'); close(fh);
save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_leftovers_' num2str(tag_figure)], [25 35]);

% =========================================================================
% Compare Leftovers to ICA Templates & Generate Final Report
% =========================================================================
templatesICA = load_ictemplateweights(DATA);
maskChanBlink = ismember({DATA.chanlocs(:).labels}, cfg.ica.blinkchans);

% Evaluation Thresholds
blinkThresholdCorr   = 0.45; % Spatial correlation threshold
blinkThreshold2Tstat = 4.0;  % Statistical significance threshold
blinkThresholdPtP    = 3.5;  % Physical amplitude threshold (uV)
blinkThresholdContr  = 0.3;  % Continuous correlation threshold (Pearson r)

% Initialise flags
flag_redo = false;
flag_CorrTstat = NaN; flag_CorrPtP = NaN;
flag_MeanTstat = NaN; flag_MeanPtP = NaN; flag_ContCorr = NaN;
corrMatTstat  = NaN; corrMatPtP  = NaN;
meanFrontalBlinkLeftoverTstat = NaN;

% 1. Calculate Metrics if data is available
if ~isnan(tstatMean) && NTRL_1 > num_min
    % Extract the raw weights
    blinkWeights = templatesICA.Blinkweights0;
    blinkLeftover = stats.tstat';

    % --- NEW: Ensure Template Polarity is Frontal-Positive ---
    % Because ICA signs are arbitrary, we flip the template so that the
    % average of your 3 core blink channels is always positive.
    if mean(blinkWeights(maskChanBlink)) < 0
        blinkWeights = -blinkWeights;
    end
    % ---------------------------------------------------------

    % Normalise and apply Laplacian as before
    Blinkweights0Norm = estimate_invlaplacian(blinkWeights ./ norm(blinkWeights), DATA.chanlocs, 1);
    BlinkAmplitudeTstatNorm = estimate_invlaplacian(blinkLeftover ./ norm(blinkLeftover), DATA.chanlocs, 1);

    % % Use the C-bundle (Anterior Mask) as discussed
    % chan_labels = {DATA.chanlocs(:).labels};
    % maskAnterior = startsWith(chan_labels, 'C');

    % Restrict the correlation to the anterior half of the head
    anterior_Tstat = BlinkAmplitudeTstatNorm;
    anterior_Weights = Blinkweights0Norm;

    % --- NEW: Strictly Positive Correlation ---
    % We removed abs() here. Now, only a positive match (blink leakage)
    % triggers a fail. An inverted match (like posterior activity) is ignored.
    corrMatTstat = corr(anterior_Tstat, anterior_Weights);

    % Evaluate Flags: Using a strictly positive threshold
    flag_CorrTstat = corrMatTstat > blinkThresholdCorr;

    meanFrontalBlinkLeftoverTstat = mean(stats.tstat(maskChanBlink));
    flag_MeanTstat = meanFrontalBlinkLeftoverTstat > blinkThreshold2Tstat;

    if ~isnan(meanFrontalBlinkLeftoverPtP) && NTRL_2 > num_min
        % Repeat flipping logic for PtP if necessary, though usually
        % the template flip above is sufficient for the whole block.
        BlinkPtPNorm = estimate_invlaplacian(BlinkPtP ./ norm(BlinkPtP), DATA.chanlocs, 1);
        anterior_PtP = BlinkPtPNorm;

        corrMatPtP = corr(anterior_PtP, anterior_Weights);

        % Evaluate Flags
        flag_CorrPtP  = corrMatPtP > blinkThresholdCorr;
        flag_MeanPtP  = meanFrontalBlinkLeftoverPtP > blinkThresholdPtP;
        flag_ContCorr = meanFrontalCorr > blinkThresholdContr;

        % Final REDO Logic
        flag_redo = (flag_CorrTstat | flag_CorrPtP) & (flag_MeanTstat | flag_MeanPtP | flag_ContCorr);
    end
end

% Save flags to DATA structure
DATA.ALSUTRECHT.leftovers.flag_redo     = flag_redo;
DATA.ALSUTRECHT.leftovers.flagCorrTstat = flag_CorrTstat;
DATA.ALSUTRECHT.leftovers.flagCorrPtP   = flag_CorrPtP;
DATA.ALSUTRECHT.leftovers.flagMeanTstat = flag_MeanTstat;
DATA.ALSUTRECHT.leftovers.flagMeanPtP   = flag_MeanPtP;

% =========================================================================
% 2. Console and Log File Reporting
% =========================================================================
fprintf('\n=========================================================\n');
fprintf('FINAL BLINK LEFTOVER QUALITY ASSURANCE REPORT\n');
fprintf('=========================================================\n');

if isnan(tstatMean) || NTRL_1 <= num_min
    fprintf('STATUS: [WARNING] - Insufficient blinks detected to run evaluation.\n');
elseif flag_redo
    fprintf('STATUS: [FAIL] - Significant blink leakage detected. ICA redo recommended.\n');
else
    fprintf('STATUS: [PASS] - Blink cleaning is within acceptable physiological limits.\n');
end

fprintf('\n--- Spatial Template Matching ---\n');
fprintf('%s T-Stat Map Correlation : %1.2f (Threshold: %1.2f)\n', get_status_str(flag_CorrTstat), corrMatTstat, blinkThresholdCorr);
fprintf('%s PtP Map Correlation    : %1.2f (Threshold: %1.2f)\n', get_status_str(flag_CorrPtP), corrMatPtP, blinkThresholdCorr);

fprintf('\n--- Amplitude & Statistical Leftovers (Frontal) ---\n');
fprintf('%s Average T-Statistic    : %1.2f (Threshold: %1.2f)\n', get_status_str(flag_MeanTstat), meanFrontalBlinkLeftoverTstat, blinkThreshold2Tstat);
fprintf('%s Average Peak-to-Peak   : %1.2f uV (Threshold: %1.2f uV)\n', get_status_str(flag_MeanPtP), meanFrontalBlinkLeftoverPtP, blinkThresholdPtP);
fprintf('%s Continuous Correlation : r = %1.2f (Threshold: %1.2f)\n', get_status_str(flag_ContCorr), meanFrontalCorr, blinkThresholdContr);
fprintf('=========================================================\n\n');

% Print identical summary to the subject log file
fid = DATA.ALSUTRECHT.subject.fid;
fprintf(fid, '\n---------------------------------------------------------\n');
fprintf(fid, 'Leftovers: Final Blink Quality Assurance Report\n');
fprintf(fid, '---------------------------------------------------------\n');
if flag_redo
    fprintf(fid, 'STATUS: [FAIL] - Significant blink leakage detected.\n');
else
    fprintf(fid, 'STATUS: [PASS] - Blink cleaning acceptable.\n');
end
fprintf(fid, 'T-Stat Map Correlation : %1.2f\n', corrMatTstat);
fprintf(fid, 'PtP Map Correlation    : %1.2f\n', corrMatPtP);
fprintf(fid, 'Average T-Statistic    : %1.2f\n', meanFrontalBlinkLeftoverTstat);
fprintf(fid, 'Average Peak-to-Peak   : %1.2f uV\n', meanFrontalBlinkLeftoverPtP);
fprintf(fid, 'Continuous Correlation : r = %1.2f\n', meanFrontalCorr);

end

% =========================================================================
% Helper functions
% =========================================================================
function multiBlink = detect_multiblinks(eyeBlinksEpochs, overlap)
eyeBlinksEpochs(:, 1) = eyeBlinksEpochs(:, 1) + overlap;
eyeBlinksEpochs(:, 2) = eyeBlinksEpochs(:, 2) - overlap;

NTRL = size(eyeBlinksEpochs, 1);
multiBlink = false(NTRL, 1);

starts = eyeBlinksEpochs(:, 1);
ends   = eyeBlinksEpochs(:, 2);

for i = 1:NTRL
    overlaps = (starts(i) <= ends) & (ends(i) >= starts);
    overlaps(i) = false;

    if any(overlaps)
        multiBlink(i) = true;
    end
end
end

function str = get_status_str(flag)
% Helper function to generate clean PASS/FAIL status tags
if isnan(flag)
    str = '[N/A]';
elseif flag
    str = '[FAIL]';
else
    str = '[PASS]';
end
end