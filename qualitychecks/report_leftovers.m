function [DATA, flag_redo] = report_leftovers(DATA, tag_figure, cfg)

fprintf('\n================================\n');
fprintf('Detecting leftovers\n');
fprintf('================================\n');

% =========================================================================
% Muscle
% =========================================================================
fprintf('\n--------------------------------\n');
fprintf('Muscle activity leftovers evaluation\n');
fprintf('--------------------------------\n');

emgSlopeThreshold = cfg.emg.slope_threshold_1;
emgSlopeDuration  = cfg.emg.slope_time;

% Estimate log-log power spectra
slopesChannelsxEpochs = detect_emg(DATA, cfg);
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
% Eye
% =========================================================================
fprintf('\n--------------------------------\n');
fprintf('Eye blink leftovers evaluation\n');
fprintf('--------------------------------\n');

% Minimum number of blinks for stats below
num_trial_min = 5;

% Select only EEG
chaneeg = strcmp({DATA.chanlocs.type}, 'EEG');
dataeeg = DATA.data(chaneeg, :);

% Define the broad anterior search zone using the C-bundle
% mask_frontal = ismember({DATA.chanlocs(:).labels}, cfg.ica.blinkchans);
% mask_frontal = startsWith({DATA.chanlocs.labels}, 'C');

% Extract the X-coordinates for all channels (EEGLAB: +X is towards the nose)
x_coords = [DATA.chanlocs.X];

% Define the frontal zone as the anterior third of the head
% (i.e., any channel sitting further forward than half the maximum forward distance)
frontal_cutoff = max(x_coords) * 0.5;
mask_frontal = (x_coords > frontal_cutoff);
indx_frontal = find(mask_frontal);

% =========================================================================
% Eye (Short window)
% =========================================================================
fprintf('\n#### Eye blink leftovers evaluation (Short window) ####\n');

blink_duration = 150; % Yields a 300 ms total window (-150 to +150 ms)
blink_iqr = 3;
[~, eyeBlinksEpochs, BlinkMaxLatency, dataeog, ~, threshold] = detect_veog(DATA, blink_duration, blink_iqr, cfg.figure.visible);

if ~isempty(eyeBlinksEpochs)
    multiBlink = detect_multiblinks(eyeBlinksEpochs, 0);
    eyeBlinksEpochs(multiBlink, :) = [];
    NTRL_1 = size(eyeBlinksEpochs, 1);
else
    NTRL_1 = 0;
end

if NTRL_1 > num_trial_min
    % Force consistent length based on sampling rate (2 x 150 ms = 300 ms window)
    smp_window = round((blink_duration * 2 / 1000) * DATA.srate);
    L = smp_window + 1;

    % Time vector from -150 to +150 ms (retained for downstream plotting)
    timeBlink1_1 = linspace(-blink_duration, blink_duration, L);

    dataeegepoched_1 = NaN(NCHANEEG, L, NTRL_1);

    % Restored and initialised for plotting
    dataeogepoched_1 = NaN(L, NTRL_1);

    for i = 1:NTRL_1
        idx_start = eyeBlinksEpochs(i,1);
        idx_end   = idx_start + smp_window;

        % Safety boundary check
        if idx_start >= 1 && idx_end <= size(dataeeg, 2)
            dataeegepoched_1(:, :, i) = dataeeg(:, idx_start:idx_end);

            % Restored EOG extraction
            dataeogepoched_1(:, i)    = dataeog(idx_start:idx_end);
        end
    end

    % 1. Local Detrending (Baseline)
    % A 300 ms window lacks a true resting baseline. We use the outer 10%
    % (first 5% and last 5%) strictly as a local detrending measure.
    idx_baseline = (timeBlink1_1 < timeBlink1_1(round(L * 0.05))) | (timeBlink1_1 > timeBlink1_1(round(L * 0.95)));

    % Baseline correct EEG
    baseline_means = mean(dataeegepoched_1(:, idx_baseline, :), 2, 'omitnan');
    dataeegepoched_1 = dataeegepoched_1 - baseline_means;

    % Restored EOG baseline correction
    eog_baseline_means = mean(dataeogepoched_1(idx_baseline, :), 1, 'omitnan');
    dataeogepoched_1 = dataeogepoched_1 - eog_baseline_means;

    % 2. Target the peak for the T-test
    % Extract a narrow 50 ms window (-25 to +25 ms) centered exactly on the blink peak.
    idx_peak = (timeBlink1_1 >= -25 & timeBlink1_1 <= 25);
    peak_amplitudes = squeeze(mean(dataeegepoched_1(:, idx_peak, :), 2, 'omitnan'));

    % 3. Compute T-test across trials (tests if the peak amplitude != 0)
    [~, ~, ~, stats] = ttest(peak_amplitudes');


    % 4. Apply the spatial masks
    anterior_tstats  = stats.tstat(mask_frontal);
    posterior_tstats = stats.tstat(~mask_frontal);

    % Extract the top 5 highest t-statistics STRICTLY from the frontal channels
    meanFrontalTstat = mean(maxk(anterior_tstats, 5));

    % Calculate the mean T-statistic for the rest of the scalp
    meanPosteriorTstat = mean(posterior_tstats, 'omitnan');

    % Calculate the Spatial Gradient for the T-test
    SpatialGradient_Tstat = meanFrontalTstat - meanPosteriorTstat;

else
    fprintf('Warning: No data (N = %d) to make an estimate of blink leftovers (Short Window).\n', NTRL_1);
    stats.tstat = NaN;
    meanFrontalTstat = NaN;
end

% fprintf('#### Eye blink leftovers (Short window) ####\n');
%
% % Detect eye blinks
% blink_duration = 150;
% blink_iqr = 3;
% [~, eyeBlinksEpochs, BlinkMaxLatency, dataeog, ~, threshold] = detect_veog(DATA, blink_duration, blink_iqr, cfg.figure.visible);
%
% % Find and remove multi-blinks
% if ~isempty(eyeBlinksEpochs)
%     multiBlink = detect_multiblinks(eyeBlinksEpochs, 0);
%     eyeBlinksEpochs(multiBlink, :) = [];
%     NTRL_1 = size(eyeBlinksEpochs, 1);
% else
%     NTRL_1 = 0;
% end
%
% if NTRL_1 > num_trial_min
%     L = mode(diff(eyeBlinksEpochs')) + 1;
%     timeBlink0 = (0:L-1) ./ DATA.srate * 1000;
%     timeBlink1_1 = timeBlink0 - blink_duration;
%
%     dataeegepoched_1 = NaN(NCHANEEG, L, NTRL_1);
%     dataeogepoched_1 = NaN(L, NTRL_1);
%
%     for i = 1:NTRL_1
%         dataeegepoched_1(:, :, i) = dataeeg(:, eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
%         dataeogepoched_1(:, i)    = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
%     end
%
%     baselineTime = timeBlink0(end) * [0.05 0.95];
%     timesel = timeBlink0 < baselineTime(1) | timeBlink0 > baselineTime(2);
%
%     % Common-average and baseline correction
%     dataeegepoched_1 = dataeegepoched_1 - mean(dataeegepoched_1, 1);
%     dataeegepoched_1 = dataeegepoched_1 - mean(dataeegepoched_1(:, timesel, :), 2);
%     dataeogepoched_1 = dataeogepoched_1 - mean(dataeogepoched_1(timesel, :), 1);
%
%     % T- test if the blink peak is different from zero
%     dataeegepoched_stat = squeeze(mean(dataeegepoched_1, 2));
%     [~, ~, ~, stats] = ttest(dataeegepoched_stat');
%     meanFrontalTstat = mean(maxk(stats.tstat, 5));
% else
%     fprintf('Warning: No data (N = %d) to make an estimate of blink leftovers (Short Window).\n', NTRL_1);
%     stats.tstat = NaN;
%     meanFrontalTstat = NaN;
% end

% =========================================================================
% VEOG-EEG correlation
% =========================================================================
corr_continuous = abs(corr(dataeeg', dataeog', "Type", "Spearman"));

% Isolate correlation values for the C-bundle
anterior_corr = corr_continuous(mask_frontal);

% Extract the mean of the top 4 highest correlated anterior channels
% This catches asymmetrical tracking without diluting the metric
meanFrontalCorr_cont = mean(maxk(anterior_corr, 5));

% =========================================================================
% Eye (Long window)
% =========================================================================
fprintf('\n#### Eye blink leftovers evaluation (Long window) ####\n');

% Shorten window
blink_duration = 1000; % Yields a 2000 ms total window (-1000 to +1000 ms)
[~, eyeBlinksEpochs, ~, dataeog, ~, threshold] = detect_veog(DATA, blink_duration, blink_iqr, cfg.figure.visible);

if ~isempty(eyeBlinksEpochs)
    multiBlink = detect_multiblinks(eyeBlinksEpochs, 0);

    fprintf('Detected blinks: %d\n', size(eyeBlinksEpochs, 1));
    fprintf('Multi-blinks rejected (overlapping within window): %d\n', sum(multiBlink));

    eyeBlinksEpochs(multiBlink, :) = [];
    NTRL_2 = size(eyeBlinksEpochs, 1);
else
    NTRL_2 = 0;
end

if NTRL_2 > num_trial_min
    % Force consistent length based on sampling rate to prevent mode(diff) edge cases
    L = 2 * round((blink_duration / 1000) * DATA.srate);

    % timeBlink0 = (0:L-1) ./ DATA.srate * 1000;
    % timeBlink1_2 = timeBlink0 - blink_duration;

    % Create a time vector centered around the blink peak (0 ms)
    timeBlink1_2 = linspace(-blink_duration/2, blink_duration/2, L);

    dataeegepoched_2 = NaN(NCHANEEG, L, NTRL_2);
    dataeogepoched_2 = NaN(L, NTRL_2);

    for i = 1:NTRL_2
        idx_start = eyeBlinksEpochs(i, 1);
        idx_end   = eyeBlinksEpochs(i, 2) - 1;

        % Safety boundary check (ensures we don't index outside the recording)
        if idx_start >= 1 && idx_end <= size(dataeeg, 2)
            dataeegepoched_2(:, :, i) = dataeeg(:, idx_start:idx_end);
            dataeogepoched_2(:, i)    = dataeog(idx_start:idx_end);
        end
    end

    % Define time window parameters (in ms)
    win_baseline_pre  = [-500, -250];
    win_active        = [-150,  250];
    win_baseline_post = [ 250,  500];

    % Define dynamic logical windows
    idx_active = (timeBlink1_2 >= win_active(1) & timeBlink1_2 <= win_active(2));

    idx_baseline = (timeBlink1_2 >= win_baseline_pre(1)  & timeBlink1_2 <= win_baseline_pre(2)) | ...
        (timeBlink1_2 >= win_baseline_post(1) & timeBlink1_2 <= win_baseline_post(2));

    % 1. Take the absolute values of the epoched data
    abs_eeg = abs(dataeegepoched_2);

    % 2. Calculate the mean amplitude in the active and baseline windows
    mean_active   = mean(abs_eeg(:, idx_active, :), 2, 'omitnan');
    mean_baseline = mean(abs_eeg(:, idx_baseline, :), 2, 'omitnan');

    % 3. Calculate the amplitude ratio for all channels and all epochs simultaneously
    BlinkAmplitudeRatioAllEpochs = mean_active ./ mean_baseline;

    % 4. Average the ratio across all trials to get one overall value per channel
    MeanBlinkRatio = mean(BlinkAmplitudeRatioAllEpochs, 3, 'omitnan');

    % 5. Extract the top 5 worst channels based on this ratio
    [worst_ratios, worst_local_indices] = maxk(MeanBlinkRatio(mask_frontal), 5);
    ptp_worst_global_idx = indx_frontal(worst_local_indices);
    meanFrontalRatio = mean(worst_ratios);

    % Console verification
    worst_labels = {DATA.chanlocs(ptp_worst_global_idx).labels};
    fprintf('Top 5 worst channels selected for audit: %s\n', strjoin(worst_labels, ', '));
    fprintf('Mean Frontal Amplitude Ratio: %.2f (Ideal is ~1.0)\n', meanFrontalRatio);

    % =====================================================================
    % ERP-EOG correlation in the active window
    % =====================================================================
    % Baseline correct data using the two-sided baseline
    baseline_means = mean(dataeegepoched_2(:, idx_baseline, :), 2, 'omitnan');
    dataeegepoched_2 = dataeegepoched_2 - baseline_means;

    % Average across trials to get the ERP
    mean_blink_ERP = mean(dataeegepoched_2, 3, 'omitnan');
    mean_blink_ERP_avg = mean(mean_blink_ERP, 1)';

    worst_erp_mean = mean(mean_blink_ERP(ptp_worst_global_idx, idx_active), 1)';
    eog_mean_active = mean(dataeogepoched_2(idx_active, :), 2, 'omitnan');

    FrontalCorr_erp = abs(corr(worst_erp_mean, eog_mean_active, 'Type', 'Spearman', 'Rows', 'complete'));

else
    fprintf('Warning: No data (N = %d) to make an estimate of blink leftovers.\n', NTRL_2);
    MeanBlinkRatio = NaN;
    meanFrontalRatio = NaN;
    FrontalCorr_erp = NaN;
end

% fprintf('#### Eye blink leftovers (Long window) ####\n');
%
% % Detect eye blinks
% blink_duration = 2000;
% [~, eyeBlinksEpochs, ~, dataeog, ~, threshold] = detect_veog(DATA, blink_duration, blink_iqr, cfg.figure.visible);
%
% if ~isempty(eyeBlinksEpochs)
%     multiBlink = detect_multiblinks(eyeBlinksEpochs, 0);
%     fprintf('Number of detected blinks: %d.\n', size(eyeBlinksEpochs, 1));
%     fprintf('Number of detected multiple blinks within each evaluation window: %d.\n', sum(multiBlink));
%
%     eyeBlinksEpochs(multiBlink, :) = [];
%     NTRL_2 = size(eyeBlinksEpochs, 1);
% else
%     NTRL_2 = 0;
% end
%
% if NTRL_2 > num_trial_min
%     L = mode(diff(eyeBlinksEpochs')) + 1;
%     timeBlink0 = (0:L-1) ./ DATA.srate * 1000;
%     timeBlink1_2 = timeBlink0 - blink_duration;
%
%     dataeegepoched_2 = NaN(NCHANEEG, L, NTRL_2);
%     dataeogepoched_2 = NaN(L, NTRL_2);
%     for i = 1:NTRL_2
%         dataeegepoched_2(:, :, i) = dataeeg(:, eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
%         dataeogepoched_2(:, i)   = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
%     end
%
%     % Define baseline and active windows
%     [~, col_500ms]  = min(abs(timeBlink0 - 500));
%     [~, col_1500ms] = min(abs(timeBlink0 - 1500));
%     [~, col_2500ms] = min(abs(timeBlink0 - 2500));
%     [~, col_3500ms] = min(abs(timeBlink0 - 3500));
%     col_4000ms = L;
%
%     % Baseline correct data using standard windows
%     dataeegepoched_2 = dataeegepoched_2 - mean(dataeegepoched_2(:, [1:col_500ms, col_3500ms:col_4000ms], :), 2);
%
%     % Average across trials to get the ERP
%     mean_blink_ERP = mean(dataeegepoched_2, 3);
%
%     % Calculate Peak-to-Peak (PtP) in the active window across all channels
%     active_window = mean_blink_ERP(:, col_1500ms:col_2500ms);
%     BlinkPtP = max(active_window, [], 2) - min(active_window, [], 2);
%
%     % Extract the top 5 values along with their local anterior indices
%     [worst_values, worst_local_indices] = maxk(BlinkPtP(maskAnterior), 5);
%
%     % Map the local pool indices back to the absolute 128-channel indices
%     ptp_worst_global_idx = anterior_global_indices(worst_local_indices);
%
%     % Compute the target mean
%     meanFrontalPtP = mean(worst_values);
%
%     % Optional console verification to see exactly which electrodes are being traced
%     worst_labels = {DATA.chanlocs(ptp_worst_global_idx).labels};
%     fprintf('Top 5 worst channels selected for audit: %s\n', strjoin(worst_labels, ', '));
%
%     % ERP-EOG correlation
%     mean_blink_ERP = mean_blink_ERP(ptp_worst_global_idx, :);
%     mean_blink_ERP_avg = mean(mean_blink_ERP, 1)';
%     dataeogepoched_2_avg = mean(dataeogepoched_2, 2);
%     FrontalCorr_erp = abs(corr(mean_blink_ERP_avg(col_1500ms:col_2500ms), dataeogepoched_2_avg(col_1500ms:col_2500ms), "Type", "Spearman"));
%
% else
%     fprintf('Warning: No data (N = %d) to make an estimate of blink leftovers (Long Window).\n', NTRL_2);
%     BlinkPtP = NaN;
%     meanFrontalPtP = NaN;
%     FrontalCorr_erp = NaN;
% end

% =========================================================================
% Log
% =========================================================================
% fprintf('\nAverage Frontal Peak-to-Peak Leftover: %1.2f uV\n', meanFrontalPtP);
fprintf('Mean Frontal Amplitude Ratio: %.2f (Ideal is ~1.0)\n', meanFrontalRatio);
fprintf('Frontal ERP-VEOG Correlation: %1.2f\n', FrontalCorr_erp);
fprintf('Average Frontal EEG-VEOG Correlation: %1.2f\n', meanFrontalCorr_cont);

% fprintf(DATA.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
% fprintf(DATA.ALSUTRECHT.subject.fid,'Leftovers: eye blink artefacts\n');
% fprintf(DATA.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
% fprintf(DATA.ALSUTRECHT.subject.fid,'Average Frontal Peak-to-Peak Leftover: %1.2f uV\n', meanFrontalPtP);
% fprintf(DATA.ALSUTRECHT.subject.fid,'Average Frontal VEOG-EEG Correlation: r = %1.2f\n', meanFrontalCorr_cont);

DATA.ALSUTRECHT.leftovers.blink1.blinksPtP      = MeanBlinkRatio;
DATA.ALSUTRECHT.leftovers.blink1.blinksTstat    = stats.tstat;
DATA.ALSUTRECHT.leftovers.blink1.corrContinuous = corr_continuous;
DATA.ALSUTRECHT.leftovers.blink1.corrERP        = FrontalCorr_erp;

% =========================================================================
% Plot
% =========================================================================
close all hidden;
fh = figure('Visible', cfg.figure.visible);
tiledlayout(2, 3, "TileSpacing", "compact", "Padding", "compact");

chanlocs = readlocs('biosemi128_eeglab.ced');
myCmap2 = brewermap(128, '*RdBu');
myCmap3 = brewermap(128, 'Reds');

if NTRL_1 > num_trial_min
    % Short Window EOG Traces
    nexttile(1); hold on;
    plot(timeBlink1_1, 0 * ones(size(timeBlink1_1)), 'Color', 0 * [1 1 1], 'LineStyle', '-');
    plot(timeBlink1_1, dataeogepoched_1, 'LineWidth', 1.2);
    plot(timeBlink1_1, threshold * ones(size(timeBlink1_1)), 'Color', 0.5 * [1 1 1], 'LineStyle', '--');

    myCmap1 = brewermap(NTRL_1, 'YlGn');
    set(gca, 'ColorOrder', [0 0 0; myCmap1; 0.5 * [1 1 1]]);
    title(['Detected blinks, N = ' num2str(NTRL_1)]);
    pbaspect([1.618 1 1]); ylabel('EOG amplitude (\muV)');
    xlim([min(timeBlink1_1), max(timeBlink1_1)]);

    % Topoplot of T-statistics
    nexttile(4); hold on;
    maxBlinkTstat = prctile(abs(stats.tstat), 95);
    maxBlinkTstat = max(maxBlinkTstat, 2.5); % Ensure colour bar is scaled nicely
    topoplot(stats.tstat, chanlocs, 'maplimits', maxBlinkTstat * [-1 1], 'headrad', 0.5, 'colormap', myCmap2, 'whitebk', 'on', 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    title({'T-stats EEG timelocked to blinks', ['Frontal t-stat = ' num2str(round(meanFrontalTstat, 1))]});
    hcb = colorbar;
    hcb.Title.String = "T-value";
end

if NTRL_2 > num_trial_min
    % Long Window EOG Traces
    nexttile(2); hold on;
    plot(timeBlink1_2, 0 * ones(size(timeBlink1_2)), 'Color', 0 * [1 1 1], 'LineStyle', '-');
    plot(timeBlink1_2, dataeogepoched_2, 'LineWidth', 1.2);
    plot(timeBlink1_2, threshold * ones(size(timeBlink1_2)), 'Color', 0.5 * [1 1 1], 'LineStyle', '--');

    myCmap1 = brewermap(NTRL_2, 'YlGn');
    set(gca, 'ColorOrder', [0 0 0; myCmap1; 0.5 * [1 1 1]]);
    title(['Detected blinks, N = ' num2str(NTRL_2)]);
    pbaspect([1.618 1 1]); ylabel('EOG (\muV)');
    xlim([min(timeBlink1_2), max(timeBlink1_2)]);

    % Topoplot of Peak-to-Peak Amplitude
    nexttile(5);
    % Robust upper limit to ignore extreme single-channel outliers
    maxMeanBlinkRatio = prctile(MeanBlinkRatio, 95);
    % Ensure the maximum is at least 1.1 so maplimits [min max] has a valid range
    maxMeanBlinkRatio = max(maxMeanBlinkRatio, 1.1);
    topoplot(MeanBlinkRatio, chanlocs, 'maplimits', [0.9 maxMeanBlinkRatio], 'headrad', 0.5, 'colormap', myCmap3, 'whitebk', 'on', 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
    title({'Peak-to-Baseline ratio', ['Frontal P2B = ' num2str(round(meanFrontalRatio, 2))]});
    hcb = colorbar;
    hcb.Title.String = "Ratio";

    %  Frontal ERP Trace
    nexttile(3); hold on;
    plot(timeBlink1_2, mean_blink_ERP, 'Color', [0.8 0.2 0.2 0.5], 'LineWidth', 1);
    plot(timeBlink1_2, mean_blink_ERP_avg, 'k', 'LineWidth', 2);
    % Plot the active window boundaries
    xline(win_active(1), '--k');
    xline(win_active(2), '--k');
    % Plot the baseline boundaries (optional)
    xline(win_baseline_pre(1), ':k');
    xline(win_baseline_pre(2), ':k');
    xline(win_baseline_post(1), ':k');
    xline(win_baseline_post(2), ':k');
    title({'Frontal Channels ERP', ['ERP-VEOG R = ' num2str(round(FrontalCorr_erp, 1))]});
    ylabel('Amplitude (\muV)'); xlabel('Time (ms)');
    pbaspect([1.618 1 1]);
    xlim([min(timeBlink1_2), max(timeBlink1_2)]);
end

% EEG-EOG correlation
nexttile(6);
maxBlinkCorr = prctile(corr_continuous, 95);
maxBlinkCorr = max(maxBlinkCorr, 0.10); % Ensure colour bar is scaled nicely
topoplot(corr_continuous, chanlocs, 'maplimits', [0 maxBlinkCorr], 'headrad', 0.5, 'colormap', myCmap3, 'whitebk', 'on', 'electrodes', 'off', 'style', 'map', 'shading', 'interp');
title({'EEG-VEOG correlation', ['Frontal R = ' num2str(round(meanFrontalCorr_cont, 2))]});
hcb = colorbar;
hcb.Title.String = "abs(R)";

% Save figure
save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_leftovers_' num2str(tag_figure)], [30 20]);

% =========================================================================
% Compare Leftovers to ICA Templates & Generate Final Report
% =========================================================================
% Evaluation Thresholds
blinkThreshold2_Tstat = 4.0;  % Statistical significance threshold
blinkThreshold_Ratio  = 2.0;  % Peak-to-Baseline ratio threshold (Ideal is ~1.0)
blinkThreshold_corr1  = 0.2;  % Continuous correlation threshold
blinkThreshold_corr2  = 0.2;  % ERP correlation threshold

% Initialise flags
flag_MeanTstat = false;
flag_MeanRatio = false;
flag_ERPCorr   = false;

if ~isnan(meanFrontalTstat) && NTRL_1 > num_trial_min
    flag_MeanTstat = meanFrontalTstat > blinkThreshold2_Tstat;
end

if ~isnan(meanFrontalRatio) && NTRL_2 > num_trial_min
    flag_MeanRatio = meanFrontalRatio > blinkThreshold_Ratio;
end

if ~isnan(FrontalCorr_erp) && NTRL_2 > num_trial_min
    flag_ERPCorr = FrontalCorr_erp > blinkThreshold_corr2;
end

flag_ContCorr = meanFrontalCorr_cont > blinkThreshold_corr1;

% Final REDO Logic
flag_redo = flag_MeanTstat | flag_MeanRatio | flag_ContCorr | flag_ERPCorr;

% Save flags to DATA structure
DATA.ALSUTRECHT.leftovers.blink1.flag_redo     = flag_redo;
DATA.ALSUTRECHT.leftovers.blink1.flagMeanTstat = flag_MeanTstat;
DATA.ALSUTRECHT.leftovers.blink1.flagMeanRatio = flag_MeanRatio;
DATA.ALSUTRECHT.leftovers.blink1.flag_ContCorr = flag_ContCorr;
DATA.ALSUTRECHT.leftovers.blink1.flag_ERPCorr  = flag_ERPCorr;

% =========================================================================
% Console and Log File Reporting
% =========================================================================
fprintf('\n=========================================================\n');
fprintf('FINAL BLINK LEFTOVER QUALITY ASSURANCE REPORT\n');
fprintf('=========================================================\n');

if NTRL_1 <= num_trial_min && NTRL_2 <= num_trial_min
    fprintf('STATUS: [WARNING] - Insufficient blinks detected to run evaluation.\n');
elseif flag_redo
    fprintf('STATUS: [FAIL] - Significant blink leakage detected. Fixing recommended.\n');
else
    fprintf('STATUS: [PASS] - Blink cleaning is within acceptable physiological limits.\n');
end

fprintf('\n--- Leftovers ---\n');
fprintf('%s Average T-Statistic      : %1.2f (Threshold: %1.2f)\n', get_status_str(flag_MeanTstat), meanFrontalTstat, blinkThreshold2_Tstat);
fprintf('%s Average Peak-to-Baseline : %1.2f (Threshold: %1.2f uV)\n', get_status_str(flag_MeanRatio), meanFrontalRatio, blinkThreshold_Ratio);
fprintf('%s Continuous Correlation   : %1.2f (Threshold: %1.2f)\n', get_status_str(flag_ContCorr), meanFrontalCorr_cont, blinkThreshold_corr1);
fprintf('%s ERP Correlation          : %1.2f (Threshold: %1.2f)\n', get_status_str(flag_ERPCorr), FrontalCorr_erp, blinkThreshold_corr2);
fprintf('=========================================================\n\n');

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