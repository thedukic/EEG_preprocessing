function EEG = detect_ic_bad(EEG, EXT, EMG, cfg)
% DETECT_IC_BAD Master function to categorise and flag artifactual ICs.

% Define thresholds
cfg.thresh.ext          = 3;
cfg.thresh.temp_blink   = 0.95;
cfg.thresh.temp_saccade = 0.95;
cfg.thresh.temp_heart   = 0.95;
cfg.thresh.ctps_pk      = 20;
cfg.thresh.emg_r        = 5;
cfg.thresh.channel      = 20;

fprintf('\n================================\n');
fprintf('Detecting bad ICs\n');
fprintf('================================\n');

% Double-check ICA presence
EEG = eeg_checkset(EEG, 'ica');
EEG.icaact = [];

% Precompute shared ICA variables
data_ica      = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind, :);
icawinv       = EEG.icawinv;
% icawinvSmooth = estimate_invlaplacian(icawinv, EEG.chanlocs, 1);
templatesICA  = load_ictemplateweights(EEG);
num_ics       = size(data_ica, 1);

%% 1. Global ICLabel Execution
fprintf('\n--------------------------------\n');
fprintf('ICLabel\n');
fprintf('--------------------------------\n');

EEG = iclabel(EEG);
EEG = pop_icflag(EEG, cfg.ica.iclabel);

EEG.ALSUTRECHT.ica.ICLabel.bics = EEG.reject.gcompreject;
EEG.ALSUTRECHT.ica.ICLabel.clss = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel.pvec, EEG.ALSUTRECHT.ica.ICLabel.cvec] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2);

%% 2. Accumulate Evidence via Subfunctions
EEG = detect_ic_eye(EEG, EXT, cfg, data_ica, icawinv, templatesICA, num_ics);
EEG = detect_ic_heart(EEG, EXT, cfg, data_ica, icawinv, templatesICA, num_ics);
EEG = detect_ic_muscle(EEG, cfg, data_ica, num_ics);
EEG = detect_ic_channel(EEG, cfg, num_ics);

%% 3. Final Combination Logic
fprintf('\n--------------------------------\n');
fprintf('Combining detected ICs\n');
fprintf('--------------------------------\n');

ICsMostLikelyEye     = EEG.ALSUTRECHT.ica.eye.bics;
ICsMostLikelyMuscle  = EEG.ALSUTRECHT.ica.muscle.bics;
ICsMostLikelyComplex = ICsMostLikelyEye & ICsMostLikelyMuscle;
ICsMostLikelyEye(ICsMostLikelyComplex) = false;
ICsMostLikelyMuscle(ICsMostLikelyComplex) = false;

ICsMostLikelyChannel = EEG.ALSUTRECHT.ica.channel.bics;
ICsMostLikelyChannelWrong = ICsMostLikelyChannel & (ICsMostLikelyEye | ICsMostLikelyMuscle | ICsMostLikelyComplex);
ICsMostLikelyChannel(ICsMostLikelyChannelWrong) = false;

ICsMostLikelyHeart = EEG.ALSUTRECHT.ica.heart.bics;
ICsMostLikelyHeart(ICsMostLikelyEye) = false;

% Validate no overlap
assert(max(sum([ICsMostLikelyEye, ICsMostLikelyMuscle, ICsMostLikelyComplex, ICsMostLikelyChannel, ICsMostLikelyHeart], 2)) == 1);

% Final Logging
EEG.ALSUTRECHT.ica.final.eye     = ICsMostLikelyEye;
EEG.ALSUTRECHT.ica.final.muscle  = ICsMostLikelyMuscle;
EEG.ALSUTRECHT.ica.final.complex = ICsMostLikelyComplex;
EEG.ALSUTRECHT.ica.final.channel = ICsMostLikelyChannel;
EEG.ALSUTRECHT.ica.final.heart   = ICsMostLikelyHeart;

EEG.ALSUTRECHT.ica.final.report = EEG.ALSUTRECHT.ica.ICLabel.cvec;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyMuscle)  = 2;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyEye)     = 3;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyHeart)   = 4;
EEG.ALSUTRECHT.ica.final.report(ICsMostLikelyChannel) = 6;

end

% =========================================================================
% SUBFUNCTION 1: EYE ARTIFACTS
% =========================================================================
function EEG = detect_ic_eye(EEG, EXT, cfg, data_ica, icawinv, templatesICA, num_ics)
% -------------------------------------------------------------------------
% Correlatioin with the external electrodes
% -------------------------------------------------------------------------
fprintf('\n--------------------------------\n');
fprintf('Eye artifacts\n');
fprintf('--------------------------------\n');

channel_veog = find(strcmp({EXT.chanlocs.labels}, 'VEOG'));
channel_heog = find(strcmp({EXT.chanlocs.labels}, 'HEOG'));
data_ext_eog = EXT.data([channel_veog channel_heog], :);

% Filter EOG (0.3 - 10 Hz)
[bh_eog, ah_eog] = butter(2, 0.3/(EEG.srate/2), 'high');
[bl_eog, al_eog] = butter(2, 10/(EEG.srate/2), 'low');

data_ica_eog = do_filteringcore(bl_eog, al_eog, data_ica, EEG.event, EEG.srate);
data_ica_eog = do_filteringcore(bh_eog, ah_eog, data_ica_eog, EEG.event, EEG.srate)';

data_ext_eog = do_filteringcore(bl_eog, al_eog, data_ext_eog, EEG.event, EEG.srate);
data_ext_eog = do_filteringcore(bh_eog, ah_eog, data_ext_eog, EEG.event, EEG.srate)';

% Correlations
corr_veog = corr(data_ica_eog, data_ext_eog(:, 1), "type", "Spearman");
corr_heog = corr(data_ica_eog, data_ext_eog(:, 2), "type", "Spearman");

corr_ext_eye = abs(zscore([corr_veog, corr_heog]));

% Evaluate VEOG (allow multiple)
bad_ic_veog = find(corr_ext_eye(:, 1) > cfg.thresh.ext);

% Evaluate HEOG (only max 1)
bad_ic_heog = find(corr_ext_eye(:, 2) > cfg.thresh.ext);
if length(bad_ic_heog) > 1
    [~, mostlikely] = max(corr_ext_eye(bad_ic_heog, 2));
    bad_ic_heog = bad_ic_heog(mostlikely);
end

% Log Eye Correlations
EEG.ALSUTRECHT.ica.eye.corr_eye.corr = single(corr_ext_eye);
EEG.ALSUTRECHT.ica.eye.corr_eye.bics = [bad_ic_veog(:); bad_ic_heog(:)];
EEG.ALSUTRECHT.ica.eye.corr_eye.cvec = [ones(length(bad_ic_veog), 1); 2*ones(length(bad_ic_heog), 1)];
EEG.ALSUTRECHT.ica.eye.corr_eye.clss = {'VEOG', 'HEOG'};

% -------------------------------------------------------------------------
% icablinkmetrics plugin
% -------------------------------------------------------------------------
 data_eog_blink = mean(EEG.data(ismember({EEG.chanlocs.labels}, cfg.ica.blinkchans), :), 1);
% data_eog_blink = mean(EEG.data(ismember({EXT.chanlocs.labels}, 'VEOG'), :), 1);

EEG.icaact = data_ica;

try
    icablinkmetricsout = icablinkmetrics(EEG, 'ArtifactChannel', data_eog_blink, 'Alpha', 0.001, 'VisualizeData', 'False');
    if any(icablinkmetricsout.identifiedcomponents > 0)
        fprintf('Blink ICs (N = %d) identified using icablinkmetrics.\n', length(icablinkmetricsout.identifiedcomponents));
    else
        icablinkmetricsout.identifiedcomponents = [];
    end
catch
    fprintf('icablinkmetrics failed. Skipping...\n');
    icablinkmetricsout.identifiedcomponents = [];
    icablinkmetricsout.metrics.corr_Pvalue  = [];
    icablinkmetricsout.metrics.conv_Pvalue  = [];
    icablinkmetricsout.metrics.perc_Pvalue  = [];
end
EEG.icaact = [];

% Log icablinkmetrics
if ~isempty(icablinkmetricsout.identifiedcomponents)
    % The method ran successfully
    EEG.ALSUTRECHT.ica.eye.icablinkmetrics.bics = false(num_ics, 1);
    EEG.ALSUTRECHT.ica.eye.icablinkmetrics.bics(icablinkmetricsout.identifiedcomponents(:)) = true;
    EEG.ALSUTRECHT.ica.eye.icablinkmetrics.pval = [icablinkmetricsout.metrics.corr_Pvalue; icablinkmetricsout.metrics.conv_Pvalue; icablinkmetricsout.metrics.perc_Pvalue]';
else
    EEG.ALSUTRECHT.ica.eye.icablinkmetrics.bics = NaN(num_ics, 1);
    EEG.ALSUTRECHT.ica.eye.icablinkmetrics.pval = [];
end

% -------------------------------------------------------------------------
% Spatial Templates
% -------------------------------------------------------------------------
% Blink
corrMatB1 = abs(corr(icawinv, templatesICA.Blinkweights0, "type", "Spearman"));
corrMatB2 = abs(corr(icawinv, templatesICA.Blinkweights1, "type", "Spearman"));
badIC_blink = unique([find(corrMatB1 > cfg.thresh.temp_blink); find(corrMatB2 > cfg.thresh.temp_blink)]);

EEG.ALSUTRECHT.ica.eye.blinkTemplateCorr.pval1 = corrMatB1;
EEG.ALSUTRECHT.ica.eye.blinkTemplateCorr.pval2 = corrMatB2;
% EEG.ALSUTRECHT.ica.eye.blinkTemplateCorr.bics  = badIC_blink(:);
EEG.ALSUTRECHT.ica.eye.blinkTemplateCorr.bics = false(num_ics ,1);
EEG.ALSUTRECHT.ica.eye.blinkTemplateCorr.bics(badIC_blink) = true;

% Saccade
corrMatS1 = abs(corr(icawinv, templatesICA.Saccadeweights0,  "type", "Spearman"));
corrMatS2 = abs(corr(icawinv, templatesICA.Saccadeweights2L, "type", "Spearman"));
corrMatS3 = abs(corr(icawinv, templatesICA.Saccadeweights2R, "type", "Spearman"));
badIC_saccade = unique([find(corrMatS1 > cfg.thresh.temp_saccade); find(corrMatS2 > cfg.thresh.temp_saccade); find(corrMatS3 > cfg.thresh.temp_saccade)]);

EEG.ALSUTRECHT.ica.eye.saccadeTemplateCorr.pval1 = corrMatS1;
EEG.ALSUTRECHT.ica.eye.saccadeTemplateCorr.pval2 = corrMatS2;
EEG.ALSUTRECHT.ica.eye.saccadeTemplateCorr.pval3 = corrMatS3;

EEG.ALSUTRECHT.ica.eye.saccadeTemplateCorr.bics = false(num_ics ,1);
EEG.ALSUTRECHT.ica.eye.saccadeTemplateCorr.bics(badIC_saccade) = true;

% -------------------------------------------------------------------------
% EyeCatch
% -------------------------------------------------------------------------
eyeDetector = eyeCatch;     % create an object from the class. Once you made an object it can
% be used for multiple detections (much faster than creating an
% object each time).

[eyeIC, similarity, scalpmapObj] = eyeDetector.detectFromEEG(EEG);
% eyeIC                          % display the IC numbers for eye ICs
% scalpmapObj.plot(eyeIC)        % plot eye ICs

EEG.ALSUTRECHT.ica.eye.EyeCatch.bics = eyeIC(:);
EEG.ALSUTRECHT.ica.eye.EyeCatch.pval = similarity(:);

% -------------------------------------------------------------------------
% Combine
% -------------------------------------------------------------------------
% *Eye ICs (ICLabel)
% ICsMostLikelyEyeICLabel = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 3);
ICsMostLikelyEyeICLabel = EEG.ALSUTRECHT.ica.ICLabel.bics & EEG.ALSUTRECHT.ica.ICLabel.cvec == 3;

% *Eye ICs (EyeCatch)
ICsMostLikelyEyeEyeCatch = EEG.ALSUTRECHT.ica.eye.EyeCatch.bics;

% *Blink ICs
m = EEG.ALSUTRECHT.ica.eye.corr_eye.bics(EEG.ALSUTRECHT.ica.eye.corr_eye.cvec == 1); % VEOG
blink1 = false(num_ics, 1);
blink1(m) = true;
blink2 = EEG.ALSUTRECHT.ica.eye.icablinkmetrics.bics;
blink3 = EEG.ALSUTRECHT.ica.eye.blinkTemplateCorr.bics;
ICsMostLikelyBlink = [blink1(:), blink2(:), blink3(:)];
[ICsMostLikelyBlink, confidence_score] = determine_likelihood(ICsMostLikelyBlink);

% *Saccades ICs
m = EEG.ALSUTRECHT.ica.eye.corr_eye.bics(EEG.ALSUTRECHT.ica.eye.corr_eye.cvec == 2); % HEOG
saccade1 = false(num_ics, 1);
saccade1(m) = true;
saccade2 = EEG.ALSUTRECHT.ica.eye.saccadeTemplateCorr.bics;
ICsMostLikelySaccade = [saccade1(:), saccade2(:)];
[ICsMostLikelySaccade, confidence_score] = determine_likelihood(ICsMostLikelySaccade);

% Construct the Eye Evidence Matrix
% Stack the results from the various eye detection methods side-by-side
% Using double preserves any NaNs from failed methods
ICsMostLikelyEyeMatrix = [double(ICsMostLikelyBlink(:)), ...
    double(ICsMostLikelySaccade(:)), ...
    double(ICsMostLikelyEyeICLabel(:)), ...
    double(ICsMostLikelyEyeEyeCatch(:))];

% Strip Columns That Failed Entirely
% If a method completely failed to execute and returned nothing but NaNs,
% we drop the entire column so it does not interfere with the matrix evaluation
failed_method_columns = all(isnan(ICsMostLikelyEyeMatrix), 1);
ICsMostLikelyEyeMatrix(:, failed_method_columns) = [];

% Run the "Any" Column Evaluation
% If an independent component is flagged as an eye artifact by ANY of your active,
% successfully executed pipelines, we mark it as true
ICsMostLikelyEye = any(ICsMostLikelyEyeMatrix == 1, 2);

% Log the evidence matrix
EEG.ALSUTRECHT.ica.eye.evidence_matrix = [];
EEG.ALSUTRECHT.ica.eye.confidence      = [];
EEG.ALSUTRECHT.ica.eye.bics            = ICsMostLikelyEye;

end

% =========================================================================
% SUBFUNCTION 2: HEART ARTIFACTS
% =========================================================================
function EEG = detect_ic_heart(EEG, EXT, cfg, data_ica, icawinv, templatesICA, num_ics)

% -------------------------------------------------------------------------
fprintf('\n--------------------------------\n');
fprintf('Heart artifacts\n');
fprintf('--------------------------------\n');

% -------------------------------------------------------------------------
% ICLabel
% -------------------------------------------------------------------------
% bad_ic_iclabel = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 4);
bad_ic_iclabel = EEG.ALSUTRECHT.ica.ICLabel.bics & EEG.ALSUTRECHT.ica.ICLabel.cvec == 4;
bad_ic_iclabel = double(bad_ic_iclabel);

% -------------------------------------------------------------------------
% External
% -------------------------------------------------------------------------
channel_ecg = find(strcmp({EXT.chanlocs.labels}, 'ECG'));

if isempty(channel_ecg)
    warning('ECG signal not recorded. Heart IC correlations skipped.');
    ecg_type = 'None';
    pulse_estimate = NaN;
    % Must be NaN vector, not [], to maintain matrix structure for the NaN mask
    bad_ic_ext   = NaN(num_ics, 1);
    bad_ic_ctps  = NaN(num_ics, 1);
    corr_ext_ecg = NaN(num_ics, 1);
    V  = NaN(num_ics, 1);
    pK = NaN(num_ics, 1);
else
    data_ext_ecg = EXT.data(channel_ecg, :);
    ecg_type = 'ECG';

    % Filter ECG (10 - 20 Hz)
    [bh_ecg, ah_ecg] = butter(2, 10/(EEG.srate/2), 'high');
    [bl_ecg, al_ecg] = butter(2, 20/(EEG.srate/2), 'low');

    data_ica_ecg = do_filteringcore(bl_ecg, al_ecg, data_ica, EEG.event, EEG.srate);
    data_ica_ecg = do_filteringcore(bh_ecg, ah_ecg, data_ica_ecg, EEG.event, EEG.srate)';

    data_ext_ecg = do_filteringcore(bl_ecg, al_ecg, data_ext_ecg, EEG.event, EEG.srate);
    data_ext_ecg = do_filteringcore(bh_ecg, ah_ecg, data_ext_ecg, EEG.event, EEG.srate)';

    % Correlation
    corr_ecg     = corr(data_ica_ecg.^2, abs(data_ext_ecg), "type", "Spearman");
    corr_ext_ecg = abs(zscore(corr_ecg));
    bad_ic_ext   = double(corr_ext_ecg > cfg.thresh.ext);

    % CTPS
    [ecg_mask, ecg_epoch, ~, ~, ~, pulse_estimate] = detect_ecg(EXT, [-200 250], ecg_type, cfg.figure.visible);
    if ~isnan(ecg_mask)
        [V, pK] = my_ctps(data_ica, ecg_epoch, EEG.event, EEG.srate);
        bad_ic_ctps = double(pK >= cfg.thresh.ctps_pk);
    else
        V  = NaN(num_ics, 1);
        pK = NaN(num_ics, 1);
        bad_ic_ctps = NaN(num_ics, 1);
    end
end

% Log ECG Correlation
EEG.ALSUTRECHT.ica.heart.corr_heart.corr = corr_ext_ecg;
EEG.ALSUTRECHT.ica.heart.corr_heart.bics = bad_ic_ext;

% Log CTPS
EEG.ALSUTRECHT.ica.heart.ctps.v     = V;
EEG.ALSUTRECHT.ica.heart.ctps.pk    = pK;
EEG.ALSUTRECHT.ica.heart.ctps.bics  = bad_ic_ctps;
EEG.ALSUTRECHT.subject.pulsEstimate = pulse_estimate;

% -------------------------------------------------------------------------
% Templates
% -------------------------------------------------------------------------
corrMatH1 = abs(corr(icawinv, templatesICA.Heartweights0, "type", "Spearman"));
bad_ic_temp_1 = double(any(corrMatH1 > cfg.thresh.temp_heart, 2));

if ~isnan(templatesICA.Heartweights1)
    corrMatH2 = abs(corr(icawinv, templatesICA.Heartweights1, "type", "Spearman"));
    bad_ic_temp_2 = double(corrMatH2 > cfg.thresh.temp_heart);
else
    corrMatH2 = NaN(num_ics, 1);
    bad_ic_temp_2 = NaN(num_ics, 1);
end

if ~isnan(templatesICA.Heartweights2)
    corrMatH3 = abs(corr(icawinv, templatesICA.Heartweights2, "type", "Spearman"));
    bad_ic_temp_3 = double(corrMatH3 > cfg.thresh.temp_heart);
else
    corrMatH3 = NaN(num_ics, 1);
    bad_ic_temp_3 = NaN(num_ics, 1);
end

bad_ic_temp = [bad_ic_temp_1, bad_ic_temp_2, bad_ic_temp_3];

EEG.ALSUTRECHT.ica.heart.heartTemplateCorr.pval1 = corrMatH1;
EEG.ALSUTRECHT.ica.heart.heartTemplateCorr.pval2 = corrMatH2;
EEG.ALSUTRECHT.ica.heart.heartTemplateCorr.pval3 = corrMatH3;
EEG.ALSUTRECHT.ica.heart.heartTemplateCorr.bics  = bad_ic_temp;

% -------------------------------------------------------------------------
% Combine
% -------------------------------------------------------------------------
evidence_matrix_heart = [bad_ic_iclabel, bad_ic_ext, bad_ic_ctps, bad_ic_temp];
[final_bad_heart_ics, confidence_score] = determine_likelihood(evidence_matrix_heart);

% Log the evidence matrix
EEG.ALSUTRECHT.ica.heart.evidence_matrix = evidence_matrix_heart;
EEG.ALSUTRECHT.ica.heart.confidence      = confidence_score;
EEG.ALSUTRECHT.ica.heart.bics            = final_bad_heart_ics;

end

% =========================================================================
% SUBFUNCTION 3: MUSCLE ARTIFACTS
% =========================================================================
function EEG = detect_ic_muscle(EEG, cfg, data_ica, num_ics)

fprintf('\n--------------------------------\n');
fprintf('Muscle artifacts\n');
fprintf('--------------------------------\n');

% -------------------------------------------------------------------------
% 1. ICLabel
% -------------------------------------------------------------------------
% Class 2 corresponds to Muscle/EMG in ICLabel
bad_ic_iclabel = EEG.ALSUTRECHT.ica.ICLabel.bics & EEG.ALSUTRECHT.ica.ICLabel.cvec == 2;
bad_ic_iclabel = double(bad_ic_iclabel);

% -------------------------------------------------------------------------
% 2. Power Slopes (Broad and High)
% -------------------------------------------------------------------------
options.Freq_to_compute       = [1 100];
options.muscleFreqEx          = 50 + 2*[-1 1];
options.muscleFreq1           = cfg.ica.emgSlopeFreq1;
options.muscleFreq2           = cfg.ica.emgSlopeFreq2;
options.muscleSlopeThreshold1 = cfg.ica.emgSlopeThreshold1;
options.muscleSlopeThreshold2 = cfg.ica.emgSlopeThreshold2;

[pow, frefAll] = pwelch(data_ica', size(data_ica, 2), [], size(data_ica, 2), EEG.srate);
pow = pow';
frefAll = frefAll';

freq = options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
fftBins = zeros(size(pow, 1), size(freq, 2));
for i_freq = 1:length(freq)
    [~, index1] = min(abs(frefAll-((freq(1,i_freq)-0.25))));
    [~, index2] = min(abs(frefAll-((freq(1,i_freq)+0.25))));
    fftBins(:, i_freq) = mean(pow(:, index1:index2), 2);
end

slope_muscle = NaN(num_ics, 2);
for i_ic = 1:num_ics
    % Broad Window
    if ~isempty(options.muscleFreq1)
        [~, fin1] = min(abs(options.muscleFreq1(1) - freq));
        [~, fin2] = min(abs(options.muscleFreq1(2) - freq));
        freqHz_broad = freq(fin1:fin2);
        freqPow_broad = fftBins(i_ic, fin1:fin2);
    else
        freqHz_broad = freq;
        freqPow_broad = fftBins(i_ic, :);
    end
    if ~isempty(options.muscleFreqEx)
        [~, fex1] = min(abs(options.muscleFreqEx(1) - freqHz_broad));
        [~, fex2] = min(abs(options.muscleFreqEx(2) - freqHz_broad));
        if fex1 <= length(freqHz_broad)
            freqHz_broad(fex1:fex2) = [];
            freqPow_broad(fex1:fex2) = [];
        end
    end
    p_broad = polyfit(log10(freqHz_broad), log10(freqPow_broad), 1);

    % High Window
    [~, fhi1] = min(abs(options.muscleFreq2(1) - freq));
    [~, fhi2] = min(abs(options.muscleFreq2(2) - freq));
    freqHz_high = freq(fhi1:fhi2);
    freqPow_high = fftBins(i_ic, fhi1:fhi2);
    if ~isempty(options.muscleFreqEx)
        [~, fex1] = min(abs(options.muscleFreqEx(1) - freqHz_high));
        [~, fex2] = min(abs(options.muscleFreqEx(2) - freqHz_high));
        if fex1 <= length(freqHz_high)
            freqHz_high(fex1:fex2) = [];
            freqPow_high(fex1:fex2) = [];
        end
    end
    p_high = polyfit(log10(freqHz_high), log10(freqPow_high), 1);

    slope_muscle(i_ic, :) = [p_broad(1), p_high(1)];
end

bad_ic_slope_broad = double(slope_muscle(:, 1) > options.muscleSlopeThreshold1);
bad_ic_slope_high  = double(slope_muscle(:, 2) > options.muscleSlopeThreshold2);

% -------------------------------------------------------------------------
% 3. Spectral Metrics
% -------------------------------------------------------------------------
[neg_integral, slope_low, slope_discrepancy] = estimate_spectral_metrics(fftBins, freq, [1 70], [2 40]);
r = neg_integral ./ abs(slope_low);
bad_ic_ratio = double(r > cfg.thresh.emg_r);

% -------------------------------------------------------------------------
% Combine Evidence
% -------------------------------------------------------------------------
% Matrix Columns: [ICLabel, Broad Slope, High Slope, Spectral Ratio]
evidence_matrix_muscle = [bad_ic_iclabel(:), bad_ic_slope_broad(:), bad_ic_slope_high(:), bad_ic_ratio(:)];
[final_bad_muscle_ics, confidence_score] = determine_likelihood(evidence_matrix_muscle);

% -------------------------------------------------------------------------
% Log Results
% -------------------------------------------------------------------------
EEG.ALSUTRECHT.ica.muscle.slope             = slope_muscle;
EEG.ALSUTRECHT.ica.muscle.neg_integral      = neg_integral(:);
EEG.ALSUTRECHT.ica.muscle.slope_discrepancy = slope_discrepancy(:);

% Log the new matrix metrics
EEG.ALSUTRECHT.ica.muscle.evidence_matrix = evidence_matrix_muscle;
EEG.ALSUTRECHT.ica.muscle.confidence      = confidence_score;
EEG.ALSUTRECHT.ica.muscle.bics            = final_bad_muscle_ics;

end

% =========================================================================
% SUBFUNCTION 4: CHANNEL ARTIFACTS
% =========================================================================
function EEG = detect_ic_channel(EEG, cfg, num_ics)
fprintf('\n--------------------------------\n');
fprintf('Channel artifacts\n');
fprintf('--------------------------------\n');

[spatialSmoothness, bad_ic] = estimate_spatialsmoothnes(EEG, cfg.thresh.channel);

% -------------------------------------------------------------------------
% Combine Evidence
% -------------------------------------------------------------------------
% channel1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 6);
channel1 = EEG.ALSUTRECHT.ica.ICLabel.bics & EEG.ALSUTRECHT.ica.ICLabel.cvec == 6;
channel2 = false(num_ics, 1);
channel2(bad_ic) = true;

evidence_matrix_channel = [channel1, channel2];
[final_bad_channel_ics, confidence_score] = determine_likelihood(evidence_matrix_channel);

% Log Channel
EEG.ALSUTRECHT.ica.channel.spatialSmoothness.pval = spatialSmoothness;
EEG.ALSUTRECHT.ica.channel.spatialSmoothness.bics = bad_ic(:);

% Log the new matrix metrics
EEG.ALSUTRECHT.ica.channel.evidence_matrix = evidence_matrix_channel;
EEG.ALSUTRECHT.ica.channel.confidence      = confidence_score;
EEG.ALSUTRECHT.ica.channel.bics            = final_bad_channel_ics;

end

% ========================================================================
% SUBFUNCTION 5: LIKELY ARTIFACTS
% =========================================================================
function [final_bad_ics, confidence_score] = determine_likelihood(evidence_matrix)
% DETERMINE_LIKELIHOOD Consolidates artifact evidence across active pipelines
% using a dynamic, row-wise consensus threshold.

% Ensure double precision for accurate division
evidence_matrix = double(evidence_matrix);

% Calculate the sum of evidence for each IC, ignoring missing entries safely
evidence_sum = sum(evidence_matrix, 2, 'omitnan');

% Dynamically count exactly how many valid methods ran for EACH specific IC
% By summing the logical inversion of NaN, we get a true per-row count
num_methods_per_ic = sum(~isnan(evidence_matrix), 2);

% Handle the extreme safety case where an IC has no valid data at all
% This prevents a divide-by-zero warning resulting in NaN scores
num_methods_per_ic(num_methods_per_ic == 0) = 1;

% Compute the dynamic confidence score (ranging strictly from 0.0 to 1.0)
confidence_score = evidence_sum ./ num_methods_per_ic;

% Define the consensus acceptance rule (e.g., flagged by 50% or more of active methods)
threshold_probability = 0.50;

% Extract the absolute indices of the final rejected artifact ICs
final_bad_ics = confidence_score >= threshold_probability;

end

function [neg_integral, slope_low, slope_discrepancy] = estimate_spectral_metrics(powspctrm, freqs, interest_range, censor_range)
% ESTIMATE_SPECTRAL_METRICS Computes slope discrepancy and negative integral metrics.
%
% Inputs:
%   powspctrm      : Matrix of power spectra (num_channels x freqs)
%   freqs          : Vector of frequencies (1 x freqs or freqs x 1)
%   interest_range : 1x2 vector of frequencies to include (e.g., [1 70])
%   censor_range   : 1x2 vector of frequencies to exclude (e.g., [3 30])
%
% Outputs:
%   slope_discrepancy : Vector (num_channels x 1) of internal consistency differences
%   neg_integral      : Vector (num_channels x 1) of the total negative residual area

% Ensure column vector for frequencies
freqs = freqs(:);
num_channels = size(powspctrm, 1);

log_freqs = log10(freqs);
log_power_spectra = log10(powspctrm);

% --- FREQUENCY MASKING ---
interest_idx1 = freqs >= interest_range(1) & freqs <= interest_range(2);
interest_idx2 = freqs < censor_range(1)    | freqs > censor_range(2);
interest_idx3 = freqs < 49 | freqs > 51;
interest_idx4 = freqs < 99 | freqs > 101;
interest_idx  = interest_idx1 & interest_idx2 & interest_idx3 & interest_idx4;

freqs_interest = freqs(interest_idx);
log_freqs_select = log_freqs(interest_idx);
log_power_select = log_power_spectra(:, interest_idx);

% Identify the low/high frequency 'islands' within the filtered data
idx_low  = find(freqs_interest <= censor_range(1));
idx_high = find(freqs_interest >= censor_range(2));

% Fallback if data is entirely outside the censored range
if isempty(idx_low) && ~isempty(idx_high)
    idx_low  = find(freqs_interest <= freqs_interest(15));
    idx_high = find(freqs_interest >= freqs_interest(end-15));
end

% --- ALLOCATION & MODEL FITTING ---
slope_low = NaN(num_channels, 1);
slope_discrepancy = NaN(num_channels, 1);
ap_fit = NaN(num_channels, length(freqs));

for i_channel = 1:num_channels
    power_channel = log_power_select(i_channel, :);

    % Global fit to calculate the aperiodic baseline
    mdl = fitlm(log_freqs_select, power_channel(:));
    intercept = mdl.Coefficients.Estimate(1);
    slope     = mdl.Coefficients.Estimate(2);

    ap_fit_log = intercept + (slope * log_freqs');
    ap_fit(i_channel, :) = 10.^ap_fit_log;

    % Sub-segment fits for internal consistency
    mdl_low  = fitlm(log_freqs_select(idx_low), power_channel(idx_low));
    mdl_high = fitlm(log_freqs_select(idx_high), power_channel(idx_high));

    slope_low(i_channel, 1) = mdl_low.Coefficients.Estimate(2);
    slope_discrepancy(i_channel, 1) = abs(mdl_low.Coefficients.Estimate(2) - mdl_high.Coefficients.Estimate(2));
end

% --- NEGATIVE RESIDUAL ESTIMATION ---
periodic_estimate = log_power_spectra - log10(ap_fit);
periodic_estimate = periodic_estimate(:, interest_idx1);
neg_mask          = periodic_estimate < 0;

% Area of 'impossible' power per channel
neg_integral = sum(abs(periodic_estimate .* neg_mask), 2);

end