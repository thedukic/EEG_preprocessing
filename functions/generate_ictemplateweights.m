function generate_ictemplateweights(EEG, EMG, EXT, cfg)

fprintf('\n================================\n');
fprintf('Generating individualised ICA artifact templates\n');
fprintf('================================\n');

% -------------------------------------------------------------------------
% Merge blocks
DATA = merge_electrodetypes(EEG, EMG, EXT);

% -------------------------------------------------------------------------
% Eye blinks
fprintf('\n--------------------------------\n');
fprintf('Detecting eye blinks\n');
fprintf('--------------------------------\n');

trIQRblink = 3;
winBlink = 150;
[VEOGmask, VEOGepochs, VEOGmaxLatency, VEOGdata, VEOGEEGdata, treshold] = detect_veog(DATA, winBlink, trIQRblink, cfg.figure.visible);

% Eye saccades
fprintf('\n--------------------------------\n');
fprintf('Detecting eye saccades\n');
fprintf('--------------------------------\n');

trIQRsaccade = 1;
winSaccade = 150;
[HEOGmask, HEOGepochs, HEOGmaxLatency, HEOGdata, HEOGEEGdata, treshold] = detect_heog(DATA, winSaccade, trIQRsaccade, cfg.figure.visible);

% Heart beats
fprintf('\n--------------------------------\n');
fprintf('Detecting heartbeats\n');
fprintf('--------------------------------\n');

ECGtype = 'ecg';
winHeart = 30;
[ECGmask, ECGepochs, ECGlatency, ECGdata, ECGEEGdata, pulsEstimate] = detect_ecg(DATA, winHeart, cfg.figure.visible);

% -------------------------------------------------------------------------
% Setup availability flags to clean up downstream logic
% -------------------------------------------------------------------------
has_VEOG = any(VEOGmask);
if islogical(ECGmask)
    has_ECG  = any(ECGmask);
end

if iscell(HEOGepochs)
    assert(length(HEOGepochs) == 2);
    has_HEOG_L = any(HEOGmask{1});
    has_HEOG_R = any(HEOGmask{2});
    has_HEOG   = has_HEOG_L || has_HEOG_R;
else
    has_HEOG   = any(HEOGmask);
    has_HEOG_L = false;
    has_HEOG_R = false;
end

% -------------------------------------------------------------------------
% Estimate individualised templates & correlations
% -------------------------------------------------------------------------
fprintf('\n--------------------------------\n');
fprintf('Estimating individualised templates\n');
fprintf('--------------------------------\n');

% Pre-allocate outputs with zeros so topoplot does not crash on empty matrices
C_VEOG = NaN; A_VEOG = zeros(128,1); R_VEOGEEG = zeros(128,1);
C_ECG  = NaN; A_ECG  = zeros(128,1); R_ECGEEG  = zeros(128,1);
C_HEOG = NaN; A_HEOG = zeros(128,1); R_HEOGEEG = zeros(128,1);
C_HEOG_L = NaN; A_HEOG_L = zeros(128,1); R_HEOGEEG_L = zeros(128,1);
C_HEOG_R = NaN; A_HEOG_R = zeros(128,1); R_HEOGEEG_R = zeros(128,1);

fprintf('Estimating artifacts topoplots and covariance matrices...\n');

if has_VEOG
    [C_VEOG, A_VEOG] = estimate_params(DATA, VEOGmask);
    R_VEOGEEG = estimate_extcorrelation(DATA, VEOGmask, 'veog');
end

if has_ECG
    [C_ECG, A_ECG] = estimate_params(DATA, ECGmask);
    if strcmpi(ECGtype, 'ecg')
        R_ECGEEG = estimate_extcorrelation(DATA, ECGmask, 'ecg');
    end
end

if iscell(HEOGepochs)
    if has_HEOG_L
        [C_HEOG_L, A_HEOG_L] = estimate_params(DATA, HEOGmask{1});
        R_HEOGEEG_L = estimate_extcorrelation(DATA, HEOGmask{1}, 'heog');
    end
    if has_HEOG_R
        [C_HEOG_R, A_HEOG_R] = estimate_params(DATA, HEOGmask{2});
        R_HEOGEEG_R = estimate_extcorrelation(DATA, HEOGmask{2}, 'heog');
    end
else
    if has_HEOG
        [C_HEOG, A_HEOG] = estimate_params(DATA, HEOGmask);
        R_HEOGEEG = estimate_extcorrelation(DATA, HEOGmask, 'heog');
    end
end

% -------------------------------------------------------------------------
% Compare estimates with group average templates
% -------------------------------------------------------------------------
fprintf('Comparing the estimates with the group average templates...\n');

load('wBlink.mat',   'Blinkweights');
load('wSaccade.mat', 'Saccadeweights');
load('wHeart.mat',   'Heartweights');

% Pre-allocate correlations
R_VEOG1 = NaN; R_VEOG2 = NaN;
R_ECG1  = NaN; R_ECG2  = NaN;
R_HEOG1 = NaN; R_HEOG2 = NaN;
R_HEOG1_L = NaN; R_HEOG1_R = NaN;
R_HEOG2_L = NaN; R_HEOG2_R = NaN;

if has_VEOG
    [R_VEOG1, ~] = corr(Blinkweights, A_VEOG(1:128));
    [R_VEOG2, ~] = corr(Blinkweights, R_VEOGEEG);
    R_VEOG1 = round(R_VEOG1, 2);
    R_VEOG2 = round(R_VEOG2, 2);
end

if has_ECG
    [R_ECG1, ~] = corr(Heartweights, A_ECG(1:128));
    if strcmpi(ECGtype, 'ecg')
        [R_ECG2, ~] = corr(Heartweights, R_ECGEEG);
    end
    R_ECG1 = round(R_ECG1, 2);
    R_ECG2 = round(R_ECG2, 2);
end

if iscell(HEOGepochs)
    if has_HEOG_L
        [R_HEOG1_L, ~] = corr(Saccadeweights, A_HEOG_L(1:128));
        [R_HEOG2_L, ~] = corr(Saccadeweights, R_HEOGEEG_L);
        R_HEOG1_L = round(R_HEOG1_L, 2);
        R_HEOG2_L = round(R_HEOG2_L, 2);
    end
    if has_HEOG_R
        [R_HEOG1_R, ~] = corr(Saccadeweights, A_HEOG_R(1:128));
        [R_HEOG2_R, ~] = corr(Saccadeweights, R_HEOGEEG_R);
        R_HEOG1_R = round(R_HEOG1_R, 2);
        R_HEOG2_R = round(R_HEOG2_R, 2);
    end
else
    if has_HEOG
        [R_HEOG1, ~] = corr(Saccadeweights, A_HEOG(1:128));
        [R_HEOG2, ~] = corr(Saccadeweights, R_HEOGEEG);
        R_HEOG1 = round(R_HEOG1, 2);
        R_HEOG2 = round(R_HEOG2, 2);
    end
end

% -------------------------------------------------------------------------
% Calculate safe sample sizes for plotting
% -------------------------------------------------------------------------
if has_VEOG, N_VEOG = size(VEOGepochs, 1); else, N_VEOG = 0; end
if has_ECG,  N_ECG = size(ECGepochs, 1);   else, N_ECG = 0;  end

if iscell(HEOGepochs)
    if has_HEOG_L, N_HEOG1 = size(HEOGepochs{1}, 1); else, N_HEOG1 = 0; end
    if has_HEOG_R, N_HEOG2 = size(HEOGepochs{2}, 1); else, N_HEOG2 = 0; end
else
    if has_HEOG, N_HEOG = size(HEOGepochs, 1); else, N_HEOG = 0; end
end

% -------------------------------------------------------------------------
fprintf('Plotting...\n');

% Plot 1
fh = figure('Visible', cfg.figure.visible);
tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");
% 1
mytopoplot(A_VEOG(1:128),[],['avgVEOG: N = ' num2str(N_VEOG) ', R = ' num2str(R_VEOG1)],nexttile(1));
mytopoplot(A_HEOG_L(1:128) + A_HEOG_R(1:128),[],['avgHEOG: N = ' num2str(N_HEOG1+N_HEOG1) ', R = ' num2str(R_HEOG1_L) '/' num2str(R_HEOG1_R)],nexttile(2));
if ~strcmpi(ECGtype, 'none')
    R_ECG1 = mean(R_ECG1);
    mytopoplot(A_ECG(1:128),[],['avgECG: N = ' num2str(N_ECG) ', R = ' num2str(R_ECG1)],nexttile(3));
end
% 2
mytopoplot(R_VEOGEEG,[],['corrVEOG: N = ' num2str(N_VEOG) ', R = ' num2str(R_VEOG2)],nexttile(4)); colorbar;
mytopoplot(R_HEOGEEG_L+R_HEOGEEG_R,[],['corrHEOG: N = ' num2str(N_HEOG1+N_HEOG1) ', R = ' num2str(R_HEOG2_L) '/' num2str(R_HEOG2_R)],nexttile(5)); colorbar;
if strcmpi(ECGtype, 'ecg')
    mytopoplot(R_ECGEEG,[],['corrECG: N = ' num2str(N_ECG) ', R = ' num2str(mean(R_ECG2))],nexttile(6)); colorbar;
end
% 3
mytopoplot(Blinkweights,[],'templateVEOG',nexttile(7));
mytopoplot(Saccadeweights,[],'templateHEOG',nexttile(8));
mytopoplot(mean(Heartweights,2),[],'templateECG',nexttile(9));

% Save
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_templates1'], [20 20]);

% Plot 2
if iscell(HEOGepochs)
    fh = figure('Visible', cfg.figure.visible); tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact");
    mytopoplot(A_HEOG_L(1:128),[],['avgHEOG1: N = ' num2str(N_HEOG1) ', R = ' num2str(R_HEOG1_L)],nexttile); axis tight;
    mytopoplot(A_HEOG_R(1:128),[],['avgHEOG2: N = ' num2str(N_HEOG2) ', R = ' num2str(R_HEOG1_R)],nexttile); axis tight;
    mytopoplot(R_HEOGEEG_L,[],['corrHEOG1: N = ' num2str(N_HEOG1) ', R = ' num2str(R_HEOG2_L)],nexttile); axis tight;
    mytopoplot(R_HEOGEEG_R,[],['corrHEOG2: N = ' num2str(N_HEOG2) ', R = ' num2str(R_HEOG2_L)],nexttile); axis tight;
    plotX=20; plotY=20;
else
    fh = figure; tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
    mytopoplot(A_HEOG_L(1:128),[],['avgHEOG: N = ' num2str(N_HEOG) ', R = ' num2str(R_HEOG1)],nexttile); axis tight;
    mytopoplot(A_HEOG_L(1:128),[],['corrHEOG: N = ' num2str(N_HEOG) ', R = ' num2str(R_HEOG2)],nexttile); axis tight;
    plotX=20; plotY=10;
end

% Save
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_templates2'], [plotX plotY]);

% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
wArtifacts = [];

wArtifacts.Blinkmask        = VEOGmask;
wArtifacts.Heartmask        = ECGmask;
if iscell(HEOGepochs)
    wArtifacts.SaccademaskL = HEOGmask{1};
    wArtifacts.SaccademaskR = HEOGmask{2};
    wArtifacts.SaccadeepochsL = HEOGepochs{1};
    wArtifacts.SaccadeepochsR = HEOGepochs{2};
else
    wArtifacts.SaccademaskL = HEOGmask;
    wArtifacts.SaccademaskR = [];
    wArtifacts.SaccadeepochsL = HEOGepochs;
    wArtifacts.SaccadeepochsR = [];
end

wArtifacts.Blinkepochs1     = VEOGepochs;
wArtifacts.Heartepochs1     = ECGepochs;

wArtifacts.Blinkweights1    = A_VEOG;
wArtifacts.Blinkweights2    = R_VEOGEEG;
wArtifacts.Heartweights1    = A_ECG;
wArtifacts.Heartweights2    = R_ECGEEG;
wArtifacts.Saccadeweights1L = A_HEOG_L;
wArtifacts.Saccadeweights1R = A_HEOG_R;
wArtifacts.Saccadeweights2L = R_HEOGEEG_L;
wArtifacts.Saccadeweights2R = R_HEOGEEG_R;

wArtifacts.Blinkcov         = C_VEOG;
wArtifacts.Heartcov         = C_ECG;
wArtifacts.SaccadecovL      = C_HEOG_L;
wArtifacts.SaccadecovR      = C_HEOG_R;

wArtifacts.Blinkcorr        = R_VEOG1;
wArtifacts.Heartcorr        = R_ECG1;
wArtifacts.Saccadecorr1L    = R_HEOG1_L;
wArtifacts.Saccadecorr1R    = R_HEOG1_R;
wArtifacts.Saccadecorr2L    = R_HEOG2_L;
wArtifacts.Saccadecorr2R    = R_HEOG2_R;

fileName = [EEG(1).ALSUTRECHT.subject.id '_' EEG(1).ALSUTRECHT.subject.visit '_' EEG(1).ALSUTRECHT.subject.task '_wartifacts.mat'];
save(fullfile(EEG(1).ALSUTRECHT.subject.data, fileName), "wArtifacts");

fprintf('Done!\n');

end

% =========================================================================
% Helper fuctions
% =========================================================================

function [C, A] = estimate_params(DATA, epochs)
y = DATA.data(:, epochs);
C = y * y';
A = mean(y, 2);
end

% function [DATA_tmp, ECGtype] = make_bipolar0(DATA_block)
% DATA_tmp = DATA_block;
%
% eog_labels = {'VEOG', 'HEOG'};
%
% for i = 1:length(eog_labels)
%     current_label = eog_labels{i};
%     mask = find(contains({DATA_tmp.chanlocs.labels}, current_label));
%
%     if length(mask) == 2
%         % fprintf('Found 2 %s electrodes. Creating bipolar channel...\n', current_label);
%         % DATA_tmp.data(mask(1),:) = DATA_tmp.data(mask(1),:) - DATA_tmp.data(mask(2),:);
%         % DATA_tmp.chanlocs(mask(1)).labels = current_label;
%     elseif isscalar(mask)
%         fprintf('%s appears to be already bipolarised (only 1 channel found).\n', current_label);
%         DATA_tmp.chanlocs(mask(1)).labels = current_label;
%     else
%         % warning('Expected 1 or 2 %s electrodes, but found %d.', current_label, length(mask));
%     end
% end
%
% maskECG = find(contains({DATA_tmp.chanlocs.labels}, 'ECG'));
%
% if length(maskECG) == 2
%     % fprintf('Found 2 ECG electrodes. Making one bipolar ECG channel.\n');
%     % ECGtype = 'ECG';
%     % DATA_tmp.data(maskECG(1),:) = DATA_tmp.data(maskECG(1),:) - DATA_tmp.data(maskECG(2),:);
%     % DATA_tmp.chanlocs(maskECG(1)).labels = 'ECG';
% elseif isscalar(maskECG)
%     fprintf('ECG appears to be already bipolarised (only 1 channel found).\n');
%     ECGtype = 'ECG';
%     DATA_tmp.chanlocs(maskECG(1)).labels = 'ECG';
% else
%     % maskEMG = find(strcmp({DATA_tmp.chanlocs.type}, 'EMG'));
%     %
%     % if ~isempty(maskEMG)
%     %     fprintf('No ECG recorded. Using average of %d EMG electrodes as ECG proxy.\n', length(maskEMG));
%     %     ECGtype = 'EMG';
%     %     target_idx = maskEMG(1);
%     %     DATA_tmp.data(target_idx, :) = mean(DATA_tmp.data(maskEMG, :), 1);
%     %     DATA_tmp.chanlocs(target_idx).labels = 'ECG';
%     % else
%     %     fprintf('No ECG or EMG electrodes found. QRS detection will be skipped.\n');
%     %     ECGtype = 'None';
%     % end
% end
% end