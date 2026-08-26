function generate_ictemplateweights_old(EEG, EMG, EXT, cfg)

% winBlink = 150;
% winSaccades = 200;
% winECG = 30;

fprintf('\n================================\n');
fprintf('Generating individualised ICA artifact templates\n');
fprintf('================================\n');

% -------------------------------------------------------------------------
% % Merge blocks
% DATA_tmp1 = merge_eeglabblocks(DATA);
%
% % Make external electrodes bipolar
% [DATA_tmp2, ECGtype] = make_bipolar0(DATA_tmp1);
%
% % % Interpolate bad electrodes temporarily
% % % All channels are needed for generating templates
% % DATA_tmp2.ALSUTRECHT.badchaninfo.badElectrodes = DATA(1).ALSUTRECHT.badchaninfo.flatElectrodes;
% % DATA_tmp2 = do_channelinterp(DATA_tmp2,'spherical');

% -------------------------------------------------------------------------
% Merge blocks
DATA_tmp1 = merge_electrodetypes(EEG, EMG, EXT);

% Make external electrodes bipolar
[DATA_tmp2, ECGtype] = make_bipolar0(DATA_tmp1);

% % Interpolate bad electrodes temporarily
% % All channels are needed for generating templates
% DATA_tmp2.ALSUTRECHT.badchaninfo.badElectrodes = DATA(1).ALSUTRECHT.badchaninfo.flatElectrodes;
% DATA_tmp2 = do_channelinterp(DATA_tmp2,'spherical');

% -------------------------------------------------------------------------
% Eye blinks
fprintf('\n--------------------------------\n');
fprintf('Detecting eye blinks\n');
fprintf('--------------------------------\n');

trIQRblink = 3;
winBlink = 150;
[VEOGmask, VEOGepochs, VEOGmaxLatency, VEOGdata, VEOGEEGdata, treshold] = detect_veog(DATA_tmp2, winBlink, trIQRblink, cfg.figure.visible);

% Eye saccades
fprintf('\n--------------------------------\n');
fprintf('Detecting eye saccades\n');
fprintf('--------------------------------\n');

trIQRsaccade = 1;
winSaccade = 150;
[HEOGmask, HEOGepochs, HEOGmaxLatency, HEOGdata, HEOGEEGdata, treshold] = detect_heog(DATA_tmp2, winSaccade, trIQRsaccade, cfg.figure.visible);

% Heart beats
fprintf('\n--------------------------------\n');
fprintf('Detecting heartbeats\n');
fprintf('--------------------------------\n');

if ~strcmpi(ECGtype, 'none')
    winHeart = 30;
    [ECGmask, ECGepochs, ECGlatency, ECGdata, ECGEEGdata, pulsEstimate] = detect_ecg(DATA_tmp2, winHeart, ECGtype, cfg.figure.visible);

    % Did the method detect any heartbeats?
    if ~islogical(ECGmask)
        assert(isnan(ECGmask), 'Strange this should be NaN in this case.');
        ECGtype = 'none';
    end
else
    fprintf('ECG not recorded. Skipping...\n');
    ECGdata = [];
    ECGmask = [];
    ECGepochs = [];
end

% -------------------------------------------------------------------------
% % Plot
% plot_detections(VEOGdata, VEOGepochs, 'VEOG amplitude');
% plot_detections(HEOGdata, HEOGepochs , 'HEOG amplitude');
% plot_detections(ECGdata, ECGepochs, 'ECG amplitude');
% Eye saccades
fprintf('\n--------------------------------\n');
fprintf('Estimating individualised templates\n');
fprintf('--------------------------------\n');

% -------------------------------------------------------------------------
% Estimate artifact topo/cov
% -------------------------------------------------------------------------
fprintf('Estimating artifacts topoplots and covariance matrices...\n');
% [C_VEOG, W_VEOG] = estimate_params(DATA_tmp1, VEOGepochs);
% [C_ECG, W_ECG]   = estimate_params(DATA_tmp1, ECGepochs);
%
% if iscell(HEOGepochs)
%     assert(length(HEOGepochs) == 2);
%     [C_HEOG_L, W_HEOG_L] = estimate_params(DATA_tmp1, HEOGepochs{1});
%     [C_HEOG_R, W_HEOG_R] = estimate_params(DATA_tmp1, HEOGepochs{2});
% else
%     [C_HEOG, W_HEOG] = estimate_params(DATA_tmp1, HEOGepochs);
% end

if any(VEOGmask)
    [C_VEOG, A_VEOG] = estimate_params(DATA_tmp1, VEOGmask);
else
    C_VEOG = NaN;
    A_VEOG = NaN;
end

if iscell(HEOGepochs)
    assert(length(HEOGepochs) == 2);
    if any(HEOGmask{1})
        [C_HEOG_L, A_HEOG_L] = estimate_params(DATA_tmp1, HEOGmask{1});
    else
        C_HEOG_L = NaN;
        A_HEOG_L = NaN;
    end
    if any(HEOGmask{2})
        [C_HEOG_R, A_HEOG_R] = estimate_params(DATA_tmp1, HEOGmask{2});
    else
        C_HEOG_R = NaN;
        A_HEOG_R = NaN;
    end
else
    if any(HEOGmask)
        [C_HEOG, A_HEOG] = estimate_params(DATA_tmp1, HEOGmask);
    else
        C_HEOG = NaN;
        A_HEOG = NaN;
    end
end

if ~strcmpi(ECGtype, 'none')
    [C_ECG, A_ECG] = estimate_params(DATA_tmp1, ECGmask);
else
    C_ECG = NaN;
    A_ECG = NaN;
end

% -------------------------------------------------------------------------
% Estimate EEG-EXT correlations
% -------------------------------------------------------------------------
fprintf('EEG-EXT channel correlations...\n');

if any(VEOGmask)
    R_VEOGEEG = estimate_extcorrelation(DATA_tmp2, VEOGmask, 'veog');
else
    R_VEOGEEG = NaN;
end

if iscell(HEOGepochs)
    if any(HEOGmask{1})
        R_HEOGEEG_L = estimate_extcorrelation(DATA_tmp2, HEOGmask{1}, 'heog');
    else
        R_HEOGEEG_L = NaN;
    end
    if any(HEOGmask{2})
        R_HEOGEEG_R = estimate_extcorrelation(DATA_tmp2, HEOGmask{2}, 'heog');
    else
        R_HEOGEEG_R = NaN;
    end
else
    if any(HEOGmask)
        R_HEOGEEG = estimate_extcorrelation(DATA_tmp2, HEOGmask, 'heog');
    else
        R_HEOGEEG = NaN;
    end
end

if strcmpi(ECGtype, 'ecg')
    % R_ECGEEG = estimate_extcorrelation(DATA_tmp2, [], 'ecg');
    R_ECGEEG = estimate_extcorrelation(DATA_tmp2, ECGmask, 'ecg');
else
    R_ECGEEG = NaN;
end

% -------------------------------------------------------------------------
% Load group tempaltes
% -------------------------------------------------------------------------
fprintf('Comparing the estimates with the group average templates...\n');

load('wBlink.mat',   'Blinkweights');
load('wSaccade.mat', 'Saccadeweights');
load('wHeart.mat',   'Heartweights');

if ~isnan(R_HEOGEEG)
[R_VEOG1, P_VEOG] = corr(Blinkweights, A_VEOG(1:128));
else
    R_VEOG1 = NaN;
end
if ~isnan(R_VEOGEEG)
[R_VEOG2, P_VEOG] = corr(Blinkweights, R_VEOGEEG);
else
    R_VEOG2 = NaN;
end

if ~strcmpi(ECGtype, 'none')
    [R_ECG1, P_ECG] = corr(Heartweights, A_ECG(1:128));
else
    R_ECG1 = NaN;
end
if strcmpi(ECGtype, 'ecg')
    [R_ECG2, P_ECG] = corr(Heartweights, R_ECGEEG);
else
    R_ECG2 = NaN;
end

R_VEOG1 = round(R_VEOG1,2);
R_VEOG2 = round(R_VEOG2,2);
R_ECG1  = round(R_ECG1,2);
R_ECG2  = round(R_ECG2,2);

if iscell(HEOGepochs)
    [R_HEOG1_L, P_HEOG] = corr(Saccadeweights, A_HEOG_L(1:128));
    [R_HEOG1_R, P_HEOG] = corr(Saccadeweights, A_HEOG_R(1:128));
    [R_HEOG2_L, P_HEOG] = corr(Saccadeweights, R_HEOGEEG_L);
    [R_HEOG2_R, P_HEOG] = corr(Saccadeweights, R_HEOGEEG_R);

    R_HEOG1_L = round(R_HEOG1_L,2);
    R_HEOG1_R = round(R_HEOG1_R,2);
    R_HEOG2_L = round(R_HEOG2_L,2);
    R_HEOG2_R = round(R_HEOG2_R,2);
else
    [R_HEOG1, P_HEOG] = corr(Saccadeweights, A_HEOG(1:128));
    [R_HEOG2, P_HEOG] = corr(Saccadeweights, R_HEOGEEG);
    R_HEOG1 = round(R_HEOG1,2);
    R_HEOG2 = round(R_HEOG2,2);
end

% -------------------------------------------------------------------------
fprintf('Plotting...\n');

% Plot 1
N_VEOG = size(VEOGepochs,1);
if ~strcmpi(ECGtype, 'none')
    N_ECG = size(ECGepochs,1);
end
if iscell(HEOGepochs)
    N_HEOG1 = size(HEOGepochs{1},1);
    N_HEOG2 = size(HEOGepochs{2},1);
else
    N_HEOG = size(HEOGepochs,1);
end

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

plotX=20; plotY=20;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh, fullfile(EEG(1).ALSUTRECHT.subject.figures, [EEG(1).ALSUTRECHT.subject.id '_ica_templates1']),'-dtiff','-r200');
close(fh);

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

set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh, fullfile(EEG(1).ALSUTRECHT.subject.figures, [EEG(1).ALSUTRECHT.subject.id '_ica_templates2']), '-dtiff', '-r200');
close(fh);

% -------------------------------------------------------------------------
% Save data 1
wArtifacts = [];

wArtifacts.Blinkmask        = VEOGmask;
wArtifacts.Heartmask        = ECGmask;
wArtifacts.SaccademaskL     = HEOGmask{1};
wArtifacts.SaccademaskR     = HEOGmask{2};

wArtifacts.Blinkepochs1     = VEOGepochs;
wArtifacts.Heartepochs1     = ECGepochs;
wArtifacts.SaccadeepochsL   = HEOGepochs{1};
wArtifacts.SaccadeepochsR   = HEOGepochs{2};

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

% -------------------------------------------------------------------------
% % Save data 2
% extSignals = [];
% extSignals.VEOGdata = VEOGdata;
% extSignals.HEOGdata = HEOGdata;
% extSignals.ECGdata  = ECGdata;
%
% save(fullfile(EEG(1).ALSUTRECHT.subject.figures, [EEG(1).ALSUTRECHT.subject.id '_extSignals.mat']), "extSignals");

fprintf('Done!\n');

end

% =========================================================================
% Helper fuctions
% =========================================================================
% %
% function plot_detections(data, epochs, label)
% if ~isempty(epochs) && ~any(isnan(data))
%     N = size(epochs,1);
%
%     fh = figure; hold on;
%     T = linspace(-0.4, 0.4, length(epochs(1,1):epochs(1,2)));
%
%     for i = 1:N
%         y = data(epochs(i,1):epochs(i,2));
%         plot(T, y, 'LineWidth', 1.2);
%     end
%
%     set(gca,'ColorOrder', brewermap(N,'BuGn'));
%     title(['N = ' num2str(N)]);
%     xlabel('Time a.u.'); ylabel(label);
% end
% end

function [C, A] = estimate_params(DATA, epochs)
% N = size(epochs,1);
% C = NaN(DATA.nbchan,DATA.nbchan,N);
% W = NaN(DATA.nbchan,N);
%
% for i = 1:N
%     y = DATA.data(:, epochs(i,1):epochs(i,2));
%     % y = y - mean(y,1); % avg ref
%     % y = y - mean(y,2); % baseline
%     C(:,:,i) = y * y';
%     W(:,i) = mean(y,2);
% end
% C = mean(C,3);
% W = mean(W,2);

y = DATA.data(:, epochs);
% y = y - mean(y,1); % avg ref
% y = y - mean(y,2); % baseline
C = y * y';
A = mean(y, 2);

end

% function [DATA_tmp, ECGtype] = make_bipolar0(DATA_block)
% % Make a copy
% DATA_tmp = DATA_block;
%
% % 1. Bipolar VEOG is needed for blink detection
% maskEOG = find(contains({DATA_tmp.chanlocs.labels}, 'VEOG'));
% assert(length(maskEOG) == 2);
%
% DATA_tmp.chanlocs(maskEOG(1)).labels = 'VEOG';
% DATA_tmp.data(maskEOG(1),:) = DATA_tmp.data(maskEOG(1),:) - DATA_tmp.data(maskEOG(2),:);
%
% % 2. Bipolar VEOG is needed for blink detection
% maskEOG = find(contains({DATA_tmp.chanlocs.labels}, 'HEOG'));
% assert(length(maskEOG) == 2);
%
% DATA_tmp.chanlocs(maskEOG(1)).labels = 'HEOG';
% DATA_tmp.data(maskEOG(1),:) = DATA_tmp.data(maskEOG(1),:) - DATA_tmp.data(maskEOG(2),:);
%
% % 3. Check if ECG signal is recorded
% maskECG = contains({DATA_tmp.chanlocs.labels}, 'ECG');
% if any(maskECG)
%     % Bipolar ECG is needed for QRS detection
%     fprintf('\nTwo ECG electrodes were recorded. Making one bipolar ECG channel.\n');
%     ECGtype = 'ECG';
%
%     maskECG = find(maskECG);
%     assert(length(maskECG) == 2);
%
%     DATA_tmp.chanlocs(maskECG(1)).labels = 'ECG';
%     DATA_tmp.data(maskECG(1),:) = DATA_tmp.data(maskECG(1),:) - DATA_tmp.data(maskECG(2),:);
% else
%     % Backup plan: use EMG signals as they nicely pick up ECG
%     maskEMG = strcmp({DATA_tmp.chanlocs.type},'EMG');
%     if any(maskEMG)
%         fprintf('\nTwo ECG electrodes were not recorded. Making one ECG channel as the average of EMG electrodes.\n');
%         ECGtype = 'EMG';
%
%         DATA_tmp.chanlocs(1).labels = 'ECG';
%         DATA_tmp.data(1,:) = mean(DATA_tmp.data(maskEMG,:), 1);
%     else
%         ECGtype = 'None';
%     end
% end
% end

function [DATA_tmp, ECGtype] = make_bipolar0(DATA_block)
% Make a copy
DATA_tmp = DATA_block;

% -------------------------------------------------------------------------
% 1 & 2. Handle EOG (VEOG and HEOG)
% -------------------------------------------------------------------------
eog_labels = {'VEOG', 'HEOG'};

for i = 1:length(eog_labels)
    current_label = eog_labels{i};
    mask = find(contains({DATA_tmp.chanlocs.labels}, current_label));

    if length(mask) == 2
        % Case: Unipolar electrodes found, need to bipolarize
        fprintf('Found 2 %s electrodes. Creating bipolar channel...\n', current_label);
        DATA_tmp.data(mask(1),:) = DATA_tmp.data(mask(1),:) - DATA_tmp.data(mask(2),:);
        DATA_tmp.chanlocs(mask(1)).labels = current_label;

    elseif isscalar(mask)
        % Case: Already bipolarized (or only one recorded)
        fprintf('%s appears to be already bipolarised (only 1 channel found).\n', current_label);
        DATA_tmp.chanlocs(mask(1)).labels = current_label;

    else
        warning('Expected 1 or 2 %s electrodes, but found %d.', current_label, length(mask));
    end
end

% -------------------------------------------------------------------------
% 3. Handle ECG Signal
% -------------------------------------------------------------------------
maskECG = find(contains({DATA_tmp.chanlocs.labels}, 'ECG'));

if length(maskECG) == 2
    % Standard unipolar setup
    fprintf('Found 2 ECG electrodes. Making one bipolar ECG channel.\n');
    ECGtype = 'ECG';
    DATA_tmp.data(maskECG(1),:) = DATA_tmp.data(maskECG(1),:) - DATA_tmp.data(maskECG(2),:);
    DATA_tmp.chanlocs(maskECG(1)).labels = 'ECG';

elseif isscalar(maskECG)
    % Already bipolar
    fprintf('ECG appears to be already bipolarised (only 1 channel found).\n');
    ECGtype = 'ECG';
    DATA_tmp.chanlocs(maskECG(1)).labels = 'ECG';

else
    % Backup plan: use EMG signals as they nicely pick up ECG
    maskEMG = find(strcmp({DATA_tmp.chanlocs.type}, 'EMG'));

    if ~isempty(maskEMG)
        fprintf('No ECG recorded. Using average of %d EMG electrodes as ECG proxy.\n', length(maskEMG));
        ECGtype = 'EMG';

        % We use the first EMG channel's slot to host our new ECG proxy
        % Note: You might prefer to append a new channel, but this maintains your original structure
        target_idx = maskEMG(1);
        DATA_tmp.data(target_idx, :) = mean(DATA_tmp.data(maskEMG, :), 1);
        DATA_tmp.chanlocs(target_idx).labels = 'ECG';
    else
        fprintf('No ECG or EMG electrodes found. QRS detection will be skipped.\n');
        ECGtype = 'None';
    end
end

end