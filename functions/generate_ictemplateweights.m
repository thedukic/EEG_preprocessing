function generate_ictemplateweights(EEG, EXT, cfg)

fprintf('\n================================\n');
fprintf('Generating individualised ICA artifact templates\n');
fprintf('================================\n');

% -------------------------------------------------------------------------
% 0. Ensure complete 128-channel montage via temporary interpolation
% -------------------------------------------------------------------------
% Work on a local copy so downstream pipeline stages (e.g. ICA) retain
% their true physical channel rank.

if EEG.nbchan < 128
    fprintf('Detected missing channels (%d/128 present). Interpolating for template generation...\n', EEG.nbchan);
    EEG = do_channelinterp(EEG, 'spherical');
end

% Rerefrence
EEG = do_reref(EEG, 'aRegular');

% -------------------------------------------------------------------------
% Merge blocks
% -------------------------------------------------------------------------
DATA = merge_electrodetypes(EEG, EXT);

cfg_tmp = cfg;
cfg_tmp.flt.erp.hp = [1 4];
cfg_tmp.flt.rs.hp  = [1 4];
cfg_tmp.flt.mt.hp  = [1 4];

DATA = do_filtering(DATA, 'highpass', {'eeg','ext'}, cfg_tmp);

% -------------------------------------------------------------------------
% Eye blinks
% -------------------------------------------------------------------------
fprintf('\n--------------------------------\n');
fprintf('Detecting eye blinks\n');
fprintf('--------------------------------\n');

cfg_tmp.win          = 150;
cfg_tmp.trIQR        = 3;
cfg_tmp.do_plot      = true;
cfg_tmp.plot_visible = cfg.figure.visible;
[VEOGmask, VEOGepochs, VEOGmaxLatency, VEOGdata, VEOGEEGdata, treshold] = detect_veog(DATA, cfg_tmp);

% -------------------------------------------------------------------------
% Eye saccades
% -------------------------------------------------------------------------
fprintf('\n--------------------------------\n');
fprintf('Detecting eye saccades\n');
fprintf('--------------------------------\n');

cfg_tmp.win          = 150;
cfg_tmp.trIQR        = 1;
cfg_tmp.do_plot      = true;
cfg_tmp.plot_visible = cfg.figure.visible;
[HEOGmask, HEOGepochs, HEOGmaxLatency, HEOGdata, HEOGEEGdata, treshold] = detect_heog(DATA, cfg_tmp);

% -------------------------------------------------------------------------
% Heart beats
% -------------------------------------------------------------------------
fprintf('\n--------------------------------\n');
fprintf('Detecting heartbeats\n');
fprintf('--------------------------------\n');

if ~strcmpi(DATA.ALSUTRECHT.subject.ecg, 'none')
    cfg_tmp.win          = 30;
    cfg_tmp.do_plot      = false;
    cfg_tmp.plot_visible = cfg.figure.visible;
    [ECGmask, ECGepochs, ECGlatency, ECGdata, ECGEEGdata, pulsEstimate] = detect_ecg_new(DATA, cfg_tmp);
else
    ECGmask   = NaN;
    ECGepochs = NaN;
    has_ECG   = false;
end

% -------------------------------------------------------------------------
% Setup availability flags to clean up downstream logic
% -------------------------------------------------------------------------
has_VEOG = any(VEOGmask);
if islogical(ECGmask)
    has_ECG = any(ECGmask);
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
C_VEOG   = NaN; A_VEOG   = zeros(128,1); R_VEOGEEG   = zeros(128,1);
C_ECG    = NaN; A_ECG    = zeros(128,1); R_ECGEEG    = zeros(128,1);
C_HEOG   = NaN; A_HEOG   = zeros(128,1); R_HEOGEEG   = zeros(128,1);
C_HEOG_L = NaN; A_HEOG_L = zeros(128,1); R_HEOGEEG_L = zeros(128,1);
C_HEOG_R = NaN; A_HEOG_R = zeros(128,1); R_HEOGEEG_R = zeros(128,1);

fprintf('Estimating artifacts topoplots and covariance matrices...\n');

if has_VEOG
    [C_VEOG, A_VEOG] = estimate_params(DATA, VEOGmask);
    R_VEOGEEG = estimate_extcorrelation(DATA, VEOGmask, 'veog');
end

if has_ECG
    [C_ECG, A_ECG] = estimate_params(DATA, ECGmask);
    R_ECGEEG = estimate_extcorrelation(DATA, ECGmask, 'ecg');
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
R_VEOG1   = NaN; R_VEOG2 = NaN;
R_ECG1    = NaN; R_ECG2  = NaN;
R_HEOG1   = NaN; R_HEOG2 = NaN;
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
    [R_ECG2, ~] = corr(Heartweights, R_ECGEEG);
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
if has_ECG,  N_ECG  = size(ECGepochs, 1);  else, N_ECG  = 0; end

if iscell(HEOGepochs)
    if has_HEOG_L, N_HEOG1 = size(HEOGepochs{1}, 1); else, N_HEOG1 = 0; end
    if has_HEOG_R, N_HEOG2 = size(HEOGepochs{2}, 1); else, N_HEOG2 = 0; end
else
    if has_HEOG, N_HEOG = size(HEOGepochs, 1); else, N_HEOG = 0; end
end

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------
fprintf('Plotting...\n');

if cfg.figure.plot
    % Plot 1
    fh = figure('Visible', cfg.figure.visible);
    tiledlayout(3, 3, "TileSpacing", "compact", "Padding", "compact");

    % Row 1: Average Topographies
    mytopoplot(A_VEOG(1:128), [], ['avgVEOG: N = ' num2str(N_VEOG) ', R = ' num2str(R_VEOG1)], nexttile(1));
    mytopoplot(A_HEOG_L(1:128) + A_HEOG_R(1:128), [], ['avgHEOG: N = ' num2str(N_HEOG1+N_HEOG2) ', R = ' num2str(R_HEOG1_L) '/' num2str(R_HEOG1_R)], nexttile(2));
    R_ECG1 = mean(R_ECG1);
    mytopoplot(A_ECG(1:128), [], ['avgECG: N = ' num2str(N_ECG) ', R = ' num2str(R_ECG1)], nexttile(3));

    % Row 2: Correlation Topographies
    mytopoplot(R_VEOGEEG(1:128), [], ['corrVEOG: N = ' num2str(N_VEOG) ', R = ' num2str(R_VEOG2)], nexttile(4)); colorbar;
    mytopoplot(R_HEOGEEG_L(1:128) + R_HEOGEEG_R(1:128), [], ['corrHEOG: N = ' num2str(N_HEOG1+N_HEOG2) ', R = ' num2str(R_HEOG2_L) '/' num2str(R_HEOG2_R)], nexttile(5)); colorbar;
    mytopoplot(R_ECGEEG(1:128), [], ['corrECG: N = ' num2str(N_ECG) ', R = ' num2str(mean(R_ECG2))], nexttile(6)); colorbar;

    % Row 3: Standard Templates
    mytopoplot(Blinkweights, [], 'templateVEOG', nexttile(7));
    mytopoplot(Saccadeweights, [], 'templateHEOG', nexttile(8));
    mytopoplot(mean(Heartweights, 2), [], 'templateECG', nexttile(9));

    % Save
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_templates1'], [20 20]);

    % Plot 2
    if iscell(HEOGepochs)
        fh = figure('Visible', cfg.figure.visible);
        tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "compact");
        mytopoplot(A_HEOG_L(1:128), [], ['avgHEOG1: N = ' num2str(N_HEOG1) ', R = ' num2str(R_HEOG1_L)], nexttile); axis tight;
        mytopoplot(A_HEOG_R(1:128), [], ['avgHEOG2: N = ' num2str(N_HEOG2) ', R = ' num2str(R_HEOG1_R)], nexttile); axis tight;
        mytopoplot(R_HEOGEEG_L(1:128), [], ['corrHEOG1: N = ' num2str(N_HEOG1) ', R = ' num2str(R_HEOG2_L)], nexttile); axis tight;
        mytopoplot(R_HEOGEEG_R(1:128), [], ['corrHEOG2: N = ' num2str(N_HEOG2) ', R = ' num2str(R_HEOG2_R)], nexttile); axis tight;
        plotX = 20; plotY = 20;
    else
        fh = figure('Visible', cfg.figure.visible);
        tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
        mytopoplot(A_HEOG(1:128), [], ['avgHEOG: N = ' num2str(N_HEOG) ', R = ' num2str(R_HEOG1)], nexttile); axis tight;
        mytopoplot(R_HEOGEEG(1:128), [], ['corrHEOG: N = ' num2str(N_HEOG) ', R = ' num2str(R_HEOG2)], nexttile); axis tight;
        plotX = 20; plotY = 10;
    end

    % Save
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_templates2'], [plotX plotY]);
end

% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
wArtifacts = [];

wArtifacts.Blinkmask          = VEOGmask;
wArtifacts.Heartmask          = ECGmask;
if iscell(HEOGepochs)
    wArtifacts.SaccademaskL   = HEOGmask{1};
    wArtifacts.SaccademaskR   = HEOGmask{2};
    wArtifacts.SaccadeepochsL = HEOGepochs{1};
    wArtifacts.SaccadeepochsR = HEOGepochs{2};
else
    wArtifacts.SaccademaskL   = HEOGmask;
    wArtifacts.SaccademaskR   = [];
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

file_name = [EEG(1).ALSUTRECHT.subject.id '_' EEG(1).ALSUTRECHT.subject.visit '_' EEG(1).ALSUTRECHT.subject.task '_wartifacts.mat'];
save(fullfile(EEG(1).ALSUTRECHT.subject.data, file_name), "wArtifacts");

fprintf('Done!\n');

end

% =========================================================================
% Helper functions
% =========================================================================

function [C, A] = estimate_params(DATA, epochs)
y = DATA.data(:, epochs);
C = y * y';
A = mean(y, 2);
end