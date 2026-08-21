function EEG = remove_badcomponents(EEG,cfg)
% REMOVE_BADCOMPONENTS
% Evaluates, filters, and removes artifactual Independent Components.
% Applies peripherality/motor-protection checks strictly to candidate Muscle ICs.

fprintf('\n================================\n');
fprintf('Evaluating and removing bad ICs\n');
fprintf('================================\n');

% Ensure ICA structure is validk
EEG = eeg_checkset(EEG, 'ica');
num_ica = size(EEG.icaweights, 1);

% Extract candidate bad ICs (ensure column vectors)
ICsEye     = EEG.ALSUTRECHT.ica.final.eye(:);
ICsComplex = EEG.ALSUTRECHT.ica.final.complex(:);
ICsHeart   = EEG.ALSUTRECHT.ica.final.heart(:);
ICsMuscle  = EEG.ALSUTRECHT.ica.final.muscle(:);
ICsChannel = EEG.ALSUTRECHT.ica.final.channel(:);
ICsGenBad  = EEG.ALSUTRECHT.ica.final.genbad(:);

% -------------------------------------------------------------------------
% Dynamic Variance & Rank Evaluation (BioSemi 128 Montage)
% -------------------------------------------------------------------------
floor_compvar  = 0.5;  % PVAF floor (%): ICs with < 0.5% compvar are spared past rank 25
protect_rank_n = 25;   % Top 25 components are always evaluated regardless of PVAF

if isfield(cfg, 'ica')
    if isfield(cfg.ica, 'floor_compvar');  floor_compvar  = cfg.ica.floor_compvar;  end
    if isfield(cfg.ica, 'protect_rank_n'); protect_rank_n = cfg.ica.protect_rank_n; end
end

fprintf('\nEvaluating IC importance based on EEGLAB compvar PVAF...\n');
fprintf('Threshold: vaf_compvar >= %.2f%% (Always evaluate top %d ICs)\n', ...
    floor_compvar, protect_rank_n);

% Extract vaf_compvar (calculates if not present)
vaf_compvar = EEG.ALSUTRECHT.ica.vaf_compvar(:);
num_ics  = length(vaf_compvar);
ic_ranks = (1:num_ics)';

% Spares components that have low PVAF (< 0.25%) AND are outside top 25 rank
is_negligible = (vaf_compvar < floor_compvar) & (ic_ranks > protect_rank_n);

% Track spared components
cnt_m = sum(ICsMuscle  & is_negligible);
cnt_c = sum(ICsChannel & is_negligible);
cnt_b = sum(ICsGenBad  & is_negligible);

num_emg_0     = sum(ICsMuscle);
num_channel_0 = sum(ICsChannel);
num_bad_0     = sum(ICsGenBad);

% Exclude negligible components from candidate lists
ICsMuscle(is_negligible)  = false;
ICsChannel(is_negligible) = false;
ICsGenBad(is_negligible)  = false;

fprintf('Spared low-PVAF late ICs (vaf_compvar < %.2f%% & Rank > %d):\n', floor_compvar, protect_rank_n);
fprintf('  Muscle ICs spared  : %d / %d\n', cnt_m, num_emg_0);
fprintf('  Channel ICs spared : %d / %d\n', cnt_c, num_channel_0);
fprintf('  GenBad ICs spared  : %d / %d\n', cnt_b, num_bad_0);

% -------------------------------------------------------------------------
% Check likelihood of EMG ICs
% -------------------------------------------------------------------------
ICsMuscle_tmp = ICsMuscle;
ICsMuscle_tmp(30:end) = false;

% Check
[EEG, safe_to_remove] = check_ic_peripherality(EEG, ICsMuscle_tmp, cfg);

% % Prune muscle candidates that overlap sensorimotor cortex or central scalp
% num_emg_before_periph = sum(ICsMuscle);
% num_bad_before_periph = sum(ICsGenBad);
%
% ICsMuscle_safe = ICsMuscle & safe_to_remove(:);
% ICsGenBad_safe = ICsGenBad & safe_to_remove(:);
%
% protected_emg_cnt = num_emg_before_periph - sum(ICsMuscle_safe);
% protected_bad_cnt = num_bad_before_periph - sum(ICsGenBad_safe);
%
% fprintf('Muscle ICs evaluated: %d candidate(s) -> %d verified peripheral (%d protected over motor/central areas).\n', ...
%     num_emg_before_periph, sum(ICsMuscle_safe), protected_emg_cnt);
% fprintf('GenBad ICs evaluated: %d candidate(s) -> %d verified peripheral (%d protected over motor/central areas).\n', ...
%     num_bad_before_periph, sum(ICsGenBad_safe), protected_bad_cnt);

ICsMuscle_safe = ICsMuscle;
ICsGenBad_safe = ICsGenBad;

% -------------------------------------------------------------------------
% Check ICs using dipole fitting
% -------------------------------------------------------------------------
ICsforRemove = false(num_ica, 1);
ICsforRemove(ICsEye | ICsHeart | ICsComplex | ICsMuscle_safe | ICsChannel | ICsGenBad_safe) = true;

[EEG, inside_brain, good_fits] = fit_ic_dipoles(EEG, find(ICsforRemove), cfg);

% ICsEye(inside_brain)         = false; % eyes are sometimes
% ICsComplex(inside_brain)     = false;
ICsHeart(inside_brain)       = false;
ICsMuscle_safe(inside_brain) = false;
ICsChannel(inside_brain)     = false;
ICsGenBad_safe(inside_brain) = false;

% -------------------------------------------------------------------------
% 3. Compile Master Removal Vector
% -------------------------------------------------------------------------
ICsforRemove = false(num_ica, 1);

% Non-negotiable removals (Ocular, Cardiac, Complex)
ICsforRemove(ICsEye | ICsHeart | ICsComplex) = true;

% Muscle ICs (Controlled via cfg.ica.emg)
if isfield(cfg.ica, 'emg') && cfg.ica.emg
    fprintf('Muscle components (N = %d) will be removed.\n', sum(ICsMuscle_safe));
    ICsforRemove(ICsMuscle_safe) = true;
else
    fprintf('Muscle components (N = %d) will be retained.\n', sum(ICsMuscle_safe));
end

% Channel ICs (Controlled via cfg.ica.channel)
if isfield(cfg.ica, 'channel') && cfg.ica.channel
    fprintf('Channel components (N = %d) will be removed.\n', sum(ICsChannel));
    ICsforRemove(ICsChannel) = true;
else
    fprintf('Channel components (N = %d) will be retained.\n', sum(ICsChannel));
end

% GenBad ICs (Controlled via cfg.ica.bad)
if isfield(cfg.ica, 'bad') && cfg.ica.bad
    fprintf('Generally bad components (N = %d) will be removed.\n', sum(ICsGenBad_safe));
    ICsforRemove(ICsGenBad_safe) = true;
else
    fprintf('Generally bad components (N = %d) will be retained.\n', sum(ICsGenBad_safe));
end

% Summary Report to Console
fprintf('\n--------------------------------------------------\n');
fprintf('Component Rejection Breakdown:\n');
fprintf('  Eye ICs:         %2d\n', sum(ICsEye));
fprintf('  Heart ICs:       %2d\n', sum(ICsHeart));
fprintf('  Complex ICs:     %2d\n', sum(ICsComplex));
fprintf('  Muscle ICs:      %2d (Safe Peripheral)\n', sum(ICsMuscle_safe));
fprintf('  Channel ICs:     %2d\n', sum(ICsChannel));
fprintf('  Gen. Bad ICs:    %2d\n', sum(ICsGenBad_safe));
fprintf('  --------------------\n');
fprintf('  Total Removed:   %2d / %d ICs\n', sum(ICsforRemove), num_ica);
fprintf('--------------------------------------------------\n\n');

% -------------------------------------------------------------------------
% 4. Compute True Channel Variance Removed & Diagnostic Check
% -------------------------------------------------------------------------
ch_idx   = EEG.ALSUTRECHT.ica.icachansind;
eeg_data = double(reshape(EEG.data(ch_idx, :, :), length(ch_idx), []));
ica_data = (EEG.ALSUTRECHT.ica.icaweights * EEG.ALSUTRECHT.ica.icasphere) * eeg_data;

if any(ICsforRemove)
    [~, true_pct_var_removed] = compvar(eeg_data, ica_data, EEG.ALSUTRECHT.ica.icawinv, find(ICsforRemove));
    fprintf('True scalp channel variance removed: %.2f%%\n', true_pct_var_removed);

    % Power spectrum verification plot
    check_power_cleaning(EEG, [], ICsforRemove, 'ica', cfg);

    % Back-project clean data (0 avoids interactive pop-up window)
    EEG = pop_subcomp(EEG, find(ICsforRemove), 0);
else
    fprintf('No components flagged for removal. Dataset unchanged.\n');
end

% Re-verify EEGLAB dataset structure
EEG = eeg_checkset(EEG);

% % =========================================================================
% % Must-remove ICs
% if any(ICsforRemove)
%     fprintf('-> Removing bad ICs (N = %d)...\n', sum(ICsforRemove));
%     artifactComponents(ICsforRemove, :) = dataICs(ICsforRemove, :);
% end
%
% % -------------------------------------------------------------------------
% % Muscle ICs
% if cfg.ica.emg
%     if any(ICsMostLikelyMuscle)
%         fprintf('-> Filtering muscle ICs (N = %d)...\n', sum(ICsMostLikelyMuscle));
%         [bh, ah] = butter(2, muscleFreqCutoff/(EEG.srate/2), 'high');
%         artifactComponents(ICsMostLikelyMuscle, :) = do_filteringcore(bh, ah, dataICs(ICsMostLikelyMuscle, :), EEG.event, EEG.srate);
%     else
%         fprintf('-> Muscle ICs are not found.\n');
%     end
% else
%     fprintf('-> Muscle ICs (N = %d) are kept in the EEG.\n', sum(ICsMostLikelyMuscle));
% end
%
% % -------------------------------------------------------------------------
% % Obtain channel artifact for subtraction by wavelet-thresholding
% % -> this did not work well
% % if any(ICsMostLikelyChannel)
% %     ICsMostLikelyChannel2 = find(ICsMostLikelyChannel);
% %
% %     for i = 1:length(ICsMostLikelyChannel2)
% %         artifactComponents(ICsMostLikelyChannel2(i),:) = wThresholding(dataICs(ICsMostLikelyChannel2(i),:));
% %     end
% % end
%
% % -------------------------------------------------------------------------
% % % Channel ICs
% % % Regress out channel ICs
% % if cfg.ica.channel
% %     if any(ICsMostLikelyChannel)
% %         fprintf('-> Removing channel ICs (N = %d)...\n', sum(ICsMostLikelyChannel));
% %         artifactComponents(ICsMostLikelyChannel,:) = dataICs(ICsMostLikelyChannel,:);
% %     else
% %         fprintf('-> Channel ICs are not found.\n');
% %     end
% % else
% %     fprintf('-> Channel ICs (N = %d) are kept in the EEG.\n', sum(ICsMostLikelyChannel));
% % end
% %
% % EEGNEW = EEG;
% % EEGNEW.data(chaneeg,:) = cleanEEG';
% % vis_artifacts(EEGNEW,EEG);
%
% % -------------------------------------------------------------------------
% % Remove artifact and reconstruct data
% artifactEEG = EEG.icawinv * artifactComponents;
% artifactEEG = reshape(artifactEEG, EEG.nbchan, NTPT, EEG.trials);
%
% % chaneeg = strcmp({EEG.chanlocs.type},'EEG');
% % EEGNEW = EEG;
% % EEGNEW.data(chaneeg,:) = EEG.data(chaneeg,:) - artifactEEG;
% % vis_artifacts(EEGNEW,EEG);
%
% chaneeg = strcmp({EEG.chanlocs.type}, 'EEG');
% EEG.data(chaneeg, :) = EEG.data(chaneeg, :) - artifactEEG;

% -------------------------------------------------------------------------
% 5. Audit Logging
% -------------------------------------------------------------------------
EEG.ALSUTRECHT.ica.final.removed        = ICsforRemove;
EEG.ALSUTRECHT.ica.final.all_bad_marked = ICsEye | ICsHeart | ICsComplex | ICsMuscle_safe | ICsChannel | ICsGenBad_safe;
EEG.ALSUTRECHT.ica.final.protected_emg  = ICsMuscle & ~safe_to_remove(:);

if isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'subject') && ...
        isfield(EEG.ALSUTRECHT.subject, 'fid') && ~isempty(EEG.ALSUTRECHT.subject.fid) && ...
        EEG.ALSUTRECHT.subject.fid > 0

    fid = EEG.ALSUTRECHT.subject.fid;
    fprintf(fid, '\n---------------------------------------------------------\n');
    fprintf(fid, 'ICA Bad Component Removal Summary\n');
    fprintf(fid, '---------------------------------------------------------\n');
    fprintf(fid, 'Total ICs Removed:           %d\n', sum(ICsforRemove));
    fprintf(fid, '  - Eye:                     %d\n', sum(ICsEye));
    fprintf(fid, '  - Heart:                   %d\n', sum(ICsHeart));
    fprintf(fid, '  - Complex:                 %d\n', sum(ICsComplex));
    fprintf(fid, '  - Muscle:                  %d\n', sum(ICsMuscle_safe));
    fprintf(fid, '  - Channel:                 %d\n', sum(ICsChannel));
    fprintf(fid, '  - GenBad:                  %d\n', sum(ICsGenBad_safe));
    % fprintf(fid, 'Protected Motor EMG ICs:     %d\n', protected_emg_cnt);
    % fprintf(fid, 'Protected Motor GenBad ICs:  %d\n', protected_bad_cnt);
    if any(ICsforRemove)
        fprintf(fid, 'Variance Accounted For:   %.2f%%\n', true_pct_var_removed);
    end
end

% Clear activations cache to conserve memory
EEG.icaact = [];

end

% % =========================================================================
% % Helper function
% % =========================================================================
% function plot_ic_diagnostics(EEG, ic_idx)
% % plot_ic_diagnostics(EEG, ic_idx)
% % Visualizes the Topography, PSD, and a zoomed-in Time Series of a specific IC.
%
% fprintf('Visualizing IC %d...\n', ic_idx);
%
% % 1. Reconstruct the activation for the specific IC
% % We only need to compute the time series for this one component, saving memory
% ic_act = (EEG.icaweights(ic_idx, :) * EEG.icasphere) * EEG.data(EEG.icachansind, :);
%
% % 2. Setup Figure
% fh = figure('Name', sprintf('IC %d Diagnostics', ic_idx), 'Color', 'w');
% % Use a 2x2 grid. Top row for Topo/PSD, Bottom row spanning across for Time Series
% t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
%
% % ---------------------------------------------------------------------
% % Tile 1: Topography
% nexttile(1);
% % Standard EEGLAB topoplot using the inverse weights
% topoplot(EEG.icawinv(:, ic_idx), EEG.chanlocs(EEG.icachansind), 'electrodes', 'off');
% title(sprintf('IC %d Topography', ic_idx));
% colorbar;
%
% % ---------------------------------------------------------------------
% % Tile 2: Power Spectral Density (PSD)
% nexttile(2);
% % Welch's method (2-second windows, 50% overlap for a smooth spectrum)
% window = EEG.srate * 2;
% noverlap = EEG.srate;
% [pxx, f] = pwelch(ic_act, window, noverlap, window, EEG.srate);
%
% plot(f, 10*log10(pxx), 'k', 'LineWidth', 1.5);
% xlim([1 60]); % Focus on physiological EEG range
% xlabel('Frequency (Hz)');
% ylabel('Power (10*log_{10}(\muV^2/Hz))');
% title('Power Spectral Density');
% grid on;
%
% % ---------------------------------------------------------------------
% % Tile 3: Time Series (Zoomed on Max Peak)
% nexttile(3, [1 2]); hold on;
%
% % Find the absolute maximum peak in the component
% [~, max_idx] = max(abs(ic_act));
%
% % Create a +/- 5 second window around the peak
% win_half = 5 * EEG.srate;
% idx_start = max(1, max_idx - win_half);
% idx_end   = min(length(ic_act), max_idx + win_half);
%
% t_win = (idx_start:idx_end) / EEG.srate;
% data_win = ic_act(idx_start:idx_end);
%
% plot(t_win, data_win, 'Color', [0.2 0.4 0.7], 'LineWidth', 1);
% xline(max_idx / EEG.srate, 'r--', 'Max Deflection', 'LabelVerticalAlignment', 'bottom');
%
% xlabel('Time (s)');
% ylabel('IC Activation');
% title('Time Series (10s window around maximum deflection)');
% xlim([t_win(1) t_win(end)]);
% end
