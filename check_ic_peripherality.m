function [EEG, safe_to_remove] = check_ic_peripherality(EEG, ics_plot, cfg)
% =========================================================================
% CHECK_IC_PERIPHERALITY: Automated Spatial, Spectral & Temporal IC Diagnostics
% =========================================================================
% Multi-feature classification framework designed to robustly separate
% myogenic (EMG) artefacts from genuine cortical rhythms across high-density
% BioSemi 128 recordings.
%
% -------------------------------------------------------------------------
% DIAGNOSTIC METRICS & CHARACTERISTICS:
% -------------------------------------------------------------------------
% 1. SPATIAL METRICS (Scalp Power Distribution)
%    - PEI (Periphery Energy Index):
%        Quantifies the proportion of component variance located on the outer
%        circumference of the electrode montage (polar radius >= radius_thresh).
%        * Captures: Cranial, jaw, facial, and neck muscle insertions located
%          along the scalp borders, as well as peripheral ocular artefacts.
%    - MOI (Motor Overlap Index):
%        Quantifies the proportion of component variance projecting directly onto
%        the left (D-bank / C3) and right (B-bank / C4) sensorimotor cortex.
%        * Captures: Spatial overlap with primary motor (M1) and somatosensory (S1)
%          areas. Serves as a spatial guardrail to protect genuine mu and beta rhythms.
%
% 2. SPECTRAL METRICS (Frequency Domain Profiling)
%    - Spec Ratio (Gamma-to-Beta Power Ratio: 55-95 Hz / 13-30 Hz):
%        Compares high-frequency gamma power against lower-frequency beta power.
%        * Captures: Broadband spectral flattening. Genuine neural activity
%          exhibits steep 1/f power falloff (low ratio), whereas EMG exhibits a
%          flat or elevated high-frequency plateau (high ratio).
%    - Peak Prominence (Mu / Beta Peak Prominence: 8-30 Hz):
%        Measures the distinct dB prominence of spectral peaks above the 1/f floor.
%        * Captures: True synchronized oscillatory brain rhythms (e.g. sensorimotor
%          mu or cortical beta). Prevents misclassifying genuine brain rhythms as noise.
%
% 3. TEMPORAL METRICS (Continuous Activation Dynamics)
%    - Kurtosis (Distribution Tailedness of IC Activation):
%        Measures the presence of heavy-tailed, extreme amplitude outliers.
%        * Captures: Phasic / bursty muscle artefacts (e.g. transient jaw clenches,
%          swallowing, head twitches, saccades). Tonic EMG yields near-Gaussian
%          values (~3-4.5), while bursty EMG drives kurtosis to extreme levels (>5-20+).
%    - HF Roughness (Sample-to-Sample Derivative Variance Ratio):
%        Computes var(diff(act)) / (2 * var(act)), assessing high-frequency point-to-point jitter.
%        * Captures: Tonic / continuous isometric muscle tension (e.g. postural neck
%          stiffening, subtle jaw clenching). Because tonic EMG resembles Gaussian
%          noise in distribution (low kurtosis), roughness acts as the primary safety net
%          by detecting the high-frequency temporal fuzz characteristic of unsynchronised
%          motor-unit interference.
%
% -------------------------------------------------------------------------
% DECISION LOGIC:
% -------------------------------------------------------------------------
% An IC is flagged as safe to reject if ANY of the following conditions are met:
%  1. Spectral Override: High Spec Ratio (> emg_ratio) with NO beta peak.
%     (Overrides central/motor protection when spectral EMG signature is unambiguous).
%  2. Peripheral Spatial Rim: High PEI (>= pei_thresh) with negligible MOI (<= moi_max_thresh).
%  3. Temporal EMG Pattern: Significant temporal burstiness (kurtosis) OR tonic fuzz
%     (roughness) with peripheral presence (PEI >= 0.30) and NO beta peak.
% =========================================================================

if nargin < 2; ics_plot = []; end
if nargin < 3; cfg = struct(); end

if isstruct(ics_plot)
    cfg = ics_plot;
    ics_plot = [];
end
if islogical(ics_plot); ics_plot = find(ics_plot); end

% --- Defaults ---
if ~isfield(cfg, 'radius_thresh');   cfg.radius_thresh   = 0.40; end
if ~isfield(cfg, 'pei_thresh');      cfg.pei_thresh      = 0.60; end
% if ~isfield(cfg, 'moi_max_thresh');  cfg.moi_max_thresh  = 0.08; end
% if ~isfield(cfg, 'sfi_thresh');      cfg.sfi_thresh      = 0.20; end

% Spectral defaults
if ~isfield(cfg, 'emg_ratio');       cfg.emg_ratio       = 0.25; end
if ~isfield(cfg, 'peak_prom');       cfg.peak_prom       = 2.0;  end

% Temporal defaults for automating activation inspection
if ~isfield(cfg, 'kurtosis_thresh'); cfg.kurtosis_thresh = 5.0;  end
if ~isfield(cfg, 'hf_rough_thresh'); cfg.hf_rough_thresh = 0.35; end

% if ~isfield(cfg, 'roi'); cfg.roi = struct(); end
% if ~isfield(cfg.roi, 'motor_left')
%      cfg.roi.motor_left = {'D19', 'D18', 'D14', 'D20', 'D21', 'D12', 'D11', 'D13', 'D17', 'D28', 'D27'};
% end
% if ~isfield(cfg.roi, 'motor_right')
%      cfg.roi.motor_right = {'B22', 'B21', 'B20', 'B23', 'B24', 'B31', 'B30', 'B32', 'B19', 'B18', 'B17'};
% end

% -------------------------------------------------------------------------
% 1. Data & Spatial Coordinates
% -------------------------------------------------------------------------
assert(isfield(EEG, 'icawinv') && ~isempty(EEG.icawinv), 'EEG.icawinv is empty.');

% Safe fallback for icachansind
if isfield(EEG, 'icachansind') && ~isempty(EEG.icachansind)
    data_chans = EEG.icachansind;
else
    data_chans = 1:size(EEG.icawinv, 1);
end
chanlocs = EEG.chanlocs(data_chans);

n_chan = length(chanlocs);
n_ics  = size(EEG.icawinv, 2);
labels = upper({chanlocs.labels});

if isfield(chanlocs, 'radius') && ~any(cellfun(@isempty, {chanlocs.radius}))
    idx_periphery = find([chanlocs.radius] >= cfg.radius_thresh);
else
    outer_rim_labels = { ...
        'A12', 'A13', 'A14', 'A25', 'A26', 'A27', ...
        'B8',  'B9',  'B10', 'B11', 'B14', 'B15', 'B16', 'B26', 'B27', ...
        'C7',  'C8',  'C16', 'C17', 'C29', 'C30', ...
        'D7',  'D8',  'D23', 'D24', 'D31', 'D32'};
    idx_periphery = find(ismember(labels, outer_rim_labels));
end

if isfield(cfg, 'roi') && isfield(cfg.roi, 'motor_all')
    idx_motor_all = find(ismember(labels, upper(cfg.roi.motor_all)));
elseif isfield(cfg, 'roi') && isfield(cfg.roi, 'motor_left') && isfield(cfg.roi, 'motor_right')
    idx_motor_all = find(ismember(labels, upper([cfg.roi.motor_left, cfg.roi.motor_right])));
else
    idx_motor_all = [];
end

% -------------------------------------------------------------------------
% 2. Temporal & Spectral Metrics from Continuous IC Activations
% -------------------------------------------------------------------------
if isempty(EEG.icaact)
    tmp_data = reshape(EEG.data(data_chans, :, :), length(data_chans), []);
    icaact_2d = (EEG.icaweights * EEG.icasphere) * tmp_data;
else
    icaact_2d = reshape(EEG.icaact, size(EEG.icaact, 1), []);
end

fs = EEG.srate;
% Limit nfft to prevent errors on very short datasets
win_len = min(fs, size(icaact_2d, 2));
nfft = min(1024, size(icaact_2d, 2));
[pxx, f] = pwelch(icaact_2d', win_len, floor(win_len/2), nfft, fs);
pxx = pxx'; % Transpose back to [n_ics x freqs]

idx_mu_beta = f >= 8 & f <= 30;
idx_beta    = f >= 13 & f <= 30;
idx_gamma   = f >= 55 & f <= 95;

% -------------------------------------------------------------------------
% 3. Automated Classification Engine
% -------------------------------------------------------------------------
stats = struct('ic_num', cell(1, n_ics), 'pei', 0, 'moi', 0, 'sfi', 0, ...
    'spec_ratio', 0, 'has_beta_peak', false, 'kurtosis', 0, ...
    'hf_roughness', 0, 'is_temporal_emg', false, 'safe_to_remove', false);
safe_to_remove = false(1, n_ics);

for k = 1:n_ics
    act_k = icaact_2d(k, :);

    % --- Spatial ---
    weights = double(EEG.icawinv(:, k));
    pwr = weights.^2;
    tot_pwr = sum(pwr); if tot_pwr == 0; tot_pwr = eps; end
    pei = sum(pwr(idx_periphery)) / tot_pwr;
    moi = sum(pwr(idx_motor_all)) / tot_pwr;

    % % Spatial Focalness
    % sfi = max(pwr) / tot_pwr;
    % is_focal = sfi > cfg.sfi_thresh;
    sfi = max(pwr) / tot_pwr;

    % --- Spectral ---
    p_beta  = mean(pxx(k, idx_beta));
    p_gamma = mean(pxx(k, idx_gamma));
    spec_ratio = p_gamma / (p_beta + eps);
    has_gamma = spec_ratio > cfg.emg_ratio;

    log_spec = 10 * log10(pxx(k, idx_mu_beta));
    [~, ~, ~, p] = findpeaks(log_spec);
    has_beta_peak = any(p > cfg.peak_prom);

    % --- Temporal Activation Checks ---
    kurt = kurtosis(act_k);

    diff_act = diff(act_k);
    hf_roughness = var(diff_act) / (2 * var(act_k) + eps);

    is_emg_temporal = (kurt > cfg.kurtosis_thresh) || (hf_roughness > cfg.hf_rough_thresh);

    % --- Combined Decision Rules ---
    is_safe = false;
    if has_gamma && is_emg_temporal && ~has_beta_peak
        is_safe = true;
    end

    % if (spec_ratio > cfg.emg_ratio) && ~has_beta_peak
    %     % Rule 1: Pure spectral muscle (broadband gamma/beta plateau, no brain oscillation)
    %     is_safe = true;
    %
    % elseif is_focal && (is_emg_temporal || spec_ratio > 0.18) && ~has_beta_peak
    %     % Rule 2: Hyper-focal non-peripheral EMG / single-lead pop
    %     % Catches electrodes just inside radius 0.40 with steep local gradients
    %     is_safe = true;
    %
    % elseif (pei >= cfg.pei_thresh) && (moi <= cfg.moi_max_thresh)
    %     % Rule 3: Broad outer peripheral rim artifact
    %     is_safe = true;
    %
    % elseif is_emg_temporal && (pei >= 0.25) && ~has_beta_peak
    %     % Rule 4: Tonic neck/jaw tension with moderate peripheral spread
    %     is_safe = true;
    % end

    % Store metrics
    stats(k).ic_num          = k;
    stats(k).pei             = pei;
    stats(k).moi             = moi;
    stats(k).sfi             = sfi;
    stats(k).spec_ratio      = spec_ratio;
    stats(k).has_beta_peak   = has_beta_peak;
    stats(k).kurtosis        = kurt;
    stats(k).hf_roughness    = hf_roughness;
    stats(k).is_temporal_emg = is_emg_temporal;
    stats(k).safe_to_remove  = is_safe;

    safe_to_remove(k) = is_safe;
end

% -------------------------------------------------------------------------
% 4. Console Summary Table
% -------------------------------------------------------------------------
fprintf('\n==================================================================================================\n');
fprintf('  IC Spatial, Spectral & Temporal (Burst/Roughness) Diagnostic Report\n');
fprintf('==================================================================================================\n');
fprintf('  IC# | PEI (Periph) | MOI (Total) | Spec Ratio | Beta Peak? | Kurtosis | Roughness | Action\n');
fprintf('--------------------------------------------------------------------------------------------------\n');
for k = 1:n_ics
    peak_str = 'No'; if stats(k).has_beta_peak, peak_str = 'Yes'; end
    if safe_to_remove(k)
        flag_str = '[REJECT - Muscle/Artefact]';
        % elseif ic_stats(k).moi > cfg.moi_max_thresh
        %     flag_str = '** PROTECT (Motor Strip) **';
    else
        flag_str = 'Retain (Brain / Central)';
    end
    fprintf('  %3d |    %5.1f%%    |    %5.1f%%   |   %6.2f   |    %3s     |  %7.2f |   %6.2f  | %s\n', ...
        k, stats(k).pei * 100, stats(k).moi * 100, stats(k).spec_ratio, ...
        peak_str, stats(k).kurtosis, stats(k).hf_roughness, flag_str);
end
fprintf('==================================================================================================\n\n');

% -------------------------------------------------------------------------
% 5. Visualisation & Saving
% -------------------------------------------------------------------------
if ~isempty(ics_plot)
    n_plots = length(ics_plot);
    Ncol    = 4;
    Nrow    = ceil(n_plots / 2);

    % Dynamic figure dimensions (cm) for saving
    tile_w  = 5.5;   % Width per tile
    tile_h  = 4.5;   % Height per row
    margin  = 2.0;   % Margin for padding/titles
    fig_w   = max(20, Ncol * tile_w);
    fig_h   = max(10, Nrow * tile_h + margin);
    fig_dim = [fig_w, fig_h];

    fh = figure('Color', 'w', 'Name', 'IC Diagnostics', ...
        'Position', [100, 100, 1400, min(1200, 300 * max(1, Nrow))]);
    th = tiledlayout(fh, Nrow, Ncol, 'TileSpacing', 'compact', 'Padding', 'compact');

    try
        colormap(th, brewermap([], '*RdBu'));
    catch
        try
            colormap(fh, brewermap([], '*RdBu'));
        catch
        end
    end

    for i = 1:n_plots
        k = ics_plot(i);

        % Col A: Topoplot using mytopoplot
        ax_topo = nexttile(th);
        set(ax_topo, 'Color', 'w');
        topo_title = sprintf('IC %d Topography', k);
        % mytopoplot(data, mask_channel_mark, title_char, handle_tile, clim_vals)
        mytopoplot(EEG.icawinv(:, k), [], topo_title, ax_topo);

        % Col B: Power Spectrum
        ax_spec = nexttile(th);
        set(ax_spec, 'Color', 'w');
        plot(ax_spec, f, 10*log10(pxx(k, :)), 'k', 'LineWidth', 1.5);
        xlim(ax_spec, [1, 100]); grid(ax_spec, 'on'); hold(ax_spec, 'on');
        xline(ax_spec, 13, '--b'); xline(ax_spec, 30, '--b');
        xline(ax_spec, 55, '--r'); xline(ax_spec, 95, '--r');

        % Plot detected peak(s) if present in the mu/beta range
        f_mu_beta  = f(idx_mu_beta);
        log_spec_k = 10 * log10(pxx(k, idx_mu_beta));
        [pks, locs, ~, p_proms] = findpeaks(log_spec_k);
        sig_peaks = p_proms > cfg.peak_prom;
        if any(sig_peaks)
            plot(ax_spec, f_mu_beta(locs(sig_peaks)), pks(sig_peaks), 'ro', ...
                'MarkerFaceColor', 'r', 'MarkerSize', 5);
        end

        % Recalculate status specifically for the current plot
        if safe_to_remove(k)
            plot_action = '[REJECT]';
            bg_col = [1 0.9 0.9 0.85]; % Light red background for rejected
            % elseif ic_stats(k).moi > cfg.moi_max_thresh
            %      plot_action = '[PROTECT]';
            %      bg_col = [0.9 0.9 1 0.85]; % Light blue background for protected
        else
            plot_action = '[RETAIN]';
            bg_col = [0.9 1 0.9 0.85]; % Light green background for retained
        end

        info_str = sprintf('PEI: %.1f%% | MOI: %.1f%%\nKurt: %.1f | Rough: %.2f\nAction: %s', ...
            stats(k).pei*100, stats(k).moi*100, stats(k).kurtosis, ...
            stats(k).hf_roughness, plot_action);

        text(ax_spec, 0.50, 0.75, info_str, 'Units', 'normalized', 'FontSize', 8, ...
            'BackgroundColor', bg_col, 'EdgeColor', 'k', 'HorizontalAlignment', 'center');
        title(ax_spec, sprintf('IC %d Spectrum', k), 'FontSize', 10, 'FontWeight', 'bold');
    end

    % -------------------------------------------------------------------------
    % Save Figure
    % -------------------------------------------------------------------------
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_ica_diagnostics'], fig_dim);
end

% Log
EEG.ALSUTRECHT.ica.muscle.stats = stats;

end