function DATA = generate_finalplots(DATA, NumberTrials, task, tag, cfg)
% =========================================================================
% GENERATE_FINALPLOTS: BioSemi 128 Quality Control Visualisations
% =========================================================================
% Directly uses central cfg.roi masks for all paradigm-specific overlays.
% =========================================================================

% -------------------------------------------------------------------------
% 1. Config Parsing & Defaults
% -------------------------------------------------------------------------
% % Ensure central ROI definitions exist in cfg (BioSemi 128 defaults)
% if ~isfield(cfg, 'roi'); cfg.roi = struct(); end
%
% if ~isfield(cfg.roi, 'mmn')
%     cfg.roi.mmn = {'C21', 'C20', 'C19', 'C22', 'C23', 'C25', 'C26', 'C12', 'C13', 'D1', 'D2', 'C1', 'C2'};
% end
% if ~isfield(cfg.roi, 'sart_p300')
%     cfg.roi.sart_p300 = {'A19', 'A4', 'A3', 'A2', 'A20', 'A21', 'A5', 'A18', 'A32', 'A31'};
% end
% if ~isfield(cfg.roi, 'rs_alpha')
%     cfg.roi.rs_alpha = {'A23', 'A22', 'A21', 'A24', 'A25', 'A15', 'A16', 'A17', 'A28', 'A29', 'A30', 'A14', 'A27'};
% end
% if ~isfield(cfg.roi, 'motor_left')
%     cfg.roi.motor_left = {'D19', 'D18', 'D14', 'D20', 'D21', 'D12', 'D11', 'D13', 'D17', 'D28', 'D27'};
% end
% if ~isfield(cfg.roi, 'motor_right')
%     cfg.roi.motor_right = {'B22', 'B21', 'B20', 'B23', 'B24', 'B31', 'B30', 'B32', 'B19', 'B18', 'B17'};
% end

% -------------------------------------------------------------------------
% 2. Extract Data & Metadata
% -------------------------------------------------------------------------
fprintf('\n================================\n');
fprintf('Final QC Reports: %s (BioSemi 128 Montage)\n', upper(task));
fprintf('================================\n');

if length(NumberTrials) >= 2
    Nremoved = abs(diff(NumberTrials));
    fprintf('Initial trials   : %d\n', NumberTrials(1));
    fprintf('Stepwise removed : %s (Total: %d)\n', mat2str(Nremoved(:)'), sum(Nremoved));
    fprintf('Remaining trials : %d\n', NumberTrials(end));
end

ALSnr        = DATA.ALSUTRECHT.subject.id;
path_figures = DATA.ALSUTRECHT.subject.figures;
chaneeg      = strcmp({DATA.chanlocs.type}, 'EEG');
chanemg      = strcmp({DATA.chanlocs.type}, 'EMG');
has_emg      = any(chanemg);
chan_labels  = upper({DATA.chanlocs.labels});
eeg_labels   = chan_labels(chaneeg);

% Estimate spectra
[psdspectra, freq] = estimate_power(DATA, 'preproc2');

if tag == 1
    subject = DATA.ALSUTRECHT.subject;
    path_rawpower = fullfile(subject.data, [subject.filename '_rawpower.mat']);
    
    % fh = plot_pre_post_topoplots(path_rawpower, psdspectra', freq);
    [fh, power_diff] = plot_relative_power_diagnostics(path_rawpower, psdspectra', freq);
    
    save_figure(fh, path_figures, sprintf('%s_power_post_final_%s', ALSnr, num2str(tag)), [32 15]);
    DATA.ALSUTRECHT.power_diff = power_diff;
end

% Visual styling
col_gray_trace = [0.25 0.25 0.25 0.15];
col_ga         = [0.05 0.05 0.05];
col_emg        = [0.85 0.33 0.10];

% =========================================================================
% 3. Cognitive ERPs (MMN / SART)
% =========================================================================
if any(strcmpi(task, {'MMN', 'SART'}))
    if strcmpi(task, 'SART')
        if strcmpi(DATA.ALSUTRECHT.SART.type, 'StimulusLocked')
            roi_target      = cfg.roi.sart_p300;
            roi_name        = 'Parietal Pz';
            triggers_unique = [6 3];
            cond_names      = {'Go', 'NoGo'};
            analysis_window = [300 550];
            topo_title      = 'Diff P300 (300-550 ms)';
            size_fig        = [32 9];
        else
            roi_target      = cfg.roi.motor_left;
            roi_name        = 'Central C3';
            triggers_unique = 1;
            cond_names      = {'Go'};
            analysis_window = [0 200];
            topo_title      = 'Diff (0-200 ms)';
            size_fig        = [18 9];
        end

    elseif strcmpi(task, 'MMN')
        roi_target          = cfg.roi.mmn;
        roi_name            = 'Frontal Fz/FCz';
        triggers_unique     = [12 17];
        cond_names          = {'Standard', 'Deviant'};
        analysis_window     = [100 250];
        topo_title          = 'MMN (100-250 ms)';
        size_fig            = [32 9];
    end

    if strcmpi(task, 'SART')
        [med_rt, iqr_rt, rt_pct] = estimate_rt(DATA);
    end

    % Select relevant channels
    roi_mask = ismember(eeg_labels, roi_target);

    % Extract all events and filter for only the condition triggers
    triggers_list = extract_events(DATA);

    % Compute ERP per condition
    num_conds = length(triggers_unique);
    erp_data = cell(num_conds, 1);
    n_trials_cond = zeros(num_conds, 1);

    for i_c = 1:num_conds
        c_mask = (triggers_list == triggers_unique(i_c));
        erp_data{i_c} = mean(DATA.data(chaneeg, :, c_mask), 3);
        n_trials_cond(i_c) = sum(c_mask);
    end

    total_tiles = num_conds + (num_conds == 2) + 1;
    fh = figure('Visible', cfg.figure.visible, 'Color', [1 1 1]);
    tiled_h = tiledlayout(1, total_tiles, 'TileSpacing', 'compact', 'Padding', 'compact');

    all_erp_vals = cell2mat(erp_data);
    y_lims = 1.15 * [min(all_erp_vals(:)), max(all_erp_vals(:))];
    y_lims = [min(y_lims(1), -4), max(y_lims(2), 4)];

    for i_c = 1:num_conds
        ax = nexttile(tiled_h);
        hold(ax, 'on'); box(ax, 'off');

        % Analysis window patch
        patch(ax, [analysis_window(1) analysis_window(2) analysis_window(2) analysis_window(1)], ...
            [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], ...
            [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');

        % -----------------------------------------------------------------
        % Visual RT Indicator (Shaded 25th-75th Percentile + Median Line)
        % -----------------------------------------------------------------
        is_go_stimlocked = strcmpi(task, 'SART') && ...
            strcmpi(cond_names{i_c}, 'Go') && ...
            strcmpi(DATA.ALSUTRECHT.SART.type, 'StimulusLocked') && ...
            ~isnan(med_rt);

        if is_go_stimlocked
            % Shaded patch for Q1 (25%) to Q3 (75%)
            patch(ax, [rt_pct(1) rt_pct(2) rt_pct(2) rt_pct(1)], ...
                [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], ...
                [0.85 0.85 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');

            % Vertical dashed line for Median RT
            p_rt = xline(ax, med_rt, '-.', 'Color', [0.4 0.2 0.7], 'LineWidth', 1.5, 'DisplayName', 'RT');
        end

        % Grid / Zero lines
        xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        yline(ax, 0, '-',  'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');

        % ERP Waveforms
        plot(ax, DATA.times, erp_data{i_c}, 'Color', col_gray_trace, 'LineWidth', 0.6, 'HandleVisibility', 'off');
        % p_ga  = plot(ax, DATA.times, mean(erp_data{i_c}, 1), 'Color', col_ga, 'LineWidth', 1.8, 'DisplayName', 'Grand Avg (128)');
        p_roi = plot(ax, DATA.times, mean(erp_data{i_c}(roi_mask, :), 1), 'Color', [0.00 0.45 0.74], 'LineWidth', 2.2, 'DisplayName', roi_name);

        xlim(ax, [DATA.times(1), DATA.times(end)]);
        ylim(ax, y_lims);
        xlabel(ax, 'Time (ms)'); ylabel(ax, 'Amplitude (\muV)');

        % Title displaying Median and IQR
        if is_go_stimlocked
            title_str = sprintf('%s: %s (N=%d | RT=%.0f ms, IQR=%.0f)', ...
                ALSnr, cond_names{i_c}, n_trials_cond(i_c), med_rt, iqr_rt);
        else
            title_str = sprintf('%s: %s (N=%d)', ALSnr, cond_names{i_c}, n_trials_cond(i_c));
        end
        title(ax, title_str, 'Interpreter', 'none');

        if is_go_stimlocked
            legend(ax, [p_roi, p_rt], 'Location', 'southwest', 'FontSize', 6.5);
        else
            legend(ax, p_roi, 'Location', 'southwest', 'FontSize', 6.5);
        end

        pbaspect(ax, [1.4 1 1]);
    end

    if num_conds == 2
        diff_wave = erp_data{2} - erp_data{1};
        ax = nexttile(tiled_h);
        hold(ax, 'on'); box(ax, 'off');

        % Analysis window patch
        patch(ax, [analysis_window(1) analysis_window(2) analysis_window(2) analysis_window(1)], ...
            [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], ...
            [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');

        % -----------------------------------------------------------------
        % Visual RT Indicator (Shaded SD window + Mean line)
        % -----------------------------------------------------------------
        is_sart_stimlocked = strcmpi(task, 'SART') && ...
            strcmpi(DATA.ALSUTRECHT.SART.type, 'StimulusLocked') && ...
            ~isnan(med_rt);

        if is_sart_stimlocked
            % Shaded patch for Q1 to Q3
            patch(ax, [rt_pct(1) rt_pct(2) rt_pct(2) rt_pct(1)], ...
                [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], ...
                [0.85 0.85 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');

            % Vertical dashed line for Median RT
            p_rt_diff = xline(ax, med_rt, '-.', 'Color', [0.4 0.2 0.7], 'LineWidth', 1.5, 'DisplayName', 'RT');
        end

        % Grid / Zero lines
        xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        yline(ax, 0, '-',  'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');

        plot(ax, DATA.times, diff_wave, 'Color', col_gray_trace, 'LineWidth', 0.6, 'HandleVisibility', 'off');
        % p_ga_diff  = plot(ax, DATA.times, mean(diff_wave, 1), 'Color', col_ga, 'LineWidth', 1.8, 'DisplayName', 'Grand Avg (128)');
        p_roi_diff = plot(ax, DATA.times, mean(diff_wave(roi_mask, :), 1), 'Color', [0.85 0.10 0.10], 'LineWidth', 2.2, 'DisplayName', roi_name);

        xlim(ax, [DATA.times(1), DATA.times(end)]);
        ylim(ax, y_lims);
        xlabel(ax, 'Time (ms)'); ylabel(ax, 'Amplitude (\muV)');
        title(ax, sprintf('%s: Diff (%s - %s)', ALSnr, cond_names{2}, cond_names{1}), 'Interpreter', 'none');

        if is_sart_stimlocked
            legend(ax, [p_roi_diff, p_rt_diff], 'Location', 'southwest', 'FontSize', 6.5);
        else
            legend(ax, p_roi_diff, 'Location', 'southwest', 'FontSize', 6.5);
        end

        pbaspect(ax, [1.4 1 1]);
    end

    ax_topo = nexttile(tiled_h);
    if num_conds == 2; target_topo = diff_wave; else; target_topo = erp_data{1}; end

    time_mask = (DATA.times >= analysis_window(1)) & (DATA.times <= analysis_window(2));
    topo_vals = mean(target_topo(:, time_mask), 2);

    max_c = max(abs(topo_vals));
    if max_c == 0; max_c = 1; end
    mytopoplot(topo_vals, [], topo_title, ax_topo, [-max_c max_c]);

    cb = colorbar;
    cb.Location = 'southoutside';
    cb.Label.String = 'uV';

    save_figure(fh, path_figures, sprintf('%s_erp_final_%s', ALSnr, num2str(tag)), size_fig);

    % =========================================================================
    % 4. Resting-State (RS) PSD
    % =========================================================================
elseif strcmpi(task, 'RS')

    psd_eeg = psdspectra(:, chaneeg);
    psd_db  = 10 * log10(psd_eeg + eps);
    alpha_roi_mask = ismember(eeg_labels, cfg.roi.rs_alpha);

    fh = figure('Visible', cfg.figure.visible, 'Color', [1 1 1]);
    tiled_h = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tiled_h);
    hold(ax1, 'on'); box(ax1, 'off');

    y_range = [min(psd_db(:)) - 2, max(psd_db(:)) + 2];
    y_range(1) = max(y_range(1), -40);
    patch(ax1, [8 12 12 8], [y_range(1) y_range(1) y_range(2) y_range(2)], ...
        [0.92 0.96 0.92], 'EdgeColor', 'none', 'DisplayName', 'Alpha (8-12 Hz)');
    patch(ax1, [15 30 30 15], [y_range(1) y_range(1) y_range(2) y_range(2)], ...
        [0.93 0.93 0.98], 'EdgeColor', 'none', 'DisplayName', 'Beta (15-30 Hz)');

    plot(ax1, freq, psd_db, 'Color', col_gray_trace, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot(ax1, freq, mean(psd_db, 2), 'Color', col_ga, 'LineWidth', 2.2, 'DisplayName', 'Grand Avg');

    if any(alpha_roi_mask)
        plot(ax1, freq, mean(psd_db(:, alpha_roi_mask), 2), ...
            'Color', [0.10 0.65 0.10], 'LineWidth', 2.2, 'DisplayName', 'Occipital Alpha');
    end

    xlim(ax1, [1 70]); ylim(ax1, y_range);
    xlabel(ax1, 'Frequency (Hz)'); ylabel(ax1, 'Power Spectral Density (dB/Hz)');
    title(ax1, sprintf('%s: Resting-State PSD', ALSnr), 'Interpreter', 'none');
    legend(ax1, 'Location', 'southeast', 'FontSize', 6.5);
    pbaspect(ax1, [1.5 1 1]);

    ax2 = nexttile(tiled_h);
    hold(ax2, 'on'); box(ax2, 'off');

    freq_log = log10(freq(freq >= 1));
    psd_log  = log10(psd_eeg(freq >= 1, :) + eps);

    plot(ax2, freq_log, psd_log, 'Color', col_gray_trace, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot(ax2, freq_log, mean(psd_log, 2), 'Color', [0.1 0.5 0.1], 'LineWidth', 2.2, 'DisplayName', 'Mean 1/f Spectrum');

    ticks_hz = [1 2 4 8 12 20 30 45 70];
    set(ax2, 'XTick', log10(ticks_hz), 'XTickLabel', arrayfun(@num2str, ticks_hz, 'UniformOutput', false));
    xlim(ax2, log10([1 70]));
    xlabel(ax2, 'Frequency (Hz, log scale)'); ylabel(ax2, 'log_{10}(Power)');
    title(ax2, sprintf('%s: Aperiodic Decay (1/f Check)', ALSnr), 'Interpreter', 'none');
    legend(ax2, 'Location', 'southwest', 'FontSize', 6.5);
    pbaspect(ax2, [1.5 1 1]);

    save_figure(fh, path_figures, sprintf('%s_psd_final_%s', ALSnr, num2str(tag)), [24 9]);

    % =========================================================================
    % 5. Motor Task (MT) PSD (Consuming cfg.roi.motor_left / right)
    % =========================================================================
elseif strcmpi(task, 'MT')

    psd_eeg    = psdspectra(:, chaneeg);
    psd_eeg_db = 10 * log10(psd_eeg + eps);

    % Use central cfg.roi masks
    mask_left  = ismember(eeg_labels, cfg.roi.motor_left);
    mask_right = ismember(eeg_labels, cfg.roi.motor_right);

    fh = figure('Visible', cfg.figure.visible, 'Color', [1 1 1]);
    tiled_h = tiledlayout(1, 2 + has_emg, 'TileSpacing', 'compact', 'Padding', 'compact');

    % A. Motor EEG Spectrum
    ax1 = nexttile(tiled_h);
    hold(ax1, 'on'); box(ax1, 'off');

    y_lim_eeg = [min(psd_eeg_db(:)) - 2, max(psd_eeg_db(:)) + 2];
    y_lim_eeg(1) = max(y_lim_eeg(1), -40);

    % Set transparency and y-coordinates
    patch_y = [y_lim_eeg(1) y_lim_eeg(1) y_lim_eeg(2) y_lim_eeg(2)];
    alpha_val = 0.20;

    % Mu (8-13 Hz) - Soft Blue
    patch(ax1, [8 13 13 8], patch_y, [0.12 0.47 0.71], ...
        'EdgeColor', 'none', 'FaceAlpha', alpha_val, 'DisplayName', 'Mu (8-13 Hz)');
    % Beta (13-30 Hz) - Soft Green
    patch(ax1, [13 30 30 13], patch_y, [0.20 0.63 0.35], ...
        'EdgeColor', 'none', 'FaceAlpha', alpha_val, 'DisplayName', 'Beta (13-30 Hz)');
    % Low Gamma (30-45 Hz) - Soft Coral / Amber
    patch(ax1, [30 45 45 30], patch_y, [0.85 0.37 0.10], ...
        'EdgeColor', 'none', 'FaceAlpha', alpha_val, 'DisplayName', 'Low gamma (30-45 Hz)');


    plot(ax1, freq, psd_eeg_db, 'Color', col_gray_trace, 'LineWidth', 0.6, 'HandleVisibility', 'off');
    % plot(ax1, freq, mean(psd_eeg_db, 2), 'Color', col_ga, 'LineWidth', 1.8, 'DisplayName', 'Grand Avg (128)');

    if any(mask_left)
        plot(ax1, freq, mean(psd_eeg_db(:, mask_left), 2), ...
            'Color', [0.00 0.45 0.74], 'LineWidth', 2.2, 'DisplayName', 'Left Motor (C3)');
    end

    if any(mask_right)
        plot(ax1, freq, mean(psd_eeg_db(:, mask_right), 2), ...
            'Color', [0.85 0.15 0.15], 'LineWidth', 2.2, 'DisplayName', 'Right Motor (C4)');
    end

    xlim(ax1, [1 48]); ylim(ax1, y_lim_eeg);
    xlabel(ax1, 'Frequency (Hz)'); ylabel(ax1, 'EEG Power (dB/Hz)');
    title(ax1, sprintf('%s: Sensorimotor PSD (Left vs Right)', ALSnr), 'Interpreter', 'none');
    legend(ax1, 'Location', 'southwest', 'FontSize', 6.5);
    pbaspect(ax1, [1.4 1 1]);

    % B. Surface EMG Spectrum
    if has_emg
        psd_emg    = psdspectra(:, chanemg);
        psd_emg_db = 10 * log10(psd_emg + eps);

        ax2 = nexttile(tiled_h);
        hold(ax2, 'on'); box(ax2, 'off');

        y_lim_emg = [min(psd_emg_db(:)) - 2, max(psd_emg_db(:)) + 2];
        patch(ax2, [15 45 45 15], [y_lim_emg(1) y_lim_emg(1) y_lim_emg(2) y_lim_emg(2)], ...
            [0.95 0.95 0.95], 'EdgeColor', 'none', 'DisplayName', 'CMC Drive Window');

        plot(ax2, freq, psd_emg_db, 'Color', [0.85 0.4 0.2 0.4], 'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(ax2, freq, mean(psd_emg_db, 2), 'Color', col_emg, 'LineWidth', 2.2, 'DisplayName', 'Mean Raw EMG');

        xlim(ax2, [5 60]); ylim(ax2, y_lim_emg);
        xlabel(ax2, 'Frequency (Hz)'); ylabel(ax2, 'EMG Power (dB/Hz)');
        title(ax2, sprintf('%s: Surface EMG Spectrum', ALSnr), 'Interpreter', 'none');
        legend(ax2, 'Location', 'southwest', 'FontSize', 6.5);
        pbaspect(ax2, [1.4 1 1]);
    end

    % C. Topography of Beta Band (15 to 30 Hz) Power
    ax_topo = nexttile(tiled_h);
    beta_mask = (freq >= 15) & (freq <= 30);
    beta_power = mean(psd_eeg(beta_mask, :), 1);
    beta_power_norm = beta_power ./ sum(psd_eeg, 1);

    mytopoplot(beta_power_norm, [], 'Rel. Beta Power (15-30 Hz)', ax_topo);

    cb = colorbar;
    cb.Location = 'southoutside';
    cb.Label.String = '%';

    save_figure(fh, path_figures, sprintf('%s_psd_final_%s', ALSnr, num2str(tag)), [18 9]);
end

% Diagnostic IAF Peak
plot_iaf(DATA, freq, mean(psdspectra(:, chaneeg), 2), cfg.figure.visible);

end

function plot_iaf(EEG, freq1, psdspectra1, opt_visible)
% =========================================================================
% PLOT_IAF: Individual Alpha Frequency (IAF / PAF) Diagnostic Report
% =========================================================================
% Overlays the whole-cap average power spectrum with the dedicated IAF channel
% spectrum, canonical frequency bands, and the peak alpha frequency marker.
% =========================================================================

% -------------------------------------------------------------------------
% 1. Data Extraction & Validation
% -------------------------------------------------------------------------
psdspectra2 = mean(EEG.ALSUTRECHT.pSpec.sums.muSpec, 2);
freq2       = EEG.ALSUTRECHT.pSpec.sums.freq;
paf         = EEG.ALSUTRECHT.pSpec.sums.paf;

if isnan(paf)
    psdspectra2 = [];
    freq2       = [];
end

% Ensure column vectors and compute average across channels if 2D
if size(psdspectra1, 2) > 1
    psdspectra1 = mean(psdspectra1, 2);
end
psdspectra1 = psdspectra1(:);
freq1       = freq1(:);

% Normalise spectra to unit maximum (1.0) within 1-48 Hz for visual comparison
in_range1 = freq1 >= 1 & freq1 <= 48;
max_p1    = max(psdspectra1(in_range1));
if isempty(max_p1) || max_p1 == 0; max_p1 = max(psdspectra1) + eps; end
psd_norm1 = psdspectra1 ./ max_p1;

if ~isempty(psdspectra2)
    psdspectra2 = psdspectra2(:);
    freq2       = freq2(:);
    in_range2   = freq2 >= 1 & freq2 <= 48;
    max_p2      = max(psdspectra2(in_range2));
    if isempty(max_p2) || max_p2 == 0; max_p2 = max(psdspectra2) + eps; end
    psd_norm2   = psdspectra2 ./ max_p2;
else
    psd_norm2   = [];
end

% -------------------------------------------------------------------------
% 2. Canonical Frequency Band Definitions (Contiguous Bands)
% -------------------------------------------------------------------------
bands = struct();
bands(1).name  = 'Delta'; bands(1).range = [1 4];   bands(1).col = [0.88 0.88 0.88]; % Neutral grey
bands(2).name  = 'Theta'; bands(2).range = [4 8];   bands(2).col = [0.82 0.89 0.95]; % Soft blue
bands(3).name  = 'Alpha'; bands(3).range = [8 13];  bands(3).col = [0.80 0.93 0.80]; % Target green
bands(4).name  = 'Beta';  bands(4).range = [13 30]; bands(4).col = [0.96 0.89 0.80]; % Soft amber
bands(5).name  = 'Gamma'; bands(5).range = [30 48]; bands(5).col = [0.93 0.84 0.91]; % Soft mauve

% -------------------------------------------------------------------------
% 3. Plotting
% -------------------------------------------------------------------------
fh = figure('Visible', opt_visible, 'Color', [1 1 1]);
ax = axes(fh);
hold(ax, 'on'); box(ax, 'off');

y_lim = [0 1.15];

% A. Background Frequency Band Patches
h_bands = gobjects(length(bands), 1);
for i_b = 1:length(bands)
    bx = [bands(i_b).range(1), bands(i_b).range(2), bands(i_b).range(2), bands(i_b).range(1)];
    by = [y_lim(1), y_lim(1), y_lim(2), y_lim(2)];
    h_bands(i_b) = patch(ax, bx, by, bands(i_b).col, ...
        'FaceAlpha', 0.45, 'EdgeColor', 'none', 'DisplayName', bands(i_b).name);
end

% B. Whole-Cap Grand Average Spectrum
h_avg = plot(ax, freq1, psd_norm1, 'Color', [0.25 0.25 0.25], 'LineWidth', 1.8, ...
    'DisplayName', 'Whole-Cap Average');

% C. IAF Dedicated Channel Spectrum
h_iaf = [];
if ~isempty(psd_norm2)
    h_iaf = plot(ax, freq2, psd_norm2, 'Color', [0.00 0.45 0.74], 'LineWidth', 2.2, ...
        'DisplayName', 'IAF Channel Spectrum');
end

% D. Peak Alpha Frequency Marker & Annotation
h_paf = [];
if ~isnan(paf) && paf >= 5 && paf <= 15
    if ~isempty(psd_norm2)
        [~, idx_p] = min(abs(freq2 - paf));
        paf_y = psd_norm2(idx_p);
    else
        [~, idx_p] = min(abs(freq1 - paf));
        paf_y = psd_norm1(idx_p);
    end

    % Vertical dashed drop line
    xline(ax, paf, ':', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % PAF marker dot
    h_paf = plot(ax, paf, paf_y, 'o', 'MarkerSize', 7, ...
        'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', [0.40 0 0], ...
        'LineWidth', 1.2, 'DisplayName', sprintf('PAF Peak (%.2f Hz)', paf));

    % Text callout badge
    text(ax, paf + 0.8, min(paf_y + 0.05, 1.05), sprintf('PAF = %.2f Hz', paf), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.70 0.05 0.05], ...
        'BackgroundColor', [1 1 1 0.85], 'Margin', 2);
end

% -------------------------------------------------------------------------
% 4. Axes, Titles & Legend Formatting
% -------------------------------------------------------------------------
xlim(ax, [1 48]);
ylim(ax, y_lim);
set(ax, 'XTick', [1 4 8 13 30 48]);
xlabel(ax, 'Frequency (Hz)');
ylabel(ax, 'Normalised Power (a.u.)');
pbaspect(ax, [1.618 1 1]);

if ~isnan(paf)
    title(ax, sprintf('Individual Alpha Frequency: %.2f Hz', paf));
else
    title(ax, 'Individual Alpha Frequency: Not Detected');
end

% Build clean legend using explicit graphic handles
legend_handles = h_avg;
if ~isempty(h_iaf); legend_handles = [legend_handles; h_iaf]; end
if ~isempty(h_paf); legend_handles = [legend_handles; h_paf]; end
legend_handles = [legend_handles; h_bands];

legend(ax, legend_handles, 'Location', 'eastoutside', 'FontSize', 6.5);

% -------------------------------------------------------------------------
% 5. Save Figure
% -------------------------------------------------------------------------
save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_pspectra_iaf'], [20 8]);

end


function triggers_list = extract_events(DATA)
% -------------------------------------------------------------------------
% Extract the exact time-locking trigger for each epoched trial
% -------------------------------------------------------------------------
num_trials = size(DATA.data, 3);
triggers_list = zeros(1, num_trials);

for i_tr = 1:num_trials
    % Extract latencies for all events within this specific epoch
    ep_lats = DATA.epoch(i_tr).eventlatency;
    if iscell(ep_lats)
        ep_lats = cell2mat(ep_lats);
    end

    % The time-locking stimulus/response that created the epoch is at 0 ms
    [~, zero_idx] = min(abs(ep_lats));

    % Extract the trigger value at latency 0
    if isfield(DATA.epoch(i_tr), 'eventedftype') && ~isempty(DATA.epoch(i_tr).eventedftype)
        raw_val = DATA.epoch(i_tr).eventedftype;
    else
        raw_val = DATA.epoch(i_tr).eventtype;
    end

    if iscell(raw_val)
        val = raw_val{zero_idx};
    else
        val = raw_val(zero_idx);
    end

    % Convert strings (e.g., 'condition 1' or '1') to numeric
    if ischar(val) || isstring(val)
        val = str2double(regexprep(string(val), '\D', ''));
    end

    triggers_list(i_tr) = val;
end

% Verify every trial has been matched
assert(length(triggers_list) == size(DATA.data, 3), ...
    'Trial count mismatch in DATA.epoch.');
end

function [med_rt, iqr_rt, rt_pct] = estimate_rt(DATA)
% -------------------------------------------------------------------------
% Estimate SART Response Time (Median & IQR / Percentiles)
% -------------------------------------------------------------------------
med_rt = NaN;
iqr_rt = NaN;
rt_pct = [NaN, NaN];

if isfield(DATA, 'epoch') && ~isempty(DATA.epoch)
    rt_min = 170; % Minimum plausible RT in ms
    rt_max = 550; % Maximum plausible RT in ms
    go_rts = [];

    for i_ep = 1:length(DATA.epoch)
        if isfield(DATA.epoch(i_ep), 'eventedftype') && ~isempty(DATA.epoch(i_ep).eventedftype)
            raw_ev = DATA.epoch(i_ep).eventedftype;
        else
            raw_ev = DATA.epoch(i_ep).eventtype;
        end

        if iscell(raw_ev)
            if isnumeric(raw_ev{1})
                ep_types = cell2mat(raw_ev);
            else
                ep_types = cellfun(@str2double, raw_ev);
            end
        elseif ischar(raw_ev) || isstring(raw_ev)
            ep_types = str2double(string(raw_ev));
        else
            ep_types = raw_ev;
        end

        ep_lats = DATA.epoch(i_ep).eventlatency;
        if iscell(ep_lats)
            ep_lats = cell2mat(ep_lats);
        end

        % Find button press (trigger 1) following Go stimulus (trigger 6)
        if any(ep_types == 6)
            resp_idx = find((ep_types == 1) & (ep_lats > 0), 1, 'first');
            if ~isempty(resp_idx)
                resp_lat = ep_lats(resp_idx);
                if (resp_lat >= rt_min) && (resp_lat <= rt_max)
                    go_rts(end+1) = resp_lat; %#ok<AGROW>
                end
            end
        end
    end

    if ~isempty(go_rts)
        med_rt = median(go_rts);
        iqr_rt = iqr(go_rts);
        rt_pct = prctile(go_rts, [25, 75]);
    end
end
end


function [fh, diff_db] = plot_relative_power_diagnostics(rawpower_path, psd_post, freq_post)
% PLOT_RELATIVE_POWER_DIAGNOSTICS Compares pre- and post-cleaning spectra
% across 3 rows:
%   Row 1: Relative Power Before (% of 1-45 Hz broadband)
%   Row 2: Relative Power After (% of 1-45 Hz broadband)
%   Row 3: Absolute Attenuation in Decibels (10*log10(Post / Pre))
%          with strictly zero-centred symmetric limits [-max, +max].

% -------------------------------------------------------------------------
% 1. Load Pre-Cleaning Data
% -------------------------------------------------------------------------
assert(exist(rawpower_path, 'file') == 2, 'File not found: %s', rawpower_path);
raw_dat = load(rawpower_path);

psd_pre   = raw_dat.psd_pre;
freq_pre  = raw_dat.freq_pre(:);
freq_post = freq_post(:);

n_chans = size(psd_pre, 1);
assert(size(psd_post, 1) == n_chans, ...
    'Channel mismatch: psd_pre has %d channels, psd_post has %d.', n_chans, size(psd_post, 1));

% -------------------------------------------------------------------------
% 2. Define Frequency Bands and Total Power Window
% -------------------------------------------------------------------------
f_broadband = [1 45];
idx_bb_pre  = freq_pre >= f_broadband(1)  & freq_pre <= f_broadband(2);
idx_bb_post = freq_post >= f_broadband(1) & freq_post <= f_broadband(2);

tot_pwr_pre  = trapz(freq_pre(idx_bb_pre),   psd_pre(:, idx_bb_pre),   2);
tot_pwr_post = trapz(freq_post(idx_bb_post), psd_post(:, idx_bb_post), 2);

bands = { ...
    'Delta', [1 4]; ...
    'Theta', [4 8]; ...
    'Alpha', [8 13]; ...
    'Beta',  [13 30]; ...
    'Gamma', [30 45] ...
    };
n_bands = size(bands, 1);

abs_pwr_pre  = zeros(n_chans, n_bands);
abs_pwr_post = zeros(n_chans, n_bands);
rel_pwr_pre  = zeros(n_chans, n_bands);
rel_pwr_post = zeros(n_chans, n_bands);

for i_b = 1:n_bands
    f_range = bands{i_b, 2};
    idx_pre  = freq_pre >= f_range(1)  & freq_pre <= f_range(2);
    idx_post = freq_post >= f_range(1) & freq_post <= f_range(2);

    % Absolute power integration (\muV^2)
    abs_pwr_pre(:, i_b)  = trapz(freq_pre(idx_pre),   psd_pre(:, idx_pre),   2);
    abs_pwr_post(:, i_b) = trapz(freq_post(idx_post), psd_post(:, idx_post), 2);

    % Relative power (% of 1-45 Hz total power)
    rel_pwr_pre(:, i_b)  = (abs_pwr_pre(:, i_b)  ./ tot_pwr_pre)  * 100;
    rel_pwr_post(:, i_b) = (abs_pwr_post(:, i_b) ./ tot_pwr_post) * 100;
end

% Absolute attenuation in Decibels (Post vs Pre)
diff_db = 10 * log10(abs_pwr_post ./ abs_pwr_pre);

% Absolute retention percentage (Post / Pre * 100%)
retention_pct = (abs_pwr_post ./ abs_pwr_pre) * 100;

% -------------------------------------------------------------------------
% 3. Console Integrity Diagnostic
% -------------------------------------------------------------------------
fprintf('\n--- Power Retention Audit (Median %% Across All Channels) ---\n');
for i_b = 1:n_bands
    fprintf('  %-6s (%2d-%2d Hz): Median Retention = %5.1f%% (IQR: %4.1f - %4.1f%%) | Max dB = %+5.2f dB\n', ...
        bands{i_b, 1}, bands{i_b, 2}(1), bands{i_b, 2}(2), ...
        median(retention_pct(:, i_b), 'omitnan'), ...
        prctile(retention_pct(:, i_b), 25), ...
        prctile(retention_pct(:, i_b), 75), ...
        max(diff_db(:, i_b), [], 'omitnan'));
end
fprintf('------------------------------------------------------------\n\n');

% -------------------------------------------------------------------------
% 4. Tiled Topography Plot (3 Rows x 5 Canonical Bands)
% -------------------------------------------------------------------------
fh = figure('Color', 'w', 'Position', [50 50 1600 850]);
t = tiledlayout(3, n_bands, 'TileSpacing', 'compact', 'Padding', 'normal');

for i_b = 1:n_bands
    band_name  = bands{i_b, 1};
    band_range = bands{i_b, 2};

    r_pre  = rel_pwr_pre(:, i_b);
    r_post = rel_pwr_post(:, i_b);
    r_db   = diff_db(:, i_b);

    med_ret = median(retention_pct(:, i_b), 'omitnan');
    q25_ret = prctile(retention_pct(:, i_b), 25);
    q75_ret = prctile(retention_pct(:, i_b), 75);

    % % Shared colour limits for Pre and Post rows
    % c_max = max([r_pre; r_post], [], 'omitnan');
    % c_min = min([r_pre; r_post], [], 'omitnan');
    % if c_min == c_max, c_max = c_min + 1; end
    % clim_rel = [0, c_max];

    % Shared colour limits for Pre and Post rows
    c_max_1 = max(r_pre, [], 'omitnan');
    c_max_2 = max(r_post, [], 'omitnan');
    clim_rel_1 = [0, c_max_1];
    clim_rel_2 = [0, c_max_2];

    % Strictly symmetric zero-centred limits for Difference row
    finite_db = r_db(isfinite(r_db));
    if isempty(finite_db)
        max_shift = 1;
    else
        max_shift = max(abs(finite_db), [], 'omitnan');
        if max_shift == 0, max_shift = 1; end
    end
    clim_shift = [-max_shift, max_shift];

    % Row 1: Relative Power Before
    ax1 = nexttile(i_b);
    mytopoplot(r_pre, [], '', ax1, clim_rel_1);
    colormap(ax1, brewermap([], 'Reds'));
    hcb1 = colorbar(ax1);
    hcb1.Title.String = '%';

    title(ax1, sprintf('%s\n(%d–%d Hz)', band_name, band_range(1), band_range(2)), ...
        'FontSize', 11, 'FontWeight', 'bold');

    if i_b == 1
        ylabel(ax1, 'Before (Rel %)', 'FontWeight', 'bold', 'FontSize', 12, 'Visible', 'on');
    end

    % Row 2: Relative Power After
    ax2 = nexttile(i_b + n_bands);
    mytopoplot(r_post, [], '', ax2, clim_rel_2);
    colormap(ax2, brewermap([], 'Reds'));
    hcb2 = colorbar(ax2);
    hcb2.Title.String = '%';

    if i_b == 1
        ylabel(ax2, 'After (Rel %)', 'FontWeight', 'bold', 'FontSize', 12, 'Visible', 'on');
    end

    % Row 3: Absolute Attenuation in Decibels
    ax3 = nexttile(i_b + 2 * n_bands);
    mytopoplot(r_db, false(size(r_db)), '', ax3, clim_shift);
    colormap(ax3, brewermap([], '*RdBu'));

    clim(ax3, clim_shift);
    set(ax3, 'CLim', clim_shift);

    hcb3 = colorbar(ax3);
    hcb3.Title.String = 'dB';

    title(ax3, sprintf('Ret: %.0f%%\n(IQR: %.0f–%.0f%%)', med_ret, q25_ret, q75_ret), ...
        'FontSize', 10, 'FontWeight', 'bold');

    if i_b == 1
        ylabel(ax3, 'Difference (dB)', 'FontWeight', 'bold', 'FontSize', 12, 'Visible', 'on');
    end
end

end



% function generate_finalplots(DATA, NumberTrials, task, tag, opt_visible)
%
% fprintf('\n================================\n');
% fprintf('Final reports\n');
% fprintf('================================\n');
%
% % Report
% % Nremoved = [NumberTrials(1)-NumberTrials(2), NumberTrials(2)-NumberTrials(3)];
% Nremoved = abs(diff(NumberTrials));
% fprintf('Inital trials:    %d\n',NumberTrials(1));
% fprintf('Removed trials:   %d + %d + %d + %d = %d\n', Nremoved, sum(Nremoved));
% fprintf('Remaining trials: %d\n', NumberTrials(end));
%
% % Estimate power sepctra
% [psdspectra, freq, chaneeg, chanemg] = estimate_power(DATA, 'preproc2');
% freqlog = log10(freq);
%
% % Freq labels
% freqTicks     = [2:2:10 15 20 30 50 100];
% freqTicksLog  = log10(freqTicks);
% freqTicksCell = arrayfun(@num2str, freqTicks, 'UniformOutput', false);
%
% % Extract
% ALSnr = DATA.ALSUTRECHT.subject.id;
% path_figures = DATA.ALSUTRECHT.subject.figures;
%
% % Check if there are EMG signals in the data struct
% has_emg = any(strcmpi({DATA.chanlocs.type}, 'EMG'));
%
% % =============================
% % Plot 1
% % =============================
% % Plots differe per task
% if strcmpi(task, 'MMN') || strcmpi(task, 'SART')
%     % ERP
%     triggers_list = [DATA.event.edftype];
%     if strcmpi(task, 'SART')
%         if strcmpi(DATA.ALSUTRECHT.SART.type, 'StimulusLocked')
%             triggers_unique = [3 6];
%             triggers_mask = triggers_list==triggers_unique(1) | triggers_list==triggers_unique(2);
%             num_tiles = 3;
%         elseif strcmpi(DATA.ALSUTRECHT.SART.type, 'ResponseLocked')
%             triggers_unique = 1;
%             % maskTrig = maskCond==condTrig;
%             num_tiles = 1;
%         end
%         minClim = [-8 8];
%
%     elseif strcmpi(task, 'MMN')
%         triggers_unique = [17 12];
%         triggers_mask = triggers_list==triggers_unique(1) | triggers_list==triggers_unique(2);
%         minClim = [-3 3];
%         num_tiles = 4;
%     end
%
%     % Number of conditions
%     num_triggers = length(triggers_unique);
%
%     if num_triggers == 2
%         triggers_list = triggers_list(triggers_mask);
%         assert(DATA.trials == sum(triggers_mask));
%         % condTrig = unique(maskCond);
%
%         % ERP 1/2
%         data_tmp_1 = mean(DATA.data(chaneeg, :, triggers_list==triggers_unique(1)), 3);
%         data_tmp_2 = mean(DATA.data(chaneeg, :, triggers_list==triggers_unique(2)), 3);
%
%         NumberTrials = NaN(num_triggers, 1);
%         NumberTrials(1) = sum(triggers_list == triggers_unique(1));
%         NumberTrials(2) = sum(triggers_list == triggers_unique(2));
%
%     elseif num_triggers == 1
%         % ERP 1
%         data_tmp_1 = mean(DATA.data(chaneeg, :, :), 3);
%         NumberTrials = DATA.trials;
%
%     end
%
%     % Select only EEG
%     dataCmap = brewermap(128, 'PRGn');
%     % dataCmap = brewermap(128, 'BrBG');
%
%     fh = figure('Visible', opt_visible);
%     th = tiledlayout(1, num_tiles);
%     th.TileSpacing = 'compact'; th.Padding = 'compact';
%
%     % ERP 1 & 2
%     for i = 1:num_triggers
%         if i == 1
%             data_plot = data_tmp_1;
%         else
%             data_plot = data_tmp_2;
%         end
%
%         dataClim = 1.1 * [min(data_plot(:)), max(data_plot(:))];
%         dataClim(1) = floor(dataClim(1));
%         dataClim(2) = ceil(dataClim(2));
%
%         dataClim(1) = min(dataClim(1), minClim(1));
%         dataClim(2) = max(dataClim(2), minClim(1));
%
%         th = nexttile;
%         hold on; box off;
%         xline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
%         yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
%         % 1. Plot all channels slightly thinner to reduce clutter
%         plot(DATA.times, data_plot, 'LineWidth', 0.8);
%         % 2. Plot the Grand Average (mean across all channels) on top in thick black
%         plot(DATA.times, mean(data_plot, 1), 'k', 'LineWidth', 2.5);
%
%         colororder(th,dataCmap); xlim(DATA.times([1 end])); ylim(dataClim);
%         title({ALSnr, [num2str(sum(chaneeg)) ' EEG, trig ' num2str(triggers_unique(i)) ', ' num2str(NumberTrials(i)) ' trials']});
%         pbaspect([1.618 1 1]); xlabel('Time (ms)'); ylabel('Amplitude (uV)');
%     end
%
%     % ERP difference
%     if num_triggers == 2
%         data_plot = data_tmp_1 - data_tmp_2;
%
%         dataClim = 1.1 * [min(data_plot(:)), max(data_plot(:))];
%         dataClim(1) = floor(dataClim(1));
%         dataClim(2) = ceil(dataClim(2));
%
%         dataClim(1) = min(dataClim(1), minClim(1));
%         dataClim(2) = max(dataClim(2), minClim(1));
%
%         th = nexttile;
%         hold on; box off;
%         xline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
%         yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
%         % 1. Plot all channels slightly thinner to reduce clutter
%         plot(DATA.times, data_plot, 'LineWidth', 0.8);
%         % 2. Plot the Grand Average (mean across all channels) on top in thick black
%         plot(DATA.times, mean(data_plot, 1), 'k', 'LineWidth', 2.5);
%
%         colororder(th,dataCmap); xlim(DATA.times([1 end])); ylim(dataClim);
%         title({ALSnr, [num2str(sum(chaneeg)) ' EEG, ERP difference (' num2str(triggers_unique(1)) '-' num2str(triggers_unique(2)) ')']});
%         pbaspect([1.618 1 1]); xlabel('Time (ms)'); ylabel('Amplitude (uV)');
%     end
%
%     % Add MMN topolot
%     if strcmpi(task, 'MMN')
%         [bl, al] = butter(2, 20/(DATA.srate/2), 'low'); assert(isstable(bl, al));
%         data_plot = filtfilt(bl, al, data_plot')';
%
%         data_plot_1 = mean(data_plot(:, DATA.times>100 & DATA.times<250), 2);
%         data_plot_0 = mean(data_plot(:, DATA.times<0), 2);
%
%         mytopoplot(data_plot_1 - data_plot_0, [], 'Lowpass filtered MMN (<20 Hz, 100-250 ms)', nexttile, 0.5*[-1 1]);
%     end
%
%     % plotX=30; plotY=8;
%     % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
%     % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
%     % print(fh, fullfile(subject.figures, [ALSnr '_erp_final_' num2str(tag)]), '-dtiff', '-r200'); close(fh);
%     save_figure(fh, path_figures, [ALSnr '_erp_final_' num2str(tag)], [30 8]);
%
% elseif strcmpi(task, 'RS')
%     % Resting-state
%     dataCmap = brewermap(sum(chaneeg), 'BrBG');
%
%     fh = figure('Visible', opt_visible);
%     th = tiledlayout(1, 2);
%     th.TileSpacing = 'compact'; th.Padding = 'compact';
%
%     % Regular plot
%     th = nexttile;
%     hold on; box off;
%     % plot(freq,psdspectra,'LineWidth',1.1);
%     semilogy(freq, psdspectra, 'LineWidth', 0.8);
%     hold on;
%     semilogy(freq, mean(psdspectra, 2), 'k', 'LineWidth', 2.5); % Add the thick average line
%
%     colororder(th,dataCmap); xlim([freq(1), 60]); % ylim(dataClim);
%     title({ALSnr, [num2str(sum(chaneeg)) ' EEG, ' num2str(NumberTrials(3)) ' trials']});
%     pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('Power');
%
%     % Log-log plot
%     th = nexttile;
%     hold on; box off;
%
%     data_plot   = log10(psdspectra);
%     dataClim    = 1.1 * [min(data_plot(:)), max(data_plot(:))];
%     dataClim(1) = floor(dataClim(1));
%     dataClim(2) = ceil(dataClim(2));
%     plot(freqlog, data_plot, 'LineWidth', 1.1);
%
%     colororder(th, dataCmap); xlim(freqlog([1 end])); ylim(dataClim);
%     title({ALSnr, [num2str(sum(chaneeg)) ' EEG, ' num2str(NumberTrials(3)) ' trials']});
%     pbaspect([1.618 1 1]); xlabel('log_{10}(Frequency) (Hz)'); ylabel('log_{10}(Power)');
%     xticks(freqTicksLog); xticklabels(freqTicksCell);
%
%     % plotX=20; plotY=8;
%     % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
%     % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
%     % print(fh, fullfile(subject.figures, [ALSnr '_pspectra_final_' num2str(tag)]), '-dtiff', '-r200'); close(fh);
%     save_figure(fh, path_figures, [ALSnr '_pspectra_final_' num2str(tag)], [30 8]);
%
% elseif strcmpi(task, 'MT') && has_emg
%     dataCmap1 = brewermap(sum(chaneeg), 'BrBG');
%     dataCmap2 = brewermap(sum(chanemg), 'PRGn');
%
%     fh = figure('Visible', opt_visible);
%     th = tiledlayout(2, 2);
%     th.TileSpacing = 'compact'; th.Padding = 'compact';
%
%     % EEG
%     th = nexttile;
%     hold on; box off;
%     % plot(freq,psdspectra(:,chaneeg),'LineWidth',1.1);
%     semilogy(freq, psdspectra(:, chaneeg), 'LineWidth', 0.8);
%     hold on;
%     semilogy(freq, mean(psdspectra(:, chaneeg), 2), 'k', 'LineWidth', 2.5); % Add the thick average line
%
%     colororder(th,dataCmap1); xlim([freq(1), 60]); % ylim(dataClim);
%     title({ALSnr, [num2str(sum(chaneeg)) ' EEG, ' num2str(NumberTrials(3)) ' trials']});
%     pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('Power');
%
%     % EMG
%     th = nexttile;
%     hold on; box off;
%     % plot(freq,psdspectra(:,chanemg),'LineWidth',1.1);
%     semilogy(freq, psdspectra(:, chanemg), 'LineWidth', 0.8);
%     hold on;
%     semilogy(freq, mean(psdspectra(:, chanemg), 2), 'k', 'LineWidth', 2.5); % Add the thick average line
%
%     colororder(th,dataCmap2); xlim(freq([1 end])); % ylim(dataClim);
%     title({ALSnr, [num2str(sum(chanemg)) ' EMG, ' num2str(NumberTrials(3)) ' trials']});
%     pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('Power');
%
%     % log(EEG)
%     th = nexttile;
%     hold on; box off;
%
%     data_plot   = log10(psdspectra(:, chaneeg));
%     dataClim    = 1.1 * [min(data_plot(:)), max(data_plot(:))];
%     dataClim(1) = floor(dataClim(1));
%     dataClim(2) = ceil(dataClim(2));
%     plot(freqlog,data_plot,'LineWidth',1.1);
%
%     colororder(th,dataCmap1); xlim(freqlog([1 end])); ylim(dataClim);
%     title({ALSnr, [num2str(sum(chaneeg)) ' EEG, ' num2str(NumberTrials(3)) ' trials']});
%     pbaspect([1.618 1 1]); xlabel('log_{10}(Frequency) (Hz)'); ylabel('log_{10}(Power)');
%     xticks(freqTicksLog); xticklabels(freqTicksCell);
%
%     % log(EMG)
%     th = nexttile;
%     hold on; box off;
%
%     data_plot   = log10(psdspectra(:, chanemg));
%     dataClim    = 1.1 * [min(data_plot(:)), max(data_plot(:))];
%     dataClim(1) = floor(dataClim(1));
%     dataClim(2) = ceil(dataClim(2));
%     plot(freqlog,data_plot,'LineWidth',1.1);
%
%     colororder(th,dataCmap2); xlim(freqlog([1 end])); ylim(dataClim);
%     title({ALSnr, [num2str(sum(chanemg)) ' EMG, ' num2str(NumberTrials(3)) ' trials']});
%     pbaspect([1.618 1 1]); xlabel('log_{10}(Frequency) (Hz)'); ylabel('log_{10}(Power)');
%     xticks(freqTicksLog); xticklabels(freqTicksCell);
%
%     % plotX=20; plotY=14;
%     % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
%     % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
%     % print(fh, fullfile(subject.figures, [ALSnr '_pspectra_final_' num2str(tag)]), '-dtiff', '-r200'); close(fh);
%     save_figure(fh, path_figures, [ALSnr '_pspectra_final_' num2str(tag)], [20 14]);
%
% end
%
% % =============================
% % Plot 2
% % =============================
% plot_iaf(DATA, freq, mean(psdspectra(:, chaneeg), 2), opt_visible);
%
% % fprintf('Done!\n');
%
% end