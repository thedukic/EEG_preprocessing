function generate_finalplots(DATA, NumberTrials, thisTask, thisTag, cfg)
% =========================================================================
% GENERATE_FINALPLOTS: BioSemi 128 Quality Control Visualisations
% =========================================================================
% Directly uses central cfg.roi masks for all paradigm-specific overlays.
% =========================================================================

% -------------------------------------------------------------------------
% 1. Config Parsing & Defaults
% -------------------------------------------------------------------------
if nargin < 5; cfg = struct(); end

% Figure visibility
if isfield(cfg, 'figure') && isfield(cfg.figure, 'visible')
    optVisible = cfg.figure.visible;
elseif isfield(cfg, 'visible')
    optVisible = cfg.visible;
else
    optVisible = 'on';
end

% Ensure central ROI definitions exist in cfg (BioSemi 128 defaults)
if ~isfield(cfg, 'roi'); cfg.roi = struct(); end

if ~isfield(cfg.roi, 'mmn')
    cfg.roi.mmn = {'C21', 'C20', 'C19', 'C22', 'C23', 'C25', 'C26', 'C12', 'C13', 'D1', 'D2', 'C1', 'C2'};
end
if ~isfield(cfg.roi, 'sart_p300')
    cfg.roi.sart_p300 = {'A19', 'A4', 'A3', 'A2', 'A20', 'A21', 'A5', 'A18', 'A32', 'A31'};
end
if ~isfield(cfg.roi, 'rs_alpha')
    cfg.roi.rs_alpha = {'A23', 'A22', 'A21', 'A24', 'A25', 'A15', 'A16', 'A17', 'A28', 'A29', 'A30', 'A14', 'A27'};
end
if ~isfield(cfg.roi, 'motor_left')
    cfg.roi.motor_left = {'D19', 'D18', 'D14', 'D20', 'D21', 'D12', 'D11', 'D13', 'D17', 'D28', 'D27'};
end
if ~isfield(cfg.roi, 'motor_right')
    cfg.roi.motor_right = {'B22', 'B21', 'B20', 'B23', 'B24', 'B31', 'B30', 'B32', 'B19', 'B18', 'B17'};
end

% -------------------------------------------------------------------------
% 2. Extract Data & Metadata
% -------------------------------------------------------------------------
fprintf('\n================================\n');
fprintf('Final QC Reports: %s (BioSemi 128 Montage)\n', upper(thisTask));
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

% Visual styling
col_gray_trace = [0.25 0.25 0.25 0.15];
col_ga         = [0.05 0.05 0.05];
col_emg        = [0.85 0.33 0.10];

% =========================================================================
% 3. Cognitive ERPs (MMN / SART)
% =========================================================================
if any(strcmpi(thisTask, {'MMN', 'SART'}))

    triggers_list = [DATA.event.edftype];

    if strcmpi(thisTask, 'SART')
        roi_target = cfg.roi.sart_p300;
        roi_name   = 'Parietal Pz ROI (A-Bank)';
        if strcmpi(DATA.ALSUTRECHT.SART.type, 'StimulusLocked')
            triggers_unique = [3 6];
            cond_names = {'Frequent', 'Target (No-Go)'};
        else
            triggers_unique = 1;
            cond_names = {'Response-Locked'};
        end
        analysis_window = [300 500];
        topo_title = 'P300 Topo (300-500 ms)';

    elseif strcmpi(thisTask, 'MMN')
        roi_target = cfg.roi.mmn;
        roi_name   = 'Frontal Fz/FCz ROI (C-Bank)';
        triggers_unique = [17 12];
        cond_names = {'Standard', 'Deviant'};
        analysis_window = [100 220];
        topo_title = 'MMN Diff (<20 Hz, 100-220 ms)';
    end

    num_conds = length(triggers_unique);
    roi_mask  = ismember(eeg_labels, roi_target);
    if ~any(roi_mask); roi_mask = true(size(eeg_labels)); end

    erp_data = cell(num_conds, 1);
    n_trials_cond = zeros(num_conds, 1);
    for i_c = 1:num_conds
        c_mask = (triggers_list == triggers_unique(i_c));
        erp_data{i_c} = mean(DATA.data(chaneeg, :, c_mask), 3);
        n_trials_cond(i_c) = sum(c_mask);
    end

    total_tiles = num_conds + (num_conds == 2) + 1;
    fh = figure('Visible', optVisible, 'Color', [1 1 1]);
    tiled_h = tiledlayout(1, total_tiles, 'TileSpacing', 'compact', 'Padding', 'compact');

    all_erp_vals = cell2mat(erp_data);
    y_lims = 1.15 * [min(all_erp_vals(:)), max(all_erp_vals(:))];
    y_lims = [min(y_lims(1), -4), max(y_lims(2), 4)];

    for i_c = 1:num_conds
        ax = nexttile(tiled_h);
        hold(ax, 'on'); box(ax, 'off');

        patch(ax, [analysis_window(1) analysis_window(2) analysis_window(2) analysis_window(1)], ...
            [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], ...
            [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');

        xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        yline(ax, 0, '-',  'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');

        plot(ax, DATA.times, erp_data{i_c}, 'Color', col_gray_trace, 'LineWidth', 0.6, 'HandleVisibility', 'off');
        p_ga  = plot(ax, DATA.times, mean(erp_data{i_c}, 1), 'Color', col_ga, 'LineWidth', 1.8, 'DisplayName', 'Grand Avg (128)');
        p_roi = plot(ax, DATA.times, mean(erp_data{i_c}(roi_mask, :), 1), 'Color', [0.00 0.45 0.74], 'LineWidth', 2.2, 'DisplayName', roi_name);

        xlim(ax, [DATA.times(1), DATA.times(end)]);
        ylim(ax, y_lims);
        xlabel(ax, 'Time (ms)'); ylabel(ax, 'Amplitude (\muV)');
        title(ax, sprintf('%s: %s (N=%d)', ALSnr, cond_names{i_c}, n_trials_cond(i_c)), 'Interpreter', 'none');
        legend(ax, [p_ga, p_roi], 'Location', 'northwest', 'FontSize', 8);
        pbaspect(ax, [1.4 1 1]);
    end

    if num_conds == 2
        diff_wave = erp_data{2} - erp_data{1};
        ax = nexttile(tiled_h);
        hold(ax, 'on'); box(ax, 'off');

        patch(ax, [analysis_window(1) analysis_window(2) analysis_window(2) analysis_window(1)], ...
            [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], ...
            [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');

        xline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
        yline(ax, 0, '-',  'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');

        plot(ax, DATA.times, diff_wave, 'Color', col_gray_trace, 'LineWidth', 0.6, 'HandleVisibility', 'off');
        plot(ax, DATA.times, mean(diff_wave, 1), 'Color', col_ga, 'LineWidth', 1.8, 'DisplayName', 'Grand Avg (128)');
        plot(ax, DATA.times, mean(diff_wave(roi_mask, :), 1), 'Color', [0.85 0.10 0.10], 'LineWidth', 2.2, 'DisplayName', [roi_name ' Diff']);

        xlim(ax, [DATA.times(1), DATA.times(end)]);
        ylim(ax, y_lims);
        xlabel(ax, 'Time (ms)'); ylabel(ax, 'Amplitude (\muV)');
        title(ax, sprintf('%s: Diff (%s - %s)', ALSnr, cond_names{2}, cond_names{1}), 'Interpreter', 'none');
        legend(ax, 'Location', 'northwest', 'FontSize', 8);
        pbaspect(ax, [1.4 1 1]);
    end

    ax_topo = nexttile(tiled_h);
    if num_conds == 2; target_topo = diff_wave; else; target_topo = erp_data{1}; end
    time_mask = (DATA.times >= analysis_window(1)) & (DATA.times <= analysis_window(2));
    topo_vals = mean(target_topo(:, time_mask), 2);

    if strcmpi(thisTask, 'MMN')
        [bl, al] = butter(2, 20 / (DATA.srate / 2), 'low');
        target_topo_filt = filtfilt(bl, al, double(target_topo)')';
        topo_vals = mean(target_topo_filt(:, time_mask), 2);
    end

    max_c = max(abs(topo_vals));
    if max_c == 0; max_c = 1; end
    mytopoplot(topo_vals, DATA.chanlocs(chaneeg), topo_title, ax_topo, [-max_c max_c]);

    save_figure(fh, path_figures, sprintf('%s_erp_final_%s', ALSnr, num2str(thisTag)), [32 9]);
    if strcmpi(optVisible, 'off'); close(fh); end

    % =========================================================================
    % 4. Resting-State (RS) PSD
    % =========================================================================
elseif strcmpi(thisTask, 'RS')

    psd_eeg = psdspectra(:, chaneeg);
    psd_db  = 10 * log10(psd_eeg + eps);
    alpha_roi_mask = ismember(eeg_labels, cfg.roi.rs_alpha);

    fh = figure('Visible', optVisible, 'Color', [1 1 1]);
    tiled_h = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile(tiled_h);
    hold(ax1, 'on'); box(ax1, 'off');

    y_range = [min(psd_db(:)) - 2, max(psd_db(:)) + 2];
    patch(ax1, [8 12 12 8], [y_range(1) y_range(1) y_range(2) y_range(2)], ...
        [0.92 0.96 0.92], 'EdgeColor', 'none', 'DisplayName', 'Alpha (8-12 Hz)');
    patch(ax1, [15 30 30 15], [y_range(1) y_range(1) y_range(2) y_range(2)], ...
        [0.93 0.93 0.98], 'EdgeColor', 'none', 'DisplayName', 'Beta (15-30 Hz)');

    plot(ax1, freq, psd_db, 'Color', col_gray_trace, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    plot(ax1, freq, mean(psd_db, 2), 'Color', col_ga, 'LineWidth', 2.2, 'DisplayName', 'Grand Avg (128)');

    if any(alpha_roi_mask)
        plot(ax1, freq, mean(psd_db(:, alpha_roi_mask), 2), ...
            'Color', [0.10 0.65 0.10], 'LineWidth', 2.2, 'DisplayName', 'Occipital Alpha (A-Bank)');
    end

    xlim(ax1, [1 70]); ylim(ax1, y_range);
    xlabel(ax1, 'Frequency (Hz)'); ylabel(ax1, 'Power Spectral Density (dB/Hz)');
    title(ax1, sprintf('%s: Resting-State PSD', ALSnr), 'Interpreter', 'none');
    legend(ax1, 'Location', 'southeast', 'FontSize', 8);
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
    legend(ax2, 'Location', 'southwest', 'FontSize', 8);
    pbaspect(ax2, [1.5 1 1]);

    save_figure(fh, path_figures, sprintf('%s_psd_final_%s', ALSnr, num2str(thisTag)), [24 9]);
    if strcmpi(optVisible, 'off'); close(fh); end

    % =========================================================================
    % 5. Motor Task (MT) PSD (Consuming cfg.roi.motor_left / right)
    % =========================================================================
elseif strcmpi(thisTask, 'MT')

    psd_eeg    = psdspectra(:, chaneeg);
    psd_eeg_db = 10 * log10(psd_eeg + eps);

    % Use central cfg.roi masks
    mask_left  = ismember(eeg_labels, cfg.roi.motor_left);
    mask_right = ismember(eeg_labels, cfg.roi.motor_right);

    fh = figure('Visible', optVisible, 'Color', [1 1 1]);
    tiled_h = tiledlayout(1, 2 + has_emg, 'TileSpacing', 'compact', 'Padding', 'compact');

    % A. Motor EEG Spectrum
    ax1 = nexttile(tiled_h);
    hold(ax1, 'on'); box(ax1, 'off');

    y_lim_eeg = [min(psd_eeg_db(:)) - 2, max(psd_eeg_db(:)) + 2];
    patch(ax1, [15 30 30 15], [y_lim_eeg(1) y_lim_eeg(1) y_lim_eeg(2) y_lim_eeg(2)], ...
        [0.93 0.93 0.98], 'EdgeColor', 'none', 'DisplayName', 'Beta CMC (15-30 Hz)');
    patch(ax1, [30 45 45 30], [y_lim_eeg(1) y_lim_eeg(1) y_lim_eeg(2) y_lim_eeg(2)], ...
        [0.98 0.93 0.93], 'EdgeColor', 'none', 'DisplayName', 'Piper (30-45 Hz)');

    plot(ax1, freq, psd_eeg_db, 'Color', col_gray_trace, 'LineWidth', 0.6, 'HandleVisibility', 'off');
    plot(ax1, freq, mean(psd_eeg_db, 2), 'Color', col_ga, 'LineWidth', 1.8, 'DisplayName', 'Grand Avg (128)');

    if any(mask_left)
        plot(ax1, freq, mean(psd_eeg_db(:, mask_left), 2), ...
            'Color', [0.00 0.45 0.74], 'LineWidth', 2.2, 'DisplayName', 'Left Motor (D-Bank / C3)');
    end

    if any(mask_right)
        plot(ax1, freq, mean(psd_eeg_db(:, mask_right), 2), ...
            'Color', [0.85 0.15 0.15], 'LineWidth', 2.2, 'DisplayName', 'Right Motor (B-Bank / C4)');
    end

    xlim(ax1, [1 48]); ylim(ax1, y_lim_eeg);
    xlabel(ax1, 'Frequency (Hz)'); ylabel(ax1, 'EEG Power (dB/Hz)');
    title(ax1, sprintf('%s: Sensorimotor PSD (Left vs Right)', ALSnr), 'Interpreter', 'none');
    legend(ax1, 'Location', 'southwest', 'FontSize', 8);
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
        legend(ax2, 'Location', 'southwest', 'FontSize', 8);
        pbaspect(ax2, [1.4 1 1]);
    end

    % C. Topography of Beta Band (15 to 30 Hz) Power
    ax_topo = nexttile(tiled_h);
    beta_mask = (freq >= 15) & (freq <= 30);
    beta_power = mean(psd_eeg(beta_mask, :), 1)';

    mytopoplot(beta_power, DATA.chanlocs(chaneeg), 'Beta Power (15-30 Hz)', ax_topo, []);

    save_figure(fh, path_figures, sprintf('%s_psd_final_%s', ALSnr, num2str(thisTag)), [28 9]);
    if strcmpi(optVisible, 'off'); close(fh); end
end

% Diagnostic IAF Peak
plot_iaf(DATA, freq, mean(psdspectra(:, chaneeg), 2), optVisible);

end

function plot_iaf(EEG, freq1, psdspectra1, optVisible)
% =========================================================================
% PLOT_IAF: Individual Alpha Frequency (IAF / PAF) Diagnostic Report
% =========================================================================
% Overlays the whole-cap average power spectrum with the dedicated IAF channel
% spectrum, canonical frequency bands, and the peak alpha frequency marker.
% =========================================================================

if nargin < 4 || isempty(optVisible); optVisible = 'on'; end

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
fh = figure('Visible', optVisible, 'Color', [1 1 1]);
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
legend_handles = [h_avg];
if ~isempty(h_iaf); legend_handles = [legend_handles; h_iaf]; end
if ~isempty(h_paf); legend_handles = [legend_handles; h_paf]; end
legend_handles = [legend_handles; h_bands];

legend(ax, legend_handles, 'Location', 'eastoutside', 'FontSize', 8);

% -------------------------------------------------------------------------
% 5. Save Figure
% -------------------------------------------------------------------------
if isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'subject')
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, ...
        [EEG.ALSUTRECHT.subject.id '_pspectra_iaf'], [20 8]);
end

if strcmpi(optVisible, 'off'); close(fh); end

end

% function generate_finalplots(DATA, NumberTrials, thisTask, thisTag, optVisible)
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
% if strcmpi(thisTask, 'MMN') || strcmpi(thisTask, 'SART')
%     % ERP
%     triggers_list = [DATA.event.edftype];
%     if strcmpi(thisTask, 'SART')
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
%     elseif strcmpi(thisTask, 'MMN')
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
%     fh = figure('Visible', optVisible);
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
%     if strcmpi(thisTask, 'MMN')
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
%     % print(fh, fullfile(subject.figures, [ALSnr '_erp_final_' num2str(thisTag)]), '-dtiff', '-r200'); close(fh);
%     save_figure(fh, path_figures, [ALSnr '_erp_final_' num2str(thisTag)], [30 8]);
%
% elseif strcmpi(thisTask, 'RS')
%     % Resting-state
%     dataCmap = brewermap(sum(chaneeg), 'BrBG');
%
%     fh = figure('Visible', optVisible);
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
%     % print(fh, fullfile(subject.figures, [ALSnr '_pspectra_final_' num2str(thisTag)]), '-dtiff', '-r200'); close(fh);
%     save_figure(fh, path_figures, [ALSnr '_pspectra_final_' num2str(thisTag)], [30 8]);
%
% elseif strcmpi(thisTask, 'MT') && has_emg
%     dataCmap1 = brewermap(sum(chaneeg), 'BrBG');
%     dataCmap2 = brewermap(sum(chanemg), 'PRGn');
%
%     fh = figure('Visible', optVisible);
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
%     % print(fh, fullfile(subject.figures, [ALSnr '_pspectra_final_' num2str(thisTag)]), '-dtiff', '-r200'); close(fh);
%     save_figure(fh, path_figures, [ALSnr '_pspectra_final_' num2str(thisTag)], [20 14]);
%
% end
%
% % =============================
% % Plot 2
% % =============================
% plot_iaf(DATA, freq, mean(psdspectra(:, chaneeg), 2), optVisible);
%
% % fprintf('Done!\n');
%
% end