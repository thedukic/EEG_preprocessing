function DATA = check_blink_residuals(DATA, cfg)
% CHECK_BLINK_RESIDUALS
% Evaluates residual vertical blink artifact leakage in cleaned, epoched EEG
% data using event-locked snippet averaging (butterfly plots and topographies),
% and prunes compromised epochs directly from the returned dataset.
%
% Fronto-polar evaluation clusters:
%   - Left Eye  : C29 (Fp1), C30 (AF7/AFp1)
%   - Right Eye : C16 (Fp2), C8  (AF8/AFp2)
%
% Inputs:
%   DATA : Cleaned, epoched EEGLAB dataset (scalp channels + bipolar VEOG)
%   cfg  : (Optional) Configuration structure:
%          cfg.veog_chan       : Label or index of VEOG in DATA (default: 'VEOG')
%          cfg.trIQRblink      : IQR multiplier for VEOG blink detection (default: 3.0)
%          cfg.min_phys_blink  : Minimum absolute VEOG deflection in uV (default: 50 uV)
%          cfg.max_res_uV      : Max tolerable fronto-polar residual peak in uV (default: 15 uV)
%          cfg.max_p2p_uV      : (Optional) Max tolerable peak-to-peak amplitude in uV
%          cfg.win_sec         : Half-window for ERP snippet extraction (default: 0.400 s)
%          cfg.do_plot         : Render 2x2 diagnostic summary figure (default: true)
%
% Outputs:
%   DATA : Updated EEGLAB dataset with compromised blink epochs removed and
%          diagnostics stored in DATA.ALSUTRECHT.leftovers.blink2.blink_stats

% -------------------------------------------------------------------------
% 1. Parameter Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'veog_chan'),      cfg.veog_chan      = 'VEOG'; end
if ~isfield(cfg, 'trIQRblink'),     cfg.trIQRblink     = 3.0;    end
if ~isfield(cfg, 'min_phys_blink'), cfg.min_phys_blink = 50;     end % uV
if ~isfield(cfg, 'max_res_uV'),     cfg.max_res_uV     = 20;     end % uV
if ~isfield(cfg, 'win_sec'),        cfg.win_sec        = 0.400;  end % s (+/- 400 ms)
if ~isfield(cfg, 'do_plot'),        cfg.do_plot        = true;   end

fs = DATA.srate;

% -------------------------------------------------------------------------
% 2. Extract Scalp Channels & Bipolar VEOG
% -------------------------------------------------------------------------
ch_idx = 1:min(128, DATA.nbchan);
[n_chans, n_pnts, n_trials] = size(DATA.data(ch_idx, :, :));

% Locate VEOG channel
if ischar(cfg.veog_chan) || isstring(cfg.veog_chan)
    idx_veog = find(strcmpi({DATA.chanlocs.labels}, cfg.veog_chan), 1);
else
    idx_veog = cfg.veog_chan;
end
assert(~isempty(idx_veog), 'VEOG channel [%s] not found in DATA.', string(cfg.veog_chan));

veog_3d = double(squeeze(DATA.data(idx_veog, :, :))); % [pnts x trials]
veog_2d = veog_3d(:)';

% Locate Fronto-Polar channel clusters
labels = {DATA.chanlocs(ch_idx).labels};

% Left Eye: C29 (Fp1) and C30 (AF7/AFp1)
idx_c29 = find(strcmpi(labels, 'C29') | strcmpi(labels, 'Fp1'), 1);
idx_c30 = find(strcmpi(labels, 'C30') | strcmpi(labels, 'AF7'), 1);
idx_left_cluster = [idx_c29, idx_c30];

% Right Eye: C16 (Fp2) and C8 (AF8/AFp2)
idx_c16 = find(strcmpi(labels, 'C16') | strcmpi(labels, 'Fp2'), 1);
idx_c8  = find(strcmpi(labels, 'C8')  | strcmpi(labels, 'AF8'), 1);
idx_right_cluster = [idx_c16, idx_c8];

ocular_indices = [idx_left_cluster, idx_right_cluster];
assert(~isempty(ocular_indices), 'Could not resolve fronto-polar channel clusters (C29, C30, C16, C8).');

% -------------------------------------------------------------------------
% 3. Dynamic Blink Detection in VEOG
% -------------------------------------------------------------------------
veog_detrend = veog_2d - movmedian(veog_2d, round(1.0 * fs));
veog_IQR     = iqr(veog_detrend);
veog_75P     = prctile(veog_detrend, 75);

thresh_veog   = max(cfg.min_phys_blink, veog_75P + (cfg.trIQRblink * veog_IQR));
blink_mask    = (veog_detrend > thresh_veog);
blink_mask_3d = reshape(blink_mask, n_pnts, n_trials);
blink_epochs  = any(blink_mask_3d, 1);

% -------------------------------------------------------------------------
% 4. Event-Locked Snippet Extraction & Residual Evaluation
% -------------------------------------------------------------------------
win_pts = round(cfg.win_sec * fs);
t_vec   = (-win_pts : win_pts) * (1000 / fs); % Time axis in ms
[~, t_zero_idx] = min(abs(t_vec));

% Baseline window: -350 ms to -250 ms prior to blink peak (quiescent period)
base_idx = find(t_vec >= -350 & t_vec <= -250);
if isempty(base_idx), base_idx = 1:max(2, round(0.050 * fs)); end

% Evaluation window around blink peak: -75 ms to +75 ms
search_idx = find(t_vec >= -100 & t_vec <= 100);
if isempty(search_idx), search_idx = t_zero_idx; end

all_snippets_scalp = [];
all_snippets_veog  = [];
snippet_trials     = [];
bad_blink_trials   = [];

for tr = 1:n_trials
    if blink_epochs(tr)
        v_tr = veog_3d(:, tr);
        v_tr_det = v_tr - median(v_tr);

        % Detect dominant ocular event (positive blink or downward saccade)
        [max_pos, t_pos] = max(v_tr_det);
        [max_neg, t_neg] = min(v_tr_det);

        if abs(max_pos) >= abs(max_neg)
            t_peak = t_pos;
        else
            t_peak = t_neg;
        end

        if (t_peak - win_pts >= 1) && (t_peak + win_pts <= n_pnts)
            idx_w = (t_peak - win_pts) : (t_peak + win_pts);

            % Scalp snippet across all channels
            scalp_snip = double(DATA.data(ch_idx, idx_w, tr));

            % Baseline correct relative to pre-blink period (-350 to -250 ms)
            pre_base = mean(scalp_snip(:, base_idx), 2);
            scalp_snip_corr = scalp_snip - pre_base;

            % --- Evaluate Residuals Across All 4 Fronto-Polar Channels Separately ---
            ocular_traces = scalp_snip_corr(ocular_indices, search_idx); % [4 x n_search_pts]

            % 1. Maximum absolute deflection across all individual channels
            max_res_deflection = max(abs(ocular_traces(:)));

            % 2. Peak-to-peak amplitude per channel
            p2p_per_ch = max(ocular_traces, [], 2) - min(ocular_traces, [], 2);
            max_p2p    = max(p2p_per_ch);

            is_over_or_undercleaned = max_res_deflection > cfg.max_res_uV;
            is_biphasic_artifact   = isfield(cfg, 'max_p2p_uV') && (max_p2p > cfg.max_p2p_uV);

            if is_over_or_undercleaned || is_biphasic_artifact
                bad_blink_trials(end+1) = tr; %#ok<AGROW>
            end

            all_snippets_scalp(:, :, end+1) = scalp_snip_corr; %#ok<AGROW>
            all_snippets_veog(:, :, end+1)  = v_tr(idx_w) - mean(v_tr(idx_w(base_idx))); %#ok<AGROW>
            snippet_trials(end+1)           = tr; %#ok<AGROW>
        end
    end
end

bad_blink_trials = unique(bad_blink_trials);
clean_trials     = setdiff(1:n_trials, bad_blink_trials);
n_clean_trials   = length(clean_trials);

% -------------------------------------------------------------------------
% 5. Compute Pre- and Post-Pruning Blink-Locked Averages
% -------------------------------------------------------------------------
% Pre-Pruning (All detected blink events)
if ~isempty(all_snippets_scalp)
    erp_pre_scalp = mean(all_snippets_scalp, 3);
    erp_pre_veog  = mean(all_snippets_veog, 3);
    topo_pre      = erp_pre_scalp(:, t_zero_idx);
    peak_pre_fp   = max(abs(erp_pre_scalp(ocular_indices, :)), [], 'all');
else
    erp_pre_scalp = zeros(n_chans, length(t_vec));
    erp_pre_veog  = zeros(1, length(t_vec));
    topo_pre      = zeros(n_chans, 1);
    peak_pre_fp   = 0;
end

% Post-Pruning (Only snippets from clean, retained epochs)
clean_snip_idx = find(~ismember(snippet_trials, bad_blink_trials));

if ~isempty(clean_snip_idx)
    erp_post_scalp = mean(all_snippets_scalp(:, :, clean_snip_idx), 3);
    erp_post_veog  = mean(all_snippets_veog(:, :, clean_snip_idx), 3);
    topo_post      = erp_post_scalp(:, t_zero_idx);
    peak_post_fp   = max(abs(erp_post_scalp(ocular_indices, :)), [], 'all');
else
    erp_post_scalp = zeros(n_chans, length(t_vec));
    erp_post_veog  = zeros(1, length(t_vec));
    topo_post      = zeros(n_chans, 1);
    peak_post_fp   = 0;
end

% -------------------------------------------------------------------------
% 6. Pack Statistics Structure
% -------------------------------------------------------------------------
blink_stats.t_vec             = t_vec;
blink_stats.erp_pre_scalp     = erp_pre_scalp;
blink_stats.erp_pre_veog      = erp_pre_veog;
blink_stats.erp_post_scalp    = erp_post_scalp;
blink_stats.erp_post_veog     = erp_post_veog;
blink_stats.topo_pre          = topo_pre;
blink_stats.topo_post         = topo_post;
blink_stats.peak_pre_fp       = peak_pre_fp;
blink_stats.peak_post_fp      = peak_post_fp;
blink_stats.thresh_veog       = thresh_veog;
blink_stats.n_blink_epochs    = sum(blink_epochs);
blink_stats.n_blink_snippets  = size(all_snippets_scalp, 3);
blink_stats.bad_blink_trials  = bad_blink_trials;
blink_stats.clean_trials      = clean_trials;
blink_stats.n_clean_trials    = n_clean_trials;

% -------------------------------------------------------------------------
% 7. Console Diagnostic Summary
% -------------------------------------------------------------------------
fprintf('  Total Epochs Analysed:        %d\n', n_trials);
fprintf('  Blink Epochs Detected:        %d (Threshold: %.1f uV)\n', sum(blink_epochs), thresh_veog);
fprintf('  Blink Events Extracted:       %d\n', size(all_snippets_scalp, 3));
fprintf('  -----------------------------------------------------------------\n');
fprintf('  PRE-PRUNING (All Blink Events, N = %d):\n', size(all_snippets_scalp, 3));
fprintf('    - Peak VEOG Amplitude:      %.1f uV\n', max(abs(erp_pre_veog)));
fprintf('    - Max Residual Fronto-Polar: %.2f uV\n', peak_pre_fp);
fprintf('    - Compromised Epochs (>%duV): %d / %d (%.1f%%)\n', ...
    cfg.max_res_uV, length(bad_blink_trials), sum(blink_epochs), (length(bad_blink_trials)/max(1, sum(blink_epochs)))*100);
fprintf('  -----------------------------------------------------------------\n');
fprintf('  POST-PRUNING (Clean Epochs, N = %d / %d | %.1f%% Retained):\n', ...
    n_clean_trials, n_trials, (n_clean_trials / n_trials) * 100);
fprintf('    - Clean Fronto-Polar ERP:   %.2f uV (Target: < 2.0 uV)\n', peak_post_fp);

% -------------------------------------------------------------------------
% 8. Diagnostic Figure Visualisation (2 Rows x 2 Columns)
% -------------------------------------------------------------------------
if cfg.do_plot
    fh = figure('Color', 'w', 'Position', [100, 100, 1150, 750], 'Name', 'Blink Residual Diagnostics');

    % Y-limits for ERP plots
    max_erp_val = max([max(abs(erp_pre_scalp(:))), max(abs(erp_pre_veog(:)/5)), 10]);
    y_lims = [-ceil(max_erp_val * 1.1), ceil(max_erp_val * 1.1)];

    % Symmetric colour limits for topoplots (floor at +/- 5 uV)
    min_topo_floor = 5; % uV
    max_observed   = max([max(abs(topo_pre(:))), max(abs(topo_post(:)))]);
    topo_lim_val   = max(min_topo_floor, ceil(max_observed));
    topo_lims      = [-topo_lim_val, topo_lim_val];

    % --- Row 1, Col 1: Pre-Pruning Butterfly Plot ---
    subplot(2, 2, 1);
    plot(t_vec, erp_pre_scalp', 'Color', [0.80 0.80 0.80 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    hold on;
    % Cluster averages for visualization
    erp_pre_left  = mean(erp_pre_scalp(idx_left_cluster, :), 1);
    erp_pre_right = mean(erp_pre_scalp(idx_right_cluster, :), 1);
    h_left  = plot(t_vec, erp_pre_left,  'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'DisplayName', 'Left Eye (C29, C30)');
    h_right = plot(t_vec, erp_pre_right, 'Color', [0.15 0.65 0.20], 'LineWidth', 1.8, 'DisplayName', 'Right Eye (C16, C8)');
    h_v     = plot(t_vec, erp_pre_veog / 5, '--k', 'LineWidth', 1.2, 'DisplayName', 'VEOG / 5');
    grid on; xlim([min(t_vec), max(t_vec)]); ylim(y_lims);
    xlabel('Time from Blink Peak (ms)', 'FontSize', 9);
    ylabel('Residual Voltage (\muV)', 'FontSize', 9);
    title(sprintf('Pre-Pruning Blink ERP (N = %d, Peak = %.1f \\muV)', ...
        size(all_snippets_scalp, 3), peak_pre_fp), 'FontWeight', 'bold', 'FontSize', 10);
    legend([h_left, h_right, h_v], 'Location', 'northeast', 'FontSize', 8);

    % --- Row 1, Col 2: Pre-Pruning Topography at t = 0 ms ---
    sbh2 = subplot(2, 2, 2);
    title_t1 = 'Pre-Pruning Residual (t = 0 ms)';
    render_topoplot(topo_pre, DATA.chanlocs(ch_idx), title_t1, sbh2, topo_lims);

    % --- Row 2, Col 1: Post-Pruning Butterfly Plot ---
    subplot(2, 2, 3);
    plot(t_vec, erp_post_scalp', 'Color', [0.80 0.80 0.80 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    hold on;
    erp_post_left  = mean(erp_post_scalp(idx_left_cluster, :), 1);
    erp_post_right = mean(erp_post_scalp(idx_right_cluster, :), 1);
    h_left_post  = plot(t_vec, erp_post_left,  'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'DisplayName', 'Left Eye (C29, C30)');
    h_right_post = plot(t_vec, erp_post_right, 'Color', [0.15 0.65 0.20], 'LineWidth', 1.8, 'DisplayName', 'Right Eye (C16, C8)');
    grid on; xlim([min(t_vec), max(t_vec)]); ylim(y_lims);
    xlabel('Time from Blink Peak (ms)', 'FontSize', 9);
    ylabel('Residual Voltage (\muV)', 'FontSize', 9);
    title(sprintf('Post-Pruning Clean ERP (N = %d Clean Epochs, Peak = %.1f \\muV)', ...
        length(clean_snip_idx), peak_post_fp), 'FontWeight', 'bold', 'FontSize', 10);
    legend([h_left_post, h_right_post], 'Location', 'northeast', 'FontSize', 8);

    % --- Row 2, Col 2: Post-Pruning Topography at t = 0 ms ---
    sbh4 = subplot(2, 2, 4);
    title_t2 = 'Post-Pruning Clean Residual (t = 0 ms)';
    render_topoplot(topo_post, DATA.chanlocs(ch_idx), title_t2, sbh4, topo_lims);

    drawnow;

    % Save figure if pipeline paths exist
    if isfield(DATA, 'ALSUTRECHT') && isfield(DATA.ALSUTRECHT, 'subject') && isfield(DATA.ALSUTRECHT.subject, 'figures')
        save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_leftovers_2'], [20 20]);
    end
end

% -------------------------------------------------------------------------
% 9. Log Statistics to Dataset Struct
% -------------------------------------------------------------------------
DATA.ALSUTRECHT.leftovers.blink2.blink_stats = blink_stats;

% -------------------------------------------------------------------------
% 10. Prune Compromised Blink Epochs from Dataset
% -------------------------------------------------------------------------
if ~isempty(bad_blink_trials)
    fprintf('Pruning %d compromised blink epoch(s) from dataset...\n', length(bad_blink_trials));

    if exist('pop_select', 'file') == 2
        DATA = pop_select(DATA, 'notrial', bad_blink_trials);
    else
        DATA.data   = DATA.data(:, :, clean_trials);
        DATA.trials = length(clean_trials);
        if isfield(DATA, 'epoch') && ~isempty(DATA.epoch)
            DATA.epoch = DATA.epoch(clean_trials);
        end
    end
    fprintf('Dataset successfully updated: %d clean epochs retained.\n', DATA.trials);
else
    fprintf('No epochs required pruning; all %d epochs retained.\n', DATA.trials);
end

end

% =========================================================================
% Local Helper: Topoplot Dispatcher with Colourbar
% =========================================================================
function render_topoplot(values, chanlocs, plot_title, ax_handle, map_limits)
axes(ax_handle);
if exist('mytopoplot', 'file') == 2
    mytopoplot(values, false(size(values)), plot_title, ax_handle, map_limits);
    cb = colorbar('peer', ax_handle);
    ylabel(cb, '\muV', 'FontSize', 9, 'FontWeight', 'bold');
elseif exist('topoplot', 'file') == 2
    topoplot(values, chanlocs, 'maplimits', map_limits, 'electrodes', 'off', 'style', 'map');
    title(plot_title, 'FontSize', 10, 'FontWeight', 'bold');
    cb = colorbar;
    ylabel(cb, '\muV', 'FontSize', 9, 'FontWeight', 'bold');
else
    bar(values);
    title(plot_title, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Residual Voltage (\muV)');
    xlabel('Channels');
    grid on;
end
end