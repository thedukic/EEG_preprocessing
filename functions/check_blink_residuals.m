function DATA = check_blink_residuals(DATA, tag, cfg)
% CHECK_BLINK_RESIDUALS
% Evaluates residual vertical blink artefact leakage using event-locked
% snippet averaging (butterfly plots and topographies), and prunes
% compromised epochs directly from the returned dataset.

if nargin < 2 || isempty(tag), tag = 'xxx'; end
if nargin < 3, cfg = struct(); end

if ~isfield(cfg, 'veog_chan'),      cfg.veog_chan      = 'VEOG'; end
if ~isfield(cfg, 'trIQRblink'),     cfg.trIQRblink     = 3.0;    end
if ~isfield(cfg, 'min_phys_blink'), cfg.min_phys_blink = 50;     end % uV
if ~isfield(cfg, 'max_res_uV'),     cfg.max_res_uV     = 25;     end % uV
if ~isfield(cfg, 'max_p2p_uV'),     cfg.max_p2p_uV     = 40;     end % uV
if ~isfield(cfg, 'win_sec'),        cfg.win_sec        = 0.400;  end % s
if ~isfield(cfg, 'do_plot'),        cfg.do_plot        = true;   end

% -------------------------------------------------------------------------
% Call Core Detection Engine
% -------------------------------------------------------------------------
cfg_tmp = cfg;
cfg_tmp.do_plot = false;
det = detect_blink_residuals(DATA, cfg_tmp);

ch_idx             = det.ch_idx;
n_chans            = length(ch_idx);
ocular_indices     = det.ocular_indices;
idx_left_cluster   = det.idx_left_cluster;
idx_right_cluster  = det.idx_right_cluster;
all_snippets_scalp = det.all_snippets_scalp;
all_snippets_veog  = det.all_snippets_veog;
snippet_trials     = det.snippet_trials;
bad_blink_trials   = det.bad_blink_trials;
clean_trials       = det.clean_trials;
blink_epochs       = det.blink_epochs;
thresh_veog        = det.thresh_veog;
t_vec              = det.t_vec;
t_zero_idx         = det.t_zero_idx;
is_epoched         = det.is_epoched;
n_clean_trials     = length(clean_trials);

if is_epoched
    n_trials = size(DATA.data, 3);
else
    n_trials = length(blink_epochs);
end

% -------------------------------------------------------------------------
% 5. Compute Pre- and Post-Pruning Blink-Locked Averages & Topographic Power
% -------------------------------------------------------------------------
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

% Average Channel Power: Removed vs Retained Trials
if is_epoched && ~isempty(bad_blink_trials)
    raw_rem = double(DATA.data(ch_idx, :, bad_blink_trials));
    raw_rem = raw_rem - mean(raw_rem, 2);
    topo_power_removed = mean(mean(raw_rem.^2, 2), 3);
else
    topo_power_removed = zeros(n_chans, 1);
end

if is_epoched && ~isempty(clean_trials)
    raw_clean = double(DATA.data(ch_idx, :, clean_trials));
    raw_clean = raw_clean - mean(raw_clean, 2);
    topo_power_clean = mean(mean(raw_clean.^2, 2), 3);
else
    topo_power_clean = zeros(n_chans, 1);
end

% -------------------------------------------------------------------------
% 6. Pack Statistics Structure
% -------------------------------------------------------------------------
blink_stats.t_vec              = t_vec;
blink_stats.erp_pre_scalp      = erp_pre_scalp;
blink_stats.erp_pre_veog       = erp_pre_veog;
blink_stats.erp_post_scalp     = erp_post_scalp;
blink_stats.erp_post_veog      = erp_post_veog;
blink_stats.topo_pre           = topo_pre;
blink_stats.topo_post          = topo_post;
blink_stats.topo_power_removed = topo_power_removed;
blink_stats.topo_power_clean   = topo_power_clean;
blink_stats.peak_pre_fp        = peak_pre_fp;
blink_stats.peak_post_fp       = peak_post_fp;
blink_stats.thresh_veog        = thresh_veog;
blink_stats.n_blink_epochs     = length(blink_epochs);
blink_stats.n_blink_snippets   = size(all_snippets_scalp, 3);
blink_stats.bad_blink_trials   = bad_blink_trials;
blink_stats.clean_trials       = clean_trials;
blink_stats.n_clean_trials     = n_clean_trials;

% -------------------------------------------------------------------------
% 7. Console Diagnostic Summary
% -------------------------------------------------------------------------
fprintf('Total Events/Epochs Analysed: %d\n', n_trials);
fprintf('Blink Events Detected:        %d (Threshold: %.1f uV)\n', length(blink_epochs), thresh_veog);
fprintf('Blink Snippets Extracted:     %d\n', size(all_snippets_scalp, 3));
fprintf('-----------------------------------------------------------------\n');
fprintf('PRE-PRUNING (All Blink Events, N = %d):\n', size(all_snippets_scalp, 3));
fprintf('  - Peak VEOG Amplitude:      %.1f uV\n', max(abs(erp_pre_veog)));
fprintf('  - Max Residual Fronto-Polar: %.2f uV\n', peak_pre_fp);
fprintf('  - Compromised (>%duV):       %d / %d (%.1f%%)\n', ...
    cfg.max_res_uV, length(bad_blink_trials), max(1, length(blink_epochs)), ...
    (length(bad_blink_trials)/max(1, length(blink_epochs)))*100);
fprintf('-----------------------------------------------------------------\n');
fprintf('POST-PRUNING (Clean, N = %d / %d | %.1f%% Retained):\n', ...
    n_clean_trials, n_trials, (n_clean_trials / max(1, n_trials)) * 100);
fprintf('  - Clean Fronto-Polar ERP:   %.2f uV (Target: < 2.0 uV)\n', peak_post_fp);

% -------------------------------------------------------------------------
% 8. Diagnostic Figure Visualisation (2 Rows x 3 Columns)
% -------------------------------------------------------------------------
if cfg.do_plot
    fh = figure('Color', 'w', 'Position', [80, 80, 1400, 720], 'Name', 'Blink Residual Diagnostics');

    max_erp_val = max([max(abs(erp_pre_scalp(:))), max(abs(erp_pre_veog(:)/5)), 10]);
    y_lims = [-ceil(max_erp_val * 1.1), ceil(max_erp_val * 1.1)];

    min_topo_floor = 5;
    max_observed   = max([max(abs(topo_pre(:))), max(abs(topo_post(:)))]);
    topo_lim_val   = max(min_topo_floor, ceil(max_observed));
    topo_lims      = [-topo_lim_val, topo_lim_val];

    power_lim_val  = max([max(topo_power_removed), max(topo_power_clean), 10]);
    power_lims     = [0, ceil(power_lim_val * 1.05)];

    % Row 1, Col 1: Pre-Pruning Butterfly Plot
    subplot(2, 3, 1);
    plot(t_vec, erp_pre_scalp', 'Color', [0.80 0.80 0.80 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    hold on;
    erp_pre_left  = mean(erp_pre_scalp(idx_left_cluster, :), 1);
    erp_pre_right = mean(erp_pre_scalp(idx_right_cluster, :), 1);
    h_left  = plot(t_vec, erp_pre_left,  'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'DisplayName', 'Left Eye');
    h_right = plot(t_vec, erp_pre_right, 'Color', [0.15 0.65 0.20], 'LineWidth', 1.8, 'DisplayName', 'Right Eye');
    h_v     = plot(t_vec, erp_pre_veog / 5, '--k', 'LineWidth', 1.2, 'DisplayName', 'VEOG / 5');
    grid on; xlim([min(t_vec), max(t_vec)]); ylim(y_lims);
    xlabel('Time from Blink Peak (ms)', 'FontSize', 9);
    ylabel('Residual Voltage (\muV)', 'FontSize', 9);
    title(sprintf('Pre-Pruning Blink ERP (N = %d, Peak = %.1f \\muV)', ...
        size(all_snippets_scalp, 3), peak_pre_fp), 'FontWeight', 'bold', 'FontSize', 10);
    legend([h_left, h_right, h_v], 'Location', 'south', 'FontSize', 8);

    % Row 1, Col 2: Pre-Pruning Topography at t = 0 ms
    sbh2 = subplot(2, 3, 2);
    render_topoplot(topo_pre, 'Pre-Pruning Residual (t = 0 ms)', sbh2, topo_lims, '\muV');

    % Row 1, Col 3: Removed Epochs Average Power
    sbh3 = subplot(2, 3, 3);
    render_topoplot(topo_power_removed, sprintf('Removed Power (N = %d)', length(bad_blink_trials)), ...
        sbh3, power_lims, '\muV^2');

    % Row 2, Col 1: Post-Pruning Butterfly Plot
    subplot(2, 3, 4);
    plot(t_vec, erp_post_scalp', 'Color', [0.80 0.80 0.80 0.6], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    hold on;
    erp_post_left  = mean(erp_post_scalp(idx_left_cluster, :), 1);
    erp_post_right = mean(erp_post_scalp(idx_right_cluster, :), 1);
    h_left_post  = plot(t_vec, erp_post_left,  'Color', [0.85 0.15 0.15], 'LineWidth', 2.0, 'DisplayName', 'Left Eye');
    h_right_post = plot(t_vec, erp_post_right, 'Color', [0.15 0.65 0.20], 'LineWidth', 1.8, 'DisplayName', 'Right Eye');
    grid on; xlim([min(t_vec), max(t_vec)]); 
    % ylim(y_lims);
    ylim([-15 15]);
    xlabel('Time from Blink Peak (ms)', 'FontSize', 9);
    ylabel('Residual Voltage (\muV)', 'FontSize', 9);
    title(sprintf('Post-Pruning Clean ERP (N = %d, Peak = %.1f \\muV)', ...
        length(clean_snip_idx), peak_post_fp), 'FontWeight', 'bold', 'FontSize', 10);
    legend([h_left_post, h_right_post], 'Location', 'south', 'FontSize', 8);

    % Row 2, Col 2: Post-Pruning Topography at t = 0 ms
    sbh5 = subplot(2, 3, 5);
    render_topoplot(topo_post, 'Post-Pruning Clean Residual (t = 0 ms)', sbh5, topo_lims, '\muV');

    % Row 2, Col 3: Retained Epochs Average Power
    sbh6 = subplot(2, 3, 6);
    render_topoplot(topo_power_clean, sprintf('Retained Power (N = %d)', n_clean_trials), ...
        sbh6, power_lims, '\muV^2');

    % Save figure if pipeline paths exist
    save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_leftovers_2_' num2str(tag)], [28 18]);
end

% -------------------------------------------------------------------------
% 9. Log Statistics to Dataset Struct
% -------------------------------------------------------------------------
DATA.ALSUTRECHT.leftovers.blink_stats = blink_stats;

% -------------------------------------------------------------------------
% 10. Prune Compromised Blink Epochs from Dataset (If Epoched)
% -------------------------------------------------------------------------
if is_epoched && ~isempty(bad_blink_trials)
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
elseif is_epoched
    fprintf('No epochs required pruning; all %d epochs retained.\n', DATA.trials);
end

end

% =========================================================================
% Local Helper: Topoplot Dispatcher with Dynamic Colourbar Unit
% =========================================================================
function render_topoplot(values, plot_title, ax_handle, map_limits, unit_str)
mytopoplot(values, false(size(values)), plot_title, ax_handle, map_limits);
cb = colorbar('peer', ax_handle);
ylabel(cb, unit_str, 'FontSize', 9, 'FontWeight', 'bold');
end