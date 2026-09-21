function DATA = detect_blink_residuals(DATA, cfg)
% DETECT_BLINK_RESIDUALS Evaluates residual ocular leakage on the blink-locked
% average ERP. Identifies systematic leakage (>10 uV) across both continuous
% (2D) and epoched (3D) EEGLAB datasets without dropping edge blinks.
%
% Inputs:
%    DATA : EEGLAB dataset structure (continuous 2D or epoched 3D)
%    cfg  : Optional configuration structure:
%           * cfg.veog_chan      : VEOG channel label or index (default: 'VEOG')
%           * cfg.trIQRblink     : IQR multiplier for VEOG detection (default: 3.0)
%           * cfg.min_phys_blink : Minimum VEOG deflection in uV (default: 50 uV)
%           * cfg.max_avg_res_uV : Max tolerable average ERP residual peak (default: 10 uV)
%           * cfg.max_avg_p2p_uV : Max tolerable average ERP peak-to-peak (default: 15 uV)
%           * cfg.win_sec        : Half-window duration around peak (default: 0.400 s)
%           * cfg.isolate_blinks : Discard overlapping peaks in continuous (default: true)
%           * cfg.do_plot        : Render 2-panel diagnostic figure (default: true)
%
% Outputs:
%    blink_info : Struct containing average ERP waveforms, peak residual amplitude,
%                 snippet arrays, and outlier trial indices.

% -------------------------------------------------------------------------
% 1. Parameter Validation and Defaults
% -------------------------------------------------------------------------
if nargin < 2 || isempty(cfg), cfg = struct(); end

if ~isfield(cfg, 'veog_chan'),      cfg.veog_chan      = 'VEOG'; end
if ~isfield(cfg, 'trIQRblink'),     cfg.trIQRblink     = 3.0;    end
if ~isfield(cfg, 'min_phys_blink'), cfg.min_phys_blink = 50;     end % uV
if ~isfield(cfg, 'max_avg_res_uV'), cfg.max_avg_res_uV = 10;     end % uV (Target threshold)
if ~isfield(cfg, 'max_avg_p2p_uV'), cfg.max_avg_p2p_uV = 15;     end % uV
if ~isfield(cfg, 'win_sec'),        cfg.win_sec        = 0.400;  end % s (+/- 400 ms)
if ~isfield(cfg, 'isolate_blinks'), cfg.isolate_blinks = true;   end
if ~isfield(cfg, 'do_plot'),        cfg.do_plot        = true;   end

fs         = DATA.srate;
is_epoched = (ndims(DATA.data) == 3) && (size(DATA.data, 3) > 1);

% -------------------------------------------------------------------------
% 2. Resolve EEG and VEOG Channels
% -------------------------------------------------------------------------
ch_idx = find(strcmpi({DATA.chanlocs.type}, 'EEG'));
idx_veog = find(strcmpi({DATA.chanlocs.labels}, cfg.veog_chan), 1);

assert(~isempty(idx_veog), 'VEOG channel [%s] not found in DATA.', string(cfg.veog_chan));
ch_idx = setdiff(ch_idx, idx_veog);
n_chans = length(ch_idx);

% Resolve fronto-polar cluster (BioSemi 128 and 10-20 standards)
labels = {DATA.chanlocs(ch_idx).labels};
left_eye_chans  = {'C28', 'C29', 'C30', 'C31', 'C17', 'C18', 'FP1', 'AF7'};
right_eye_chans = {'C8', 'C9', 'C15', 'C16', 'C17', 'C18', 'FP2', 'AF8'};

idx_left_cluster  = find(ismember(upper(labels), upper(left_eye_chans)));
idx_right_cluster = find(ismember(upper(labels), upper(right_eye_chans)));
ocular_indices    = unique([idx_left_cluster, idx_right_cluster]);
assert(~isempty(ocular_indices), 'Could not resolve fronto-polar channel clusters.');

% -------------------------------------------------------------------------
% 3. Timing and Snippet Windows
% -------------------------------------------------------------------------
win_pts   = round(cfg.win_sec * fs);
t_vec     = (-win_pts : win_pts) * (1000 / fs); % Time in ms
n_win_pts = length(t_vec);
[~, t_zero_idx] = min(abs(t_vec));

base_idx   = find(t_vec <= -250);
search_idx = find(t_vec >= -125 & t_vec <= 175);
if isempty(search_idx), search_idx = t_zero_idx; end

% -------------------------------------------------------------------------
% 4. Blink Event Extraction (Edge-Safe with NaN Padding)
% -------------------------------------------------------------------------
all_snippets_scalp = [];
all_snippets_veog  = [];
snippet_trials     = [];
n_snip             = 0;

if is_epoched
    % ======================== Epoched Branch =============================
    n_trials = size(DATA.data, 3);
    n_pnts   = size(DATA.data, 2);
    veog_3d  = double(squeeze(DATA.data(idx_veog, :, :)));

    veog_flat    = veog_3d(:)';
    veog_detrend = veog_flat - movmedian(veog_flat, fs);
    thresh_veog  = max(cfg.min_phys_blink, prctile(veog_detrend, 75) + (cfg.trIQRblink * iqr(veog_detrend)));

    blink_epochs = false(1, n_trials);

    for tr = 1:n_trials
        v_tr = veog_3d(:, tr);
        v_dt = v_tr - median(v_tr);

        is_blink_sample = abs(v_dt) >= thresh_veog;

        if any(is_blink_sample)
            blink_epochs(tr) = true;

            % Baseline calculated exclusively from non-blink samples of this trial
            clean_samples = ~is_blink_sample;
            if sum(clean_samples) >= round(0.150 * fs)
                epoch_base = mean(DATA.data(ch_idx, clean_samples, tr), 2);
                v_base     = median(v_tr(clean_samples));
            else
                epoch_base = median(DATA.data(ch_idx, :, tr), 2);
                v_base     = median(v_tr);
            end

            [~, t_peak] = max(abs(v_dt));

            % Bounded window mapping with NaN padding for edge blinks
            t_start = t_peak - win_pts;
            t_end   = t_peak + win_pts;

            src_start  = max(1, t_start);
            src_end    = min(n_pnts, t_end);
            dest_start = src_start - t_start + 1;
            dest_end   = dest_start + (src_end - src_start);

            % Populate snippet
            snip_scalp = NaN(n_chans, n_win_pts);
            snip_scalp(:, dest_start:dest_end) = double(DATA.data(ch_idx, src_start:src_end, tr)) - epoch_base;

            snip_veog = NaN(1, n_win_pts);
            snip_veog(dest_start:dest_end) = v_tr(src_start:src_end) - v_base;

            n_snip = n_snip + 1;
            all_snippets_scalp(:, :, n_snip) = snip_scalp; %#ok<AGROW>
            all_snippets_veog(:, :, n_snip)  = snip_veog;  %#ok<AGROW>
            snippet_trials(n_snip)           = tr;         %#ok<AGROW>
        end
    end

else
    % ======================= Continuous Branch ===========================
    n_pnts       = size(DATA.data, 2);
    veog_raw     = double(DATA.data(idx_veog, :));
    veog_detrend = veog_raw - movmedian(veog_raw, fs);

    thresh_veog = max(cfg.min_phys_blink, prctile(veog_detrend, 75) + (cfg.trIQRblink * iqr(veog_detrend)));
    min_dist    = round(0.35 * fs);

    [~, all_peaks] = findpeaks(veog_detrend, 'MinPeakHeight', thresh_veog, 'MinPeakDistance', min_dist);

    if cfg.isolate_blinks && length(all_peaks) > 1
        isolation_dist = round(cfg.win_sec * 2.0 * fs);
        is_isolated    = true(size(all_peaks));
        for i_p = 1:length(all_peaks)
            seps = abs(all_peaks - all_peaks(i_p));
            seps(i_p) = Inf;
            if any(seps < isolation_dist)
                is_isolated(i_p) = false;
            end
        end
        valid_peaks = all_peaks(is_isolated);
    else
        valid_peaks = all_peaks;
    end

    blink_epochs = valid_peaks;

    for i_b = 1:length(valid_peaks)
        pnt     = valid_peaks(i_b);
        t_start = pnt - win_pts;
        t_end   = pnt + win_pts;

        src_start  = max(1, t_start);
        src_end    = min(n_pnts, t_end);
        dest_start = src_start - t_start + 1;
        dest_end   = dest_start + (src_end - src_start);

        snip_scalp = NaN(n_chans, n_win_pts);
        raw_snip   = double(DATA.data(ch_idx, src_start:src_end));

        % Continuous baseline from pre-blink samples
        valid_base_idx = base_idx(base_idx >= dest_start & base_idx <= dest_end);
        if ~isempty(valid_base_idx)
            base_val = mean(raw_snip(:, valid_base_idx - dest_start + 1), 2);
            v_base   = mean(veog_raw(src_start + valid_base_idx(1) - dest_start));
        else
            base_val = median(raw_snip, 2);
            v_base   = median(veog_raw(src_start:src_end));
        end

        snip_scalp(:, dest_start:dest_end) = raw_snip - base_val;

        snip_veog = NaN(1, n_win_pts);
        snip_veog(dest_start:dest_end) = veog_raw(src_start:src_end) - v_base;

        n_snip = n_snip + 1;
        all_snippets_scalp(:, :, n_snip) = snip_scalp; %#ok<AGROW>
        all_snippets_veog(:, :, n_snip)  = snip_veog;  %#ok<AGROW>
        snippet_trials(n_snip)           = i_b;        %#ok<AGROW>
    end
end

% Identify and discard any snippets containing NaNs in either VEOG or Scalp
nan_mask = squeeze(any(isnan(all_snippets_scalp), [1, 2])) | squeeze(any(isnan(all_snippets_veog), [1, 2]));

if any(nan_mask)
    fprintf('Discarded %d edge-truncated snippet(s) containing NaNs.\n', sum(nan_mask));
    all_snippets_scalp(:, :, nan_mask) = [];
    all_snippets_veog(:, :, nan_mask)  = [];
    snippet_trials(nan_mask)           = [];
end

% Now the count represents genuine, complete snippets
n_snippets = size(all_snippets_scalp, 3);

% -------------------------------------------------------------------------
% 5. Evaluate Residual on the AVERAGED Blink ERP (Omit NaNs)
% -------------------------------------------------------------------------
if n_snippets > 0
    erp_scalp = mean(all_snippets_scalp, 3, 'omitnan');
    erp_veog  = mean(all_snippets_veog, 3, 'omitnan');

    % Peak residual evaluated on the fronto-polar cluster of the average waveform
    fp_erp = erp_scalp(ocular_indices, search_idx);
    max_avg_res_uV = max(abs(fp_erp(:)));
    max_avg_p2p_uV = max(max(fp_erp, [], 2) - min(fp_erp, [], 2));

    % Dataset fails only if the grand-average fronto-polar ERP exceeds threshold
    has_residue = (max_avg_res_uV > cfg.max_avg_res_uV) || (max_avg_p2p_uV > cfg.max_avg_p2p_uV);
else
    erp_scalp      = zeros(n_chans, n_win_pts);
    erp_veog       = zeros(1, n_win_pts);
    max_avg_res_uV = 0;
    max_avg_p2p_uV = 0;
    has_residue    = false;
end

% -------------------------------------------------------------------------
% 6. Trial Breakdown
% -------------------------------------------------------------------------
bad_blink_trials = [];
n_snip = 0;

% Only search for bad trials if the AVERAGE failed the 10 uV test
if has_residue && is_epoched && n_snippets > 0
    veog_template = erp_veog(search_idx);
    veog_template = veog_template / max(abs(veog_template));

    for i_s = 1:n_snippets
        tr_id = snippet_trials(i_s);
        snip  = all_snippets_scalp(ocular_indices, search_idx, i_s);

        % Ignore if snippet has missing values in the search window
        if any(isnan(snip(:))), continue; end

        tr_fp = mean(snip, 1);
        % Must show systematic positive deflection matching the blink polarity
        if (max(abs(tr_fp)) > (2 * cfg.max_avg_res_uV)) && (dot(tr_fp, veog_template) > 0)

            n_snip = n_snip + 1;
            bad_blink_trials(n_snip) = tr_id; %#ok<AGROW>
        end
    end
    bad_blink_trials = unique(bad_blink_trials);
end

if is_epoched
    n_total_trials = size(DATA.data, 3);
    clean_trials   = setdiff(1:n_total_trials, bad_blink_trials);
else
    clean_trials   = setdiff(1:length(blink_epochs), bad_blink_trials);
end

% -------------------------------------------------------------------------
% 7. Package Diagnostics Output
% -------------------------------------------------------------------------
blink_info.is_epoched         = is_epoched;
blink_info.n_snippets         = n_snippets;
blink_info.all_snippets_scalp = all_snippets_scalp;
blink_info.all_snippets_veog  = all_snippets_veog;
blink_info.erp_scalp          = erp_scalp;
blink_info.erp_veog           = erp_veog;
blink_info.max_avg_res_uV     = max_avg_res_uV;
blink_info.max_avg_p2p_uV     = max_avg_p2p_uV;
blink_info.has_residue        = has_residue;
blink_info.bad_blink_trials   = bad_blink_trials;
blink_info.clean_trials       = clean_trials;
blink_info.t_vec              = t_vec;
blink_info.t_zero_idx         = t_zero_idx;
blink_info.ch_idx             = ch_idx;
blink_info.ocular_indices     = ocular_indices;
blink_info.idx_left_cluster   = idx_left_cluster;
blink_info.idx_right_cluster  = idx_right_cluster;
blink_info.snippet_trials     = snippet_trials;
blink_info.blink_epochs       = blink_epochs;
blink_info.thresh_veog        = thresh_veog;

% Log
DATA.ALSUTRECHT.leftovers.blink_info = blink_info;

res_status = 'PASS (Clean)';
if has_residue
    res_status = 'FAIL (Residual > 10 uV)';
end

fprintf('Blink Audit (Average ERP): N = %d snippets | Peak Fronto-Polar = %.2f uV (Limit: %.1f uV) -> %s\n', ...
    n_snippets, max_avg_res_uV, cfg.max_avg_res_uV, res_status);

% -------------------------------------------------------------------------
% 8. Diagnostic Figure Visualisation
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 8. Diagnostic Figure Visualisation (Butterfly + 3 Time Topographies)
% -------------------------------------------------------------------------
if cfg.do_plot && n_snippets > 0
    fh = figure('Color', 'w', 'Position', [80 100 1500 380], 'Name', 'Blink ERP Residual Audit');
    tlo = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    line_vals = mean(erp_scalp(ocular_indices, :), 1);
    c_lim_erp = max(ceil(max(abs(line_vals))), cfg.max_avg_res_uV) + 2;

    % Panel 1: Scalp Butterfly ERP with Fronto-Polar Cluster
    ax1 = nexttile(tlo, 1); hold(ax1, 'on');
    plot(ax1, t_vec, erp_scalp', 'Color', [0.80 0.80 0.80 0.7], 'LineWidth', 0.8);
    h_fp = plot(ax1, t_vec, line_vals, 'Color', [0.85 0.20 0.10], ...
        'LineWidth', 2.0, 'DisplayName', 'Frontal Average');

    yline(ax1,  cfg.max_avg_res_uV, 'k--', 'LineWidth', 1.2, 'DisplayName', sprintf('+%.0f \\muV Limit', cfg.max_avg_res_uV));
    yline(ax1, -cfg.max_avg_res_uV, 'k--', 'LineWidth', 1.2, 'DisplayName', sprintf('-%.0f \\muV Limit', cfg.max_avg_res_uV));
    xline(ax1, -250, ':', 'Color', [0.3 0.3 0.8], 'LineWidth', 1.2);
    pbaspect(ax1, [1.618 1 1]);

    % Topoplot evaluation times
    topo_times = [0, 100, 220]; % ms (Blink peak, Re-opening, Late deflection)
    topo_colors = {[0.2 0.6 0.2], [0.8 0.5 0.1], [0.5 0.2 0.7]};
    for it = 1:length(topo_times)
        xline(ax1, topo_times(it), '--', 'Color', topo_colors{it}, 'LineWidth', 1.3);
    end

    xlabel(ax1, 'Time relative to blink peak (ms)', 'FontSize', 8);
    ylabel(ax1, 'Average Amplitude (\muV)', 'FontSize', 8);
    title(ax1, sprintf('Blink-Locked ERP (N = %d) | Peak = %.2f \\muV', n_snippets, max_avg_res_uV), ...
        'FontSize', 8, 'FontWeight', 'bold');
    box(ax1, 'off'); grid(ax1, 'on'); ylim(ax1, [-c_lim_erp c_lim_erp]);
    legend(h_fp, 'Location', 'northeast', 'Box', 'off');

    % Find sample indices closest to target latencies
    [~, topo_indices] = min(abs(t_vec - topo_times(:)), [], 2);

    % Calculate unified symmetric color limits across the 3 topographies
    all_topo_vals = erp_scalp(:, topo_indices);
    c_lim_topo = max(ceil(max(abs(all_topo_vals(:)))), cfg.max_avg_res_uV);
    topo_clim = [-c_lim_topo, c_lim_topo];

    % Panels 2-4: Scalp Topographies across time
    for it = 1:length(topo_times)
        ax = nexttile(tlo, it + 1);
        t_sample = topo_indices(it);
        topo_vals = erp_scalp(:, t_sample);

        mytopoplot(topo_vals, false(size(topo_vals)), '', ax, topo_clim);

        cb = colorbar(ax);
        ylabel(cb, '\muV', 'FontSize', 9, 'FontWeight', 'bold');
        title(ax, sprintf('t = %.0f ms', t_vec(t_sample)), ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', topo_colors{it});
    end

    % Save
    save_figure(fh, DATA.ALSUTRECHT.subject.figures, [DATA.ALSUTRECHT.subject.id '_leftovers_2_0'], [40 12]);
end
end