function [cohort_table, benchmark_stats] = benchmark_ica_cohort(data_dir, cfg)
% BENCHMARK_ICA_COHORT Loads individual subject *_ica_class.mat files, extracts
% classifications from ica.final logical masks, maps ICLabel predictions,
% harmonises schema differences for participants without ECG channels, aggregates
% recording-level variance, and computes cluster separability and spatial coherence.

if nargin < 2, cfg = struct(); end
if ~isfield(cfg, 'outdir'),       cfg.outdir       = data_dir; end
if ~isfield(cfg, 'do_plot'),      cfg.do_plot      = true;     end
if ~isfield(cfg, 'embedding'),    cfg.embedding    = 'tsne';   end
if ~isfield(cfg, 'kurtosis_cap'), cfg.kurtosis_cap = 100;      end

file_records = dir(fullfile(data_dir, '*_ica_class.mat'));
n_files = length(file_records);
assert(n_files > 0, 'No "*_ica_class.mat" files found in %s', data_dir);

fprintf('=== Benchmarking ICA Classification Across %d Files ===\n', n_files);

all_records     = cell(n_files, 1);
all_winvs       = cell(n_files, 1);
all_iclabels    = cell(n_files, 1);
all_var_records = cell(n_files, 1);

% -------------------------------------------------------------------------
% 1. Ingestion: Harvest Labels, IC Features, and Subject-Level Variance
% -------------------------------------------------------------------------
for i_sub = 1:n_files
    file_path = fullfile(file_records(i_sub).folder, file_records(i_sub).name);
    sub_id    = erase(file_records(i_sub).name, '_ica_class.mat');

    loaded = load(file_path);
    if isfield(loaded, 'ica_class')
        ica = loaded.ica_class;
    elseif isfield(loaded, 'ica')
        ica = loaded.ica;
    else
        continue;
    end

    n_ics = size(ica.icawinv, 2);
    if isempty(n_ics) || n_ics == 0, continue; end

    % --- Class Resolution via ica.final ---
    class_label = repmat("Neural", n_ics, 1);
    if isfield(ica, 'final') && isstruct(ica.final)
        fin = ica.final;
        for k = 1:n_ics
            if safe_bool(fin, 'heart', k)
                class_label(k) = "Heart";
            elseif safe_bool(fin, 'eye', k)
                class_label(k) = "Eye";
            elseif safe_bool(fin, 'muscle', k)
                class_label(k) = "Muscle";
            elseif safe_bool(fin, 'channel', k)
                class_label(k) = "Channel";
            elseif safe_bool(fin, 'complex', k)
                class_label(k) = "Complex";
            elseif safe_bool(fin, 'genbad', k)
                class_label(k) = "GenBad";
            end
        end
    end

    % --- Recording-Level Variance Summary (ica.final.var) ---
    if isfield(ica, 'final') && isfield(ica.final, 'var') && isstruct(ica.final.var)
        fv = ica.final.var;
        sub_var = table(...
            string(sub_id), ...
            safe_scalar(fv, 'brain'), ...
            safe_scalar(fv, 'muscle'), ...
            safe_scalar(fv, 'eye'), ...
            safe_scalar(fv, 'heart'), ...
            safe_scalar(fv, 'channel'), ...
            safe_scalar(fv, 'line'), ...
            safe_scalar(fv, 'other'), ...
            safe_scalar(fv, 'genbad'), ...
            'VariableNames', {'Subject', 'Brain', 'Muscle', 'Eye', 'Heart', 'Channel', 'Line', 'Other', 'GenBad'});
        all_var_records{i_sub} = sub_var;
    end

    % --- Per-IC Feature Harvesting ---
    feat_struct = struct();

    % Root-level variance metrics
    feat_struct.vaf_compvar   = extract_vec(safe_field(ica, 'vaf_compvar'), n_ics);
    feat_struct.var_ica       = extract_vec(safe_field(ica, 'var_ica'), n_ics);
    feat_struct.var_peak_chan = extract_vec(safe_field(ica, 'var_peak_chan'), n_ics);

    % Spatial metrics from mixing matrix icawinv
    winv = ica.icawinv;
    feat_struct.spatial_ipr = ((sum(winv.^2, 1).^2) ./ sum(winv.^4, 1))';
    all_winvs{i_sub} = winv;

    % Harvest domain sub-structs (corr, heart, muscle, blink, saccade, channel, bad)
    domains = {'corr', 'heart', 'muscle', 'blink', 'saccade', 'channel', 'bad'};
    for d = 1:length(domains)
        dom_name = domains{d};
        if isfield(ica, dom_name) && isstruct(ica.(dom_name))
            sub_feats = harvest_numeric_vectors(ica.(dom_name), [dom_name '_'], n_ics);
            f_names = fieldnames(sub_feats);
            for f = 1:length(f_names)
                feat_struct.(f_names{f}) = sub_feats.(f_names{f});
            end
        end
    end

    % Harvest ICLabel results (classes, cvec, pvec)
    iclabel_top_class = repmat("None", n_ics, 1);
    if isfield(ica, 'ICLabel') && isstruct(ica.ICLabel)
        icl = ica.ICLabel;
        if isfield(icl, 'cvec') && ~isempty(icl.cvec) && isfield(icl, 'classes') && ~isempty(icl.classes)
            icl_names = string(icl.classes);
            c_indices = round(double(icl.cvec(:)));
            valid_c   = c_indices >= 1 & c_indices <= length(icl_names);
            iclabel_top_class(valid_c) = icl_names(c_indices(valid_c));
        end
        if isfield(icl, 'pvec') && ~isempty(icl.pvec)
            feat_struct.ICL_Confidence = extract_vec(icl.pvec, n_ics);
        end
    end
    all_iclabels{i_sub} = iclabel_top_class;

    % Assemble single subject table
    sub_table = struct2table(feat_struct);
    sub_table.Subject       = repmat(string(sub_id), n_ics, 1);
    sub_table.IC            = (1:n_ics)';
    sub_table.AssignedClass = categorical(class_label);
    sub_table.ICLabel_Top   = categorical(iclabel_top_class);

    all_records{i_sub} = sub_table;
end

% -------------------------------------------------------------------------
% 1b. Schema Harmonisation across Varied Acquisition Montages
% -------------------------------------------------------------------------
% Find union of all feature columns across subjects
all_var_names = {};
for i_sub = 1:n_files
    if ~isempty(all_records{i_sub})
        all_var_names = union(all_var_names, all_records{i_sub}.Properties.VariableNames, 'stable');
    end
end

% Pad absent columns with NaN so subjects without ECG can be merged
for i_sub = 1:n_files
    if isempty(all_records{i_sub}), continue; end
    cur_t = all_records{i_sub};
    missing_vars = setdiff(all_var_names, cur_t.Properties.VariableNames, 'stable');
    for mv = 1:length(missing_vars)
        var_name = missing_vars{mv};
        if strcmp(var_name, 'AssignedClass')
            cur_t.(var_name) = repmat(categorical("Neural"), height(cur_t), 1);
        elseif strcmp(var_name, 'ICLabel_Top')
            cur_t.(var_name) = repmat(categorical("None"), height(cur_t), 1);
        else
            cur_t.(var_name) = nan(height(cur_t), 1);
        end
    end
    all_records{i_sub} = cur_t(:, all_var_names);
end

cohort_table = vertcat(all_records{:});
cohort_var_table = vertcat(all_var_records{:});
assert(~isempty(cohort_table), 'No valid records extracted from input files.');

% -------------------------------------------------------------------------
% 2. Dynamic Feature Audit & Missing Data Handling
% -------------------------------------------------------------------------
fprintf('\n--- Auditing Available Quantitative Metrics ---\n');

meta_cols = {'Subject', 'IC', 'AssignedClass', 'ICLabel_Top'};
cand_features = setdiff(cohort_table.Properties.VariableNames, meta_cols, 'stable');

active_features = {};
for i_f = 1:length(cand_features)
    f_name  = cand_features{i_f};
    col_val = cohort_table.(f_name);

    if ~isnumeric(col_val), continue; end

    nan_pct = mean(isnan(col_val)) * 100;
    val_var = var(col_val(~isnan(col_val)));

    if nan_pct == 100 || isempty(val_var) || val_var == 0
        fprintf('  [DISABLED]   %-28s: 100%% missing / invariant -> EXCLUDED\n', f_name);
    elseif contains(lower(f_name), {'ecg', 'heart'})
        % Explicitly handle ECG metrics where channels were uncollected or sparse
        fprintf('  [ECG-MODAL]  %-28s: %5.1f%% missing (missing lead/candidate) -> INCLUDED (0-imputed)\n', f_name, nan_pct);
        cohort_table.(f_name)(isnan(col_val)) = 0;
        active_features{end+1} = f_name; %#ok<AGROW>
    elseif nan_pct > 30
        % Other sparse candidate metrics (e.g. blink peak residual deflections)
        fprintf('  [SPARSE]     %-28s: %5.1f%% missing -> INCLUDED (0-imputed)\n', f_name, nan_pct);
        cohort_table.(f_name)(isnan(col_val)) = 0;
        active_features{end+1} = f_name; %#ok<AGROW>
    else
        % Dense feature active across subjects
        if nan_pct > 0
            fprintf('  [DENSE]      %-28s: %5.1f%% missing -> INCLUDED (median imputed)\n', f_name, nan_pct);
            cohort_table.(f_name)(isnan(col_val)) = median(col_val(~isnan(col_val)));
        else
            fprintf('  [DENSE]      %-28s: 100.0%% populated -> INCLUDED\n', f_name);
        end
        active_features{end+1} = f_name; %#ok<AGROW>
    end
end

assert(~isempty(active_features), 'No populated numerical features detected.');

feature_mat = cohort_table{:, active_features};

% Cap kurtosis if active
kurt_idx = find(contains(lower(active_features), 'kurtosis'));
if ~isempty(kurt_idx)
    for k = 1:length(kurt_idx)
        feature_mat(:, kurt_idx(k)) = min(feature_mat(:, kurt_idx(k)), cfg.kurtosis_cap);
    end
end

% -------------------------------------------------------------------------
% 3. Cluster Separability (Silhouette Analysis)
% -------------------------------------------------------------------------
labels     = cohort_table.AssignedClass;
features_z = zscore(feature_mat);

[sil_vals, ~] = silhouette(features_z, labels, 'Euclidean');
cohort_table.Silhouette = sil_vals;

classes = categories(labels);
mean_sil_per_class = zeros(length(classes), 1);
fprintf('\n--- Cluster Separation by Class (ica.final) ---\n');
for c = 1:length(classes)
    c_mask = labels == classes{c};
    mean_sil_per_class(c) = mean(sil_vals(c_mask));
    fprintf('  Class [%-12s]: N = %4d | Mean Silhouette = %+0.3f\n', ...
        classes{c}, sum(c_mask), mean_sil_per_class(c));
end
fprintf('  Overall Mean Silhouette: %+0.3f\n', mean(sil_vals));

if ~isempty(cohort_var_table)
    fprintf('\n--- Mean Variance Accounted For Across Cohort (%%) ---\n');
    fprintf('  Brain:   %5.1f ± %4.1f%%\n', mean(cohort_var_table.Brain, 'omitnan'), std(cohort_var_table.Brain, 'omitnan'));
    fprintf('  GenBad:  %5.1f ± %4.1f%%\n', mean(cohort_var_table.GenBad, 'omitnan'), std(cohort_var_table.GenBad, 'omitnan'));
    fprintf('  Muscle:  %5.1f ± %4.1f%%\n', mean(cohort_var_table.Muscle, 'omitnan'), std(cohort_var_table.Muscle, 'omitnan'));
    fprintf('  Eye:     %5.1f ± %4.1f%%\n', mean(cohort_var_table.Eye, 'omitnan'), std(cohort_var_table.Eye, 'omitnan'));
    fprintf('  Heart:   %5.1f ± %4.1f%%\n', mean(cohort_var_table.Heart, 'omitnan'), std(cohort_var_table.Heart, 'omitnan'));
    fprintf('  Channel: %5.1f ± %4.1f%%\n', mean(cohort_var_table.Channel, 'omitnan'), std(cohort_var_table.Channel, 'omitnan'));
end

% -------------------------------------------------------------------------
% 4. Dimensionality Reduction (t-SNE / PCA)
% -------------------------------------------------------------------------
fprintf('\n--- Computing 2D %s projection (%d features) ---\n', upper(cfg.embedding), length(active_features));
if strcmpi(cfg.embedding, 'tsne')
    coords_2d = tsne(features_z, 'NumDimensions', 2, 'Perplexity', 30, 'Standardize', false);
else
    [~, score] = pca(features_z);
    coords_2d = score(:, 1:2);
end
cohort_table.Dim1 = coords_2d(:, 1);
cohort_table.Dim2 = coords_2d(:, 2);

% -------------------------------------------------------------------------
% 5. Spatial Coherence Analysis
% -------------------------------------------------------------------------
spatial_stats = struct('class', {}, 'within_r_mean', {}, 'between_r_mean', {});
pool_w = horzcat(all_winvs{:});

if ~isempty(pool_w)
    w_norm = pool_w ./ sqrt(sum(pool_w.^2, 1) + eps);
    topo_labels = cohort_table.AssignedClass;

    for c = 1:length(classes)
        c_name  = classes{c};
        idx_in  = find(topo_labels == c_name);
        idx_out = find(topo_labels ~= c_name);

        if length(idx_in) > 250,  idx_in  = idx_in(randperm(length(idx_in), 250)); end
        if length(idx_out) > 250, idx_out = idx_out(randperm(length(idx_out), 250)); end

        if length(idx_in) >= 2
            R_within = abs(w_norm(:, idx_in)' * w_norm(:, idx_in));
            tri_idx  = triu(true(size(R_within)), 1);
            within_r = mean(R_within(tri_idx));

            if ~isempty(idx_out)
                R_between = abs(w_norm(:, idx_in)' * w_norm(:, idx_out));
                between_r = mean(R_between(:));
            else
                between_r = NaN;
            end

            spatial_stats(c).class          = c_name;
            spatial_stats(c).within_r_mean  = within_r;
            spatial_stats(c).between_r_mean = between_r;
        end
    end
end

% Pack outputs
benchmark_stats.active_features          = active_features;
benchmark_stats.mean_silhouette_overall  = mean(sil_vals);
benchmark_stats.mean_silhouette_by_class = table(categorical(classes), mean_sil_per_class, ...
    'VariableNames', {'Class', 'MeanSilhouette'});
benchmark_stats.spatial_stats            = spatial_stats;
benchmark_stats.variance_summary         = cohort_var_table;

% -------------------------------------------------------------------------
% 6. Visualisation (2x2 Dashboard)
% -------------------------------------------------------------------------
if ~cfg.do_plot, return; end

fh = figure('Color', 'w', 'Position', [80, 80, 1250, 800], 'Name', 'ICA Cohort Benchmark');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

n_classes = length(classes);
if n_classes <= 8
    cmap = brewermap(max(3, n_classes), 'Set2');
else
    cmap = lines(n_classes);
end

% Panel 1: 2D Embedding Space
ax1 = nexttile(t, 1);
hold(ax1, 'on');
for c = 1:n_classes
    c_idx = cohort_table.AssignedClass == classes{c};
    scatter(ax1, cohort_table.Dim1(c_idx), cohort_table.Dim2(c_idx), 28, ...
        cmap(c, :), 'filled', 'MarkerFaceAlpha', 0.65, 'DisplayName', char(classes{c}));
end
grid(ax1, 'on');
set(ax1, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(ax1, sprintf('%s Dimension 1', upper(cfg.embedding)), 'FontWeight', 'bold');
ylabel(ax1, sprintf('%s Dimension 2', upper(cfg.embedding)), 'FontWeight', 'bold');
title(ax1, sprintf('Decision Space Embedding (%d Active Features)', length(active_features)), ...
    'FontSize', 11, 'FontWeight', 'bold');
legend(ax1, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
hold(ax1, 'off');

% Panel 2: Silhouette Compactness per Class
ax2 = nexttile(t, 2);
hold(ax2, 'on');
y_offset   = 0;
yticks_vec = zeros(n_classes, 1);

for c = 1:n_classes
    c_sils = sort(sil_vals(labels == classes{c}), 'descend');
    n_k    = length(c_sils);
    if n_k > 0
        y_range = (1:n_k)' + y_offset;
        fill(ax2, [0; c_sils; 0], [y_range(1); y_range; y_range(end)], cmap(c, :), ...
            'FaceAlpha', 0.7, 'EdgeColor', 'none');
        yticks_vec(c) = y_offset + (n_k / 2);
        y_offset      = y_offset + n_k + 30;
    end
end
xline(ax2, 0, 'k--', 'LineWidth', 1);
xline(ax2, mean(sil_vals), 'r:', 'LineWidth', 1.2, 'DisplayName', 'Cohort Mean');
xlim(ax2, [-0.5, 1.0]);
ylim(ax2, [0, y_offset]);
grid(ax2, 'on');
set(ax2, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 9, 'YTick', yticks_vec, 'YTickLabel', classes);
xlabel(ax2, 'Silhouette Coefficient s(i)', 'FontWeight', 'bold');
ylabel(ax2, 'Assigned Class Clusters', 'FontWeight', 'bold');
title(ax2, sprintf('Cluster Compactness (Global Mean: %.2f)', mean(sil_vals)), ...
    'FontSize', 11, 'FontWeight', 'bold');
hold(ax2, 'off');

% Panel 3: Spatial Coherence Contrast (|r| within vs between)
ax3 = nexttile(t, 3);
if ~isempty(spatial_stats)
    c_names     = {spatial_stats.class};
    within_vec  = [spatial_stats.within_r_mean];
    between_vec = [spatial_stats.between_r_mean];

    b = bar(ax3, [within_vec; between_vec]', 'grouped');
    b(1).FaceColor = [0.20, 0.45, 0.70];
    b(2).FaceColor = [0.85, 0.30, 0.30];

    set(ax3, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 9, ...
        'XTick', 1:length(c_names), 'XTickLabel', c_names);
    xtickangle(ax3, 30);
    ylabel(ax3, 'Mean Absolute Correlation |r|', 'FontWeight', 'bold');
    ylim(ax3, [0, 1.0]);
    grid(ax3, 'on');
    legend(ax3, {'Within Class', 'Between Classes'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    title(ax3, 'Spatial Coherence from icawinv', 'FontSize', 11, 'FontWeight', 'bold');
end

% Panel 4: Cross-Method Agreement (ica.final vs ICLabel)
has_valid_icl = any(cohort_table.ICLabel_Top ~= "None");

if has_valid_icl
    % 1. Map ICLabel terminology while preserving 'Other' as a distinct category
    icl_mapped = string(cohort_table.ICLabel_Top);
    icl_mapped(icl_mapped == "Brain")         = "Neural";
    icl_mapped(icl_mapped == "Channel Noise") = "Channel";
    icl_mapped(icl_mapped == "Line Noise")    = "GenBad";
    % Do NOT map "Other" to "GenBad" -> keep it as "Other"

    % 2. Define category structure
    % Pipeline rows
    row_classes = ["Neural", "Muscle", "Eye", "Heart", "Channel", "GenBad", "Complex"];
    % ICLabel columns include "Other"
    all_classes = ["Neural", "Muscle", "Eye", "Heart", "Channel", "GenBad", "Complex", "Other"];

    pipe_cats = categorical(cohort_table.AssignedClass, all_classes);
    icl_cats  = categorical(icl_mapped, all_classes);

    % 3. Render confusion chart
    nexttile(4);
    cm = confusionchart(t, pipe_cats, icl_cats, ...
        'RowSummary', 'row-normalized', ...
        'ColumnSummary', 'column-normalized');
    % cm.Layout.Tile = 4;
    cm.Title       = 'Agreement: ica.final vs ICLabel (Preserving Other)';
    cm.FontName    = 'Helvetica';
    cm.FontSize    = 9;
elseif ~isempty(cohort_var_table)
    ax4 = nexttile(t, 4);
    var_cols = {'Brain', 'GenBad', 'Muscle', 'Eye', 'Heart', 'Channel'};
    var_mat  = cohort_var_table{:, var_cols};
    boxplot(ax4, var_mat, 'Labels', var_cols);
    grid(ax4, 'on');
    set(ax4, 'Box', 'off', 'FontName', 'Helvetica', 'FontSize', 9);
    ylabel(ax4, 'Variance Accounted For (%)', 'FontWeight', 'bold');
    title(ax4, 'Cohort Variance Distribution (ica.final.var)', 'FontSize', 11, 'FontWeight', 'bold');
else
    ax4 = nexttile(t, 4);
    text(ax4, 0.5, 0.5, 'No ICLabel predictions or variance summaries available', ...
        'HorizontalAlignment', 'center', 'Color', [0.45, 0.45, 0.45]);
    set(ax4, 'Box', 'off', 'XTick', [], 'YTick', []);
end

fig_file = fullfile(cfg.outdir, 'cohort_ica_benchmark.png');
exportgraphics(fh, fig_file, 'Resolution', 300);
fprintf('Diagnostic dashboard written to: %s\n', fig_file);
end

% =========================================================================
% Local Helpers
% =========================================================================
function tf = safe_bool(st, field_name, idx)
tf = false;
if isfield(st, field_name) && ~isempty(st.(field_name))
    val = st.(field_name);
    if length(val) >= idx
        tf = logical(val(idx));
    end
end
end

function val = safe_scalar(st, field_name)
if isfield(st, field_name) && ~isempty(st.(field_name)) && isnumeric(st.(field_name))
    val = double(st.(field_name)(1));
else
    val = NaN;
end
end

function val = safe_field(st, field_name)
if isfield(st, field_name)
    val = st.(field_name);
else
    val = [];
end
end

function vec = extract_vec(val, n_target)
if isempty(val)
    vec = nan(n_target, 1);
else
    vec = double(val(:));
    if length(vec) < n_target
        vec(end+1:n_target, 1) = NaN;
    elseif length(vec) > n_target
        vec = vec(1:n_target);
    end
end
end

function feats = harvest_numeric_vectors(st, prefix, n_target)
feats = struct();
if ~isstruct(st), return; end
f_names = fieldnames(st);

for i = 1:length(f_names)
    fname = f_names{i};
    val   = st.(fname);
    valid_key = matlab.lang.makeValidName([prefix fname]);

    if isstruct(val)
        sub_feats = harvest_numeric_vectors(val, [prefix fname '_'], n_target);
        sub_names = fieldnames(sub_feats);
        for s = 1:length(sub_names)
            feats.(sub_names{s}) = sub_feats.(sub_names{s});
        end
    elseif isnumeric(val) && numel(val) == n_target
        feats.(valid_key) = double(val(:));
    end
end
end