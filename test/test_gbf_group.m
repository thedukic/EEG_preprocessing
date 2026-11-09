function test_gbf_group

% =========================================================================
% Group-Level Diagnostic Script: Relative Band Power
% Loads pre-calculated PSD data across subjects (both pre- and post-GBF),
% evaluates relative band power, compares 41-44 Hz peak origins,
% and visualises group spectra and topographies.
% =========================================================================

base_dir = 'F:\3_PREPROCESSED_DATA\RS\ALS\T1';

% Canonical frequency bands and normalisation range
bands.Theta  = [4 8];
bands.Alpha  = [8 13];
bands.Beta   = [13 30];
bands.Gamma  = [41 44];
bands.Gamma1 = [30 48];
bands.Gamma2 = [52 100];
total_range  = [1 100];

band_names = fieldnames(bands);
n_bands    = numel(band_names);

% List subject directories starting with 'ALS'
dir_contents = dir(fullfile(base_dir, 'ALS*'));
dir_flags    = [dir_contents.isdir];
sub_folders  = dir_contents(dir_flags);

num_subject = numel(sub_folders);
fprintf('Found %d prospective subject folders.\n', num_subject);

% Pre-allocate storage
group_rel_power_band_post = NaN(num_subject, 128, n_bands);
group_rel_power_band_pre  = NaN(num_subject, 128, n_bands);
group_rel_power_all_post  = [];
group_rel_power_all_pre   = [];
group_raw_psd_post        = [];
group_raw_psd_pre         = [];

valid_sub_count = 0;
failed_subjects = {};

for s = 1:num_subject
    sub_id = sub_folders(s).name;
    file_path = fullfile(base_dir, sub_id, '2', 'data', [sub_id '_T1_RS_rawpower.mat']);

    if ~isfile(file_path)
        warning('File not found: %s', file_path);
        failed_subjects{end+1} = sub_id;
        continue;
    end

    fprintf('[%3d] Processing: %s\n', valid_sub_count + 1, sub_id);

    load(file_path, 'data_power');

    freqs        = data_power.freq_post;
    psd_post_sub = data_power.psd_post(1:128, :);
    psd_pre_sub  = data_power.psd_pre(1:128, :);

    if valid_sub_count == 0
        n_freqs = length(freqs);
        group_rel_power_all_post = NaN(num_subject, 128, n_freqs);
        group_rel_power_all_pre  = NaN(num_subject, 128, n_freqs);
        group_raw_psd_post       = NaN(num_subject, 128, n_freqs);
        group_raw_psd_pre        = NaN(num_subject, 128, n_freqs);
    end

    total_idx      = (freqs >= total_range(1)) & (freqs <= total_range(2));
    total_pwr_post = sum(psd_post_sub(:, total_idx), 2);
    total_pwr_pre  = sum(psd_pre_sub(:, total_idx), 2);

    sub_band_pwr_post = zeros(128, n_bands);
    sub_band_pwr_pre  = zeros(128, n_bands);
    for b = 1:n_bands
        b_lims = bands.(band_names{b});
        b_idx  = (freqs >= b_lims(1)) & (freqs <= b_lims(2));
        sub_band_pwr_post(:, b) = (sum(psd_post_sub(:, b_idx), 2) ./ total_pwr_post) * 100;
        sub_band_pwr_pre(:, b)  = (sum(psd_pre_sub(:, b_idx), 2)  ./ total_pwr_pre)  * 100;
    end

    valid_sub_count = valid_sub_count + 1;
    group_rel_power_band_post(valid_sub_count, :, :) = sub_band_pwr_post;
    group_rel_power_band_pre(valid_sub_count, :, :)  = sub_band_pwr_pre;
    
    group_rel_power_all_post(valid_sub_count, :, :)  = psd_post_sub ./ total_pwr_post;
    group_rel_power_all_pre(valid_sub_count, :, :)   = psd_pre_sub  ./ total_pwr_pre;

    group_raw_psd_post(valid_sub_count, :, :)        = psd_post_sub;
    group_raw_psd_pre(valid_sub_count, :, :)         = psd_pre_sub;
end

group_rel_power_band_post = group_rel_power_band_post(1:valid_sub_count, :, :);
group_rel_power_band_pre  = group_rel_power_band_pre(1:valid_sub_count, :, :);
group_rel_power_all_post  = group_rel_power_all_post(1:valid_sub_count, :, :);
group_rel_power_all_pre   = group_rel_power_all_pre(1:valid_sub_count, :, :);
group_raw_psd_post        = group_raw_psd_post(1:valid_sub_count, :, :);
group_raw_psd_pre         = group_raw_psd_pre(1:valid_sub_count, :, :);

fprintf('\nSuccessfully loaded and processed N = %d subjects.\n', valid_sub_count);
if ~isempty(failed_subjects)
    fprintf('Skipped %d folders due to missing files.\n', numel(failed_subjects));
end

% =========================================================================
% Figure 1: Diagnostic Check: 41-44 Hz Peak Origin (Pre vs Post GBF)
% =========================================================================
mean_psd_pre_ch   = squeeze(mean(group_raw_psd_pre, 1));
mean_psd_post_ch  = squeeze(mean(group_raw_psd_post, 1));
grand_raw_pre     = mean(mean_psd_pre_ch, 1);
grand_raw_post    = mean(mean_psd_post_ch, 1);

plot_zoom_idx = (freqs >= 30) & (freqs <= 60);
f_zoom        = freqs(plot_zoom_idx);

idx_42 = (freqs >= 41) & (freqs <= 44);
pwr_42_pre  = sum(group_raw_psd_pre(:, :, idx_42), 3);
pwr_42_post = sum(group_raw_psd_post(:, :, idx_42), 3);

cv_42_pre  = std(pwr_42_pre, 0, 1)  ./ max(mean(pwr_42_pre, 1), eps);
cv_42_post = std(pwr_42_post, 0, 1) ./ max(mean(pwr_42_post, 1), eps);

cmap_cv = brewermap([], 'YlOrRd');
figure('Name', sprintf('Diagnostic: 42 Hz Bump Origin (N = %d)', valid_sub_count), ...
       'Color', 'w', 'Position', [100, 150, 1100, 380]);

subplot(1, 3, 1);
plot(f_zoom, 10*log10(grand_raw_pre(plot_zoom_idx)), 'Color', [0.85 0.20 0.20], 'LineWidth', 2); hold on;
plot(f_zoom, 10*log10(grand_raw_post(plot_zoom_idx)), 'Color', [0.15 0.45 0.85], 'LineWidth', 2);
xline(40, 'k--', 'GBF Taper End (40 Hz)', 'LabelVerticalAlignment', 'bottom');
grid on; box on;
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Power (10\cdotlog_{10} \muV^2/Hz)', 'FontSize', 11);
title('Grand Mean PSD (30 to 60 Hz)', 'FontSize', 12);
legend({'Raw (Pre-GBF)', 'Cleaned (Post-GBF)'}, 'Location', 'northeast');

ax_cv_pre = subplot(1, 3, 2);
c_max_cv = max([cv_42_pre, cv_42_post]) * 1.05;
mytopoplot(cv_42_pre, [], [], ax_cv_pre, [0, c_max_cv]);
title('Pre-GBF CV (41 to 44 Hz)', 'FontSize', 11);
colormap(ax_cv_pre, cmap_cv);
c_pre = colorbar;
ylabel(c_pre, 'CV', 'FontSize', 9);

ax_cv_post = subplot(1, 3, 3);
mytopoplot(cv_42_post, [], [], ax_cv_post, [0, c_max_cv]);
title('Post-GBF CV (41 to 44 Hz)', 'FontSize', 11);
colormap(ax_cv_post, cmap_cv);
c_post = colorbar;
ylabel(c_post, 'CV', 'FontSize', 9);

% =========================================================================
% Figure 2: Group-Level Average Scalp Power Spectrum Across All Frequencies
% =========================================================================
mean_sub_psd_post = squeeze(mean(group_rel_power_all_post, 1));
grand_mean_post   = mean(mean_sub_psd_post, 1);

sd_channels_post = std(mean_sub_psd_post, 0, 1);

sub_chan_mean_post = squeeze(mean(group_rel_power_all_post, 2));
sem_subjects_post  = std(sub_chan_mean_post, 0, 1) ./ sqrt(valid_sub_count);

plot_f_idx = (freqs >= 1) & (freqs <= 85);
f_axis     = freqs(plot_f_idx);

spec_mean    = grand_mean_post(plot_f_idx);
spec_chan_sd = sd_channels_post(plot_f_idx);
spec_sub_sem = sem_subjects_post(plot_f_idx);

db_mean     = 10 * log10(spec_mean);
db_chan_up  = 10 * log10(spec_mean + spec_chan_sd);
db_chan_low = 10 * log10(max(spec_mean - spec_chan_sd, eps));

db_sem_up   = 10 * log10(spec_mean + spec_sub_sem);
db_sem_low  = 10 * log10(max(spec_mean - spec_sub_sem, eps));

figure('Name', sprintf('Group Grand Average Scalp Spectrum (N = %d)', valid_sub_count), ...
       'Color', 'w', 'Position', [200, 200, 1050, 430]);

subplot(1, 2, 1); hold on;
h_chan = fill([f_axis', fliplr(f_axis')], [db_chan_up, fliplr(db_chan_low)], ...
              [0.70 0.80 0.95], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
h_sem  = fill([f_axis', fliplr(f_axis')], [db_sem_up, fliplr(db_sem_low)], ...
              [0.15 0.45 0.85], 'FaceAlpha', 0.50, 'EdgeColor', 'none');
h_line = plot(f_axis, db_mean, 'Color', [0.10 0.30 0.70], 'LineWidth', 2);
grid on; box on;
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Relative Power (10\cdotlog_{10} fraction/Hz)', 'FontSize', 11);
title(sprintf('Linear Spectrum (N = %d)', valid_sub_count), 'FontSize', 12);
xlim([1 85]);
legend([h_line, h_sem, h_chan], ...
       {'Grand Mean', 'Group SEM (Subjects)', 'Spatial SD (Channels)'}, ...
       'Location', 'northeast', 'FontSize', 9);

subplot(1, 2, 2); hold on;
log_f = log10(f_axis);
fill([log_f', fliplr(log_f')], [db_chan_up, fliplr(db_chan_low)], ...
     [0.70 0.80 0.95], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
fill([log_f', fliplr(log_f')], [db_sem_up, fliplr(db_sem_low)], ...
     [0.15 0.45 0.85], 'FaceAlpha', 0.50, 'EdgeColor', 'none');
plot(log_f, db_mean, 'Color', [0.10 0.30 0.70], 'LineWidth', 2);
grid on; box on;
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Relative Power (10\cdotlog_{10} fraction/Hz)', 'FontSize', 11);
title('Aperiodic Trend (Log-Log Scale)', 'FontSize', 12);
xtick_freqs = [2, 4, 8, 13, 30, 50, 80];
set(gca, 'XTick', log10(xtick_freqs), 'XTickLabel', string(xtick_freqs));
xlim([log10(1) log10(85)]);

% =========================================================================
% Figure 3: Group Topographies: Mean vs Inter-Subject CV (Post-GBF)
% =========================================================================
mean_group_pwr_post = squeeze(mean(group_rel_power_band_post, 1));
std_group_pwr_post  = squeeze(std(group_rel_power_band_post, 0, 1));
cv_group_pwr_post   = std_group_pwr_post ./ max(mean_group_pwr_post, eps);

cmap_pos = brewermap([], 'YlOrRd');
figure('Name', sprintf('Group Topographies: Mean vs Subject CV (N = %d)', valid_sub_count), ...
       'Color', 'w', 'Position', [150, 100, 320 * n_bands, 580]);

for b = 1:n_bands
    b_name = band_names{b};
    b_lims = bands.(b_name);

    mean_topo = mean_group_pwr_post(:, b);
    cv_topo   = cv_group_pwr_post(:, b);

    % Row 1: Group Mean Relative Power
    ax1 = subplot(2, n_bands, b);
    c_lims_mean = [0, max(mean_topo) * 1.05];
    mytopoplot(mean_topo, [], [], ax1, c_lims_mean);
    title(sprintf('%s [%d to %d Hz]\nGroup Mean (%%)', b_name, b_lims(1), b_lims(2)), 'FontSize', 11);
    colormap(ax1, cmap_pos);
    c1 = colorbar;
    ylabel(c1, '% Total Power', 'FontSize', 9);

    % Row 2: Inter-Subject CV
    ax2 = subplot(2, n_bands, n_bands + b);
    c_lims_cv = [0, max(cv_topo) * 1.05];
    mytopoplot(cv_topo, [], [], ax2, c_lims_cv);
    title(sprintf('%s [%d to %d Hz]\nInter-Subject CV', b_name, b_lims(1), b_lims(2)), 'FontSize', 11);
    colormap(ax2, cmap_pos);
    c2 = colorbar;
    ylabel(c2, 'CV (SD / Mean)', 'FontSize', 9);
end

fprintf('Diagnostics rendering complete.\n');

end