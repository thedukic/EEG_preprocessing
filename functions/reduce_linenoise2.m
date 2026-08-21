function DATA = reduce_linenoise2(DATA, cfg)

% Length of epochs
length_epoch_sec = 4; % seconds

fprintf('\n================================\n');
fprintf('Reducing 50 Hz noise leftovers\n');
fprintf('================================\n');

% Define params
threshold_p = 0.01; % Temporal stability threshold (p values from a paired t-test across trials)
threshold_z = 5.00; % Local SNR threshold (standard deviations above spectral background)

fprintf('Number of data channels: %d\n', DATA(1).nbchan);
fprintf('P-value threshold set to: %1.4f\n', threshold_p);
fprintf('Spectral Z-Score threshold set to: %1.1f\n', threshold_z);

% Custom masked colormap for Topoplots: Pure white for non-significant, Reds for significant
% masked_cmap1 = [1 1 1; brewermap(255, 'Reds')];
masked_cmap1 = brewermap(255, 'Reds');
masked_cmap2 = brewermap(255, '*RdBu');

chaneeg = strcmp({DATA(1).chanlocs.type}, 'EEG');

num_channels = DATA(1).nbchan;
num_block = length(DATA);
badelec = cell(num_block, 1);
pval    = cell(num_block, 1);

% Figure setup (4 columns)
fh = figure('Name', 'Line Noise QA', 'Color', 'w', 'Position', [100, 100, 1400, 300 * num_block], 'Visible', cfg.figure.visible);
t = tiledlayout(num_block, 4, "TileSpacing", "compact", "Padding", "compact");

for i_block = 1:num_block
    % Epoching
    length_epoch_sample = length_epoch_sec * DATA(i_block).srate;
    num_epoch = floor(size(DATA(i_block).data, 2) / length_epoch_sample);

    if num_epoch > 10
        dataeeg = reshape(DATA(i_block).data(:, 1:num_epoch*length_epoch_sample), DATA(i_block).nbchan, length_epoch_sample, num_epoch);

        % Estimate power
        winfunc    = hann(length_epoch_sample);
        psdspectra = NaN(length_epoch_sample/2+1, DATA(i_block).nbchan, num_epoch);
        for i_trial = 1:num_epoch
            [psdspectra(:, :, i_trial), freq] = pwelch(dataeeg(:, :, i_trial)', winfunc, 0, length_epoch_sample, DATA(i_block).srate);
        end

        % Average across trials (Linear space) and Log transform
        psd_avg_linear = mean(psdspectra, 3)';
        psdspectra2a   = 10 * log10(psd_avg_linear); % Corrected trial-averaged dB power (Channels x Freq)

        % Log-transform individual trials for the statistical variance test
        psdspectra_log = 10 * log10(psdspectra);

        % Define frequency bands
        lineFreq = 50;     % Target
        noiseBW  = 0;      % Half-width of the peak
        guardBW  = 1.5;    % Gap to skip to avoid spectral leakage around the line noise
        neighBW  = 5;      % Width of the background/neighbour bins

        if i_block == 1
            freq_50 = (freq >= (lineFreq - noiseBW)) & (freq <= (lineFreq + noiseBW));

            low_band  = (freq >= (lineFreq - noiseBW - neighBW)) & ...
                (freq <= (lineFreq - noiseBW - guardBW));
            high_band = (freq >= (lineFreq + noiseBW + guardBW)) & ...
                (freq <= (lineFreq + noiseBW + neighBW));

            freq_ref = low_band | high_band;

            % Formatted Frequency Range Logger
            fprintf('\nSpectral search bands established:\n');
            fprintf('  Target Line Noise Band (50 Hz): %1.2f Hz to %1.2f Hz\n', ...
                min(freq(freq_50)), max(freq(freq_50)));
            fprintf('  Background Reference Band (Low):  %1.2f Hz to %1.2f Hz\n', ...
                min(freq(low_band)), max(freq(low_band)));
            fprintf('  Background Reference Band (High): %1.2f Hz to %1.2f Hz\n\n', ...
                min(freq(high_band)), max(freq(high_band)));
        end

        % -----------------------------------------------------------------
        % Robust Detection Logic
        % -----------------------------------------------------------------

        % 1. Temporal Stability: Paired t-test across trials (using log power for normal distribution)
        X = squeeze(mean(psdspectra_log(freq_50, :, :), 1))';
        Y = squeeze(mean(psdspectra_log(freq_ref, :, :), 1))';
        [~, pval{i_block}] = ttest(X, Y, 'tail', 'right');

        % 2. Structural Anomaly: Linear Baseline Interpolation
        z_peak = NaN(num_channels, 1);
        poly_fits = NaN(num_channels, 2); % Store the fits for plotting later

        % Extract frequency vectors for the fit
        f_bg = freq(freq_ref);
        f_target = freq(freq_50);

        for i_chan = 1:num_channels
            % Extract the 1D trial-averaged log-power for this specific channel
            power_bg = psdspectra2a(i_chan, freq_ref)';
            power_target = psdspectra2a(i_chan, freq_50)';

            % Fit a 1st-degree polynomial (linear line) to the background noise
            p = polyfit(f_bg, power_bg, 1);
            poly_fits(i_chan, :) = p; % Save for Plot 4

            % Calculate what the power at 50 Hz should be based on the line
            expected_baseline = polyval(p, f_target);

            % Calculate the raw residual (how far the peak rises above the line)
            peak_residual = mean(power_target - expected_baseline);

            % Calculate the variance (noise) of the background fit
            fit_error = power_bg - polyval(p, f_bg);
            bg_std = std(fit_error);

            % Final Z-score: Signal-to-Noise Ratio of the peak
            z_peak(i_chan) = peak_residual / bg_std;
        end

        % 3. Final Flagging
        badelec_tmp1 = pval{i_block}' < threshold_p;
        badelec_tmp2 = z_peak > threshold_z;

        % Channels must be both statistically consistent and structurally anomalous
        badelec{i_block} = find(badelec_tmp1 & badelec_tmp2);

        % -----------------------------------------------------------------
        % Visualisation: Column 1 - P-Value Masked Topoplot
        % -----------------------------------------------------------------
        nexttile((i_block-1)*4 + 1);

        topo_data_p = -log10(pval{i_block}(chaneeg))';
        % topo_data_p(pval{i_block}(chaneeg)' >= threshold_p) = 0;

        topoplot(topo_data_p, DATA(i_block).chanlocs, 'maplimits', [0 -log10(1e-6)], ...
            'headrad', 'rim', 'whitebk', 'on', 'style', 'map', 'electrodes', 'on', ...
            'emarker2', {badelec{i_block}, 'd', 'k', 8, 1}, 'shading', 'flat');

        title(sprintf('Block %d Temporal P-Value', i_block), 'FontSize', 11, 'FontWeight', 'bold');
        colormap(gca, masked_cmap1);
        hcb = colorbar;
        hcb.Title.String = "-log_{10}(P)";

        % -----------------------------------------------------------------
        % Visualisation: Column 2 - Z-Score Masked Topoplot
        % -----------------------------------------------------------------
        nexttile((i_block-1)*4 + 2);

        topo_data_z = z_peak(chaneeg)';
        % topo_data_z(z_peak(chaneeg)' <= threshold_z) = 0;

        % max_z_lim = max([threshold_z * 2, max(z_peak(chaneeg))]);
        topoplot(topo_data_z, DATA(i_block).chanlocs, 'maplimits', 5 * [-1 1], ...
            'headrad', 'rim', 'whitebk', 'on', 'style', 'map', 'electrodes', 'on', ...
            'emarker2', {badelec{i_block}, 'd', 'k', 8, 1}, 'shading', 'flat');

        title(sprintf('Block %d Spectral Z-Score', i_block), 'FontSize', 11, 'FontWeight', 'bold');
        colormap(gca, masked_cmap2);
        hcb = colorbar;
        hcb.Title.String = "Z";

        % -----------------------------------------------------------------
        % Correction Processing
        % -----------------------------------------------------------------
        if ~isempty(badelec{i_block})
            num_channel_bad = length(badelec{i_block});
            fprintf('\nBlock %d: Leftover line noise detected in %d electrode(s):\n', i_block, num_channel_bad);
            fprintf('  %-8s | %-12s | %-12s | %-12s\n', 'Channel', 'P-value', 'Z-score', 'Slope');
            fprintf('  ---------------------------------------------------\n');

            for i_bad = 1:num_channel_bad
                ch_idx   = badelec{i_block}(i_bad);
                ch_label = DATA(i_block).chanlocs(ch_idx).labels;
                p_val    = pval{i_block}(ch_idx);
                z_val    = z_peak(ch_idx);
                slope    = poly_fits(ch_idx, 1);

                fprintf('  %-8s | %-12.4e | %-12.2f | %-12.2f\n', ...
                    ch_label, p_val, z_val, slope);
            end
            fprintf('  ---------------------------------------------------\n');
            fprintf('  Applying DFT filter correction...\n\n');

            % Spectrum interpolation
            DATA(i_block).data(badelec{i_block},:) = my_dftfilter(...
                DATA(i_block).data(badelec{i_block},:), ...
                DATA(i_block).srate, ...
                [50 100], 'neighbour', [1 1], [2 2]);

            % Recompute corrected spectra
            dataeeg_fixed = reshape(DATA(i_block).data(:, 1:num_epoch*length_epoch_sample), DATA(i_block).nbchan, length_epoch_sample, num_epoch);
            psdspectra_fixed = NaN(length_epoch_sample/2+1, DATA(i_block).nbchan, num_epoch);
            for i_trial = 1:num_epoch
                [psdspectra_fixed(:, :, i_trial), ~] = pwelch(dataeeg_fixed(:, :, i_trial)', winfunc, 0, length_epoch_sample, DATA(i_block).srate);
            end

            psdspectra2b = 10 * log10(mean(psdspectra_fixed, 3)');
            Pmean = mean(pval{i_block}(badelec{i_block}));
        else
            fprintf('Block %d: Nice! No leftover line noise is found.\n', i_block);
            num_channel_bad = 0;
            Pmean = NaN;
            psdspectra2b = psdspectra2a; % Fallback for plotting
        end

        % -----------------------------------------------------------------
        % Visualisation: Column 3 - Full Spectra Density Cloud
        % -----------------------------------------------------------------
        ax3 = nexttile((i_block-1)*4 + 3); hold(ax3, 'on');

        % Plot healthy channels as a faint grey background cloud
        good_elec = setdiff(1:DATA(i_block).nbchan, badelec{i_block});
        for i_good = 1:length(good_elec)
            h = plot(ax3, freq, psdspectra2a(good_elec(i_good), :), 'Color', [0.6 0.6 0.6]);
            h.Color(4) = 0.15; % 15% opacity
        end

        % Highlight bad channels in bold red
        for i_good = 1:length(badelec{i_block})
            plot(ax3, freq, psdspectra2a(badelec{i_block}(i_good), :), 'Color', [0.8 0.2 0.2], 'LineWidth', 1.2);
        end

        grid(ax3, 'on');
        set(ax3, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'Box', 'off');
        xlim(ax3, [0 128]); ylim(ax3, [-50 40]);
        xlabel(ax3, 'Frequency (Hz)'); ylabel(ax3, 'Power (dB)');
        title(ax3, 'Pre-Correction Spectra', 'FontSize', 11, 'FontWeight', 'bold');

        % -----------------------------------------------------------------
        % Visualisation: Column 4 - Before & After Zoom with Linear Fit
        % -----------------------------------------------------------------
        ax4 = nexttile((i_block-1)*4 + 4); hold(ax4, 'on');

        if ~isempty(badelec{i_block})
            % Center the view around 0 dB relative to the background
            for i_bad = 1:length(badelec{i_block})
                ch_idx = badelec{i_block}(i_bad);

                % Baseline shifts for plotting
                baseline_a = mean(psdspectra2a(ch_idx, freq_ref));
                baseline_b = mean(psdspectra2b(ch_idx, freq_ref));

                % 1. Pre-correction (Red)
                hA = plot(ax4, freq, psdspectra2a(ch_idx, :) - baseline_a, 'Color', [0.8 0.2 0.2]);
                hA.Color(4) = 0.4;

                % 2. Post-correction (Teal)
                hB = plot(ax4, freq, psdspectra2b(ch_idx, :) - baseline_b, 'Color', [0.2 0.6 0.5], 'LineWidth', 1.2);
                hB.Color(4) = 0.8;

                % 3. Estimated Linear Fit (Dashed Grey)
                fit_line = polyval(poly_fits(ch_idx, :), freq);
                hC = plot(ax4, freq, fit_line - baseline_a, 'Color', 'k', 'LineStyle', '--');
                hC.Color(4) = 0.5; % 50% opacity for the fit line
            end
        else
            title(ax4, 'No Correction Required', 'FontSize', 11, 'FontWeight', 'bold');
            axis(ax4, 'off');
        end

        if ~isempty(badelec{i_block})
            grid(ax4, 'on');
            set(ax4, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'Box', 'off');
            xlim(ax4, [45 55]); ylim(ax4, [-10 20]);
            xlabel(ax4, 'Frequency (Hz)'); ylabel(ax4, 'Relative Power (dB)');

            if Pmean < 0.001
                P_display = '< 0.001';
            else
                P_display = sprintf('= %.3f', Pmean);
            end
            title(ax4, sprintf('Correction Zoom (N = %d, P %s)', num_channel_bad, P_display), 'FontSize', 11, 'FontWeight', 'bold');
        end
    else
        warning('Block %d: Too little data (N = %ds) in this recording block.', i_block, num_epoch*length_epoch_sec);
    end
end

% Save
plotX = 36;
plotY = max(10, num_block * 7);
save_figure(fh, DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_linenoise_2'], [plotX plotY]);

% Fix for reporting
for i_block = 1:num_block
    if ~isempty(badelec{i_block})
        pval{i_block} = pval{i_block}(badelec{i_block});
        badelec{i_block} = {DATA(i_block).chanlocs(badelec{i_block}).labels};
    else
        pval{i_block} = [];
        badelec{i_block} = {};
    end
end

% Log / Report
for i_block = 1:num_block
    DATA(i_block).ALSUTRECHT.LineNoiseCleaning2.badelec = badelec;
    DATA(i_block).ALSUTRECHT.LineNoiseCleaning2.pval    = pval;
end

fprintf(DATA(1).ALSUTRECHT.subject.fid, '\n---------------------------------------------------------\n');
fprintf(DATA(1).ALSUTRECHT.subject.fid, 'Electrodes with leftover 50 Hz noise fixed\n');
fprintf(DATA(1).ALSUTRECHT.subject.fid, '---------------------------------------------------------\n');

all_badelec = [badelec{:}];
all_pval = [pval{:}];

if isempty(all_badelec)
    fprintf(DATA(1).ALSUTRECHT.subject.fid, 'Electrodes: None\n');
    fprintf(DATA(1).ALSUTRECHT.subject.fid, 'P-values:   None\n');
else
    fprintf(DATA(1).ALSUTRECHT.subject.fid, 'Electrodes: %s\n', strjoin(all_badelec, ', '));
    str_pval = arrayfun(@(x) num2str(x, '%1.3f'), all_pval, 'UniformOutput', false);
    fprintf(DATA(1).ALSUTRECHT.subject.fid, 'P-values:   %s\n', strjoin(str_pval, ', '));
end

end