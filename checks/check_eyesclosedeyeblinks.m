function EEG = check_eyesclosedeyeblinks(EEG, EXT, cfg)

fprintf('\n================================\n');
fprintf('Detecting eye blinks in the eyes closed resting-state recording\n');
fprintf('================================\n');

% Select these
chaneeg  = strcmp({EEG.chanlocs.type},'EEG'); % & contains({EEG.chanlocs.labels},'C');
chaneog  = strcmp({EXT.chanlocs.labels},'VEOG'); assert(sum(chaneog) == 1);
chanlocs = EEG.chanlocs(chaneeg);

dataeeg  = EEG.data(chaneeg,:);
dataeog  = EXT.data(chaneog,:);

% Only EC
ec_mask = EEG.ALSUTRECHT.blockinfo.rs_mask(~EEG.ALSUTRECHT.blockinfo.eo_mask,:);
ec_mask = any(ec_mask,1);

% Check if the participant has any EC data
if any(ec_mask)
    % EC data exists
    assert(length(ec_mask) == size(dataeeg,2));

    dataeeg   = dataeeg(:, ec_mask);
    dataeogeo = dataeog(~ec_mask)';
    dataeogec = dataeog(ec_mask)';
    times     = EEG.times(ec_mask) / 1000;
    times     = times - times(1);

    %% Strong EMG
    % Temporarily filter EOG specifically for blink morphology
    [bl, al] = butter(2, 5/(EEG.srate/2), 'low'); % Keep at 5 Hz to smooth flutters
    [bh, ah] = butter(2, 0.1/(EEG.srate/2), 'high'); % Lowered from 1 Hz to 0.1 Hz

    assert(isstable(bl,al), 'Lowpass filter unstable.');
    assert(isstable(bh,ah), 'Highpass filter unstable.');

    dataeogeo = filtfilt(bl,al,dataeogeo);
    dataeogec = filtfilt(bl,al,dataeogec);

    dataeogeo = filtfilt(bh,ah,dataeogeo);
    dataeogec = filtfilt(bh,ah,dataeogec);

    % Baseline correct
    dataeogeo = dataeogeo - trimmean(dataeogeo,10);
    dataeogec = dataeogec - trimmean(dataeogec,10);

    % =========================================================================
    % Standardise copies of the data to assess shape without variance penalty
    % We use new variables (_z) so the raw uV data is preserved for blink detection
    dataeogeo_z = (dataeogeo - mean(dataeogeo)) / std(dataeogeo);
    dataeogec_z = (dataeogec - mean(dataeogec)) / std(dataeogec);

    % 1. Parameter Comparison (The Shape)
    pdEO = fitdist(dataeogeo_z, 'tLocationScale');
    pdEC = fitdist(dataeogec_z, 'tLocationScale');

    % 2. Kolmogorov-Smirnov Distance
    [~, ~, ks2stat] = kstest2(dataeogeo_z, dataeogec_z);

    % 3. Kullback-Leibler Divergence
    myXlim_z = [min([dataeogeo_z; dataeogec_z]), max([dataeogeo_z; dataeogec_z])];
    x_vals = linspace(myXlim_z(1), myXlim_z(2), 1000)';

    pdfEO = pdf(pdEO, x_vals);
    pdfEC = pdf(pdEC, x_vals);

    % Normalise to ensure they sum to 1 over the discrete grid
    pdfEO = pdfEO / sum(pdfEO);
    pdfEC = pdfEC / sum(pdfEC);

    % Prevent log(0) errors
    pdfEO = max(pdfEO, eps);
    pdfEC = max(pdfEC, eps);

    kl_div = sum(pdfEC .* log(pdfEC ./ pdfEO));

    % Combine
    ec_metrics = [pdEC.nu, ks2stat, kl_div];

    % =========================================================================
    % Console Report: Statistical Shape Comparison
    fprintf('Statistical Shape Comparison: Eyes Open vs. Eyes Closed\n');
    fprintf('1. Distribution Parameters (Standardised Data):\n');
    fprintf('   EO: mu = %5.3f, sigma = %5.3f, nu = %5.3f\n', pdEO.mu, pdEO.sigma, pdEO.nu);
    fprintf('   EC: mu = %5.3f, sigma = %5.3f, nu = %5.3f\n', pdEC.mu, pdEC.sigma, pdEC.nu);

    if pdEC.nu > pdEO.nu
        fprintf('   -> Good: EC data has a higher nu parameter, confirming fewer extreme outliers.\n');
        flagGoodFit = true;
    else
        fprintf('   -> Warning: EC data has heavier tails than EO. Inspect for large artefacts.\n');
        flagGoodFit = false;
    end

    fprintf('\n2. Kolmogorov-Smirnov (KS) Test:\n');
    fprintf('   Maximum Distance: %.4f\n', ks2stat);
    fprintf('   -> The cumulative shapes differ by a maximum of %1.1f%%.\n', ks2stat * 100);

    fprintf('\n3. Kullback-Leibler (KL) Divergence:\n');
    fprintf('   Information Loss: %.4f\n', kl_div);
    if kl_div > 0.5
        fprintf('   -> High divergence. The core structure of the baseline noise has shifted significantly.\n');
    elseif kl_div > 0.1
        fprintf('   -> Moderate divergence. There is a noticeable structural shift in the noise.\n');
    else
        fprintf('   -> Low divergence. The baseline profiles remain structurally similar.\n');
    end

    % =========================================================================
    % Plot histograms and fitted PDFs
    fh = figure('Visible', cfg.figure.visible);
    % Increased rows to 3 to accommodate the continuous plot at the bottom
    th = tiledlayout(3,2);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    % --- NEW: Automated Conclusion Summary ---
    if flagGoodFit
        if kl_div > 0.1
            status_text = 'Conclusion: Good compliance. EC has fewer blinks and distinct noise structure from EO.';
            status_color = [0.1 0.5 0.1]; % Dark green
        else
            status_text = 'Conclusion: Acceptable. EC has fewer blinks, but baseline structure is highly similar to EO.';
            status_color = [0.4 0.4 0.4]; % Grey
        end
    else
        status_text = 'Conclusion: WARNING! EC has heavier tails than EO. Participant may not have kept eyes closed.';
        status_color = [0.8 0.2 0.2]; % Red
    end

    % Add the summary as a super-title to the whole figure
    sgtitle({sprintf('Eyes Closed vs Eyes Open Comparison (KS: %.1f%%, KL: %.2f)', ks2stat*100, kl_div), status_text}, ...
        'Color', status_color, 'FontWeight', 'bold', 'FontSize', 12);
    % -----------------------------------------

    nexttile; hold on;
    hb = histogram(dataeogeo_z, 'Normalization', 'pdf');
    hb.FaceColor = 0.7*ones(1,3);
    fplot(@(x) pdf(pdEO, x), myXlim_z, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
    title(['EO VEOG (Z-scored), \nu = ' num2str(round(pdEO.nu, 2))]);
    xlabel('Standardised Amplitude (Z)'); ylabel('Probability Density'); xlim(myXlim_z);

    nexttile; hold on;
    hb = histogram(dataeogec_z, 'Normalization', 'pdf');
    hb.FaceColor = 0.7*ones(1,3);
    fplot(@(x) pdf(pdEC, x), myXlim_z, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
    title(['EC VEOG (Z-scored), \nu = ' num2str(round(pdEC.nu, 2))]);
    xlabel('Standardised Amplitude (Z)'); ylabel('Probability Density'); xlim(myXlim_z);

    % =========================================================================
    % Estimate the threshold using the EC RS raw data itself
    % We use the IQR of the EC data because the baseline noise is different
    threshold = prctile(dataeogec, 75) + 5 * iqr(dataeogec);

    % Enforce a strict, higher absolute minimum to avoid catching slow rolling eyes
    threshold = max(threshold, 100);
    fprintf('Eye blink detection using a threshold of %1.0f uV.\n', threshold);

    % Find the blinks in EC using temporal constraints
    dataDetect = dataeogec;
    min_distance_sec = 0.5; % Minimum half-second between blinks
    [qrspeaks, locs] = findpeaks(dataDetect, times, ...
        'MinPeakHeight', threshold, ...
        'MinPeakDistance', min_distance_sec);

    % =========================================================================
    NEOG = length(locs);
    if NEOG > 0
        badEpoch2 = round(locs * EEG.srate);
        badEpoch2 = [badEpoch2-80; badEpoch2+80]';
        N = length(badEpoch2(1,1):badEpoch2(1,2));

        badEpoch2(badEpoch2 < 1) = 1;
        badEpoch2(badEpoch2 > length(dataDetect)) = length(dataDetect);

        EOG = NaN(NEOG,N);
        for i = 1:NEOG
            if length(badEpoch2(i,1):badEpoch2(i,2)) == N
                EOG(i,:) = dataDetect(badEpoch2(i,1):badEpoch2(i,2));
            end
        end
        EOG(isnan(EOG(:,1)),:) = [];
        mEOG = mean(EOG,1);

        nexttile; hold on;
        F = (0:size(EOG,2)-1)./EEG.srate;

        % --- NEW: Darker Colormap Logic ---
        % Generate a larger colormap and chop off the lightest 40% of colours
        N_colors = size(EOG,1);
        full_cmap = brewermap(round(N_colors * 1.6), 'Greys');
        dark_cmap = full_cmap(end - N_colors + 1 : end, :);

        % IMPORTANT: Set ColorOrder BEFORE calling plot()
        set(gca, 'ColorOrder', dark_cmap);

        % Plot individual blinks (they will now use the dark_cmap)
        plot(F, EOG, 'LineWidth', 1.2);

        % Plot the median EOG on top in red
        plot(F, mEOG, 'Color', [0.8 0.1 0.1], 'LineWidth', 3);

        title(['N = ' num2str(NEOG) ', threshold = ' num2str(round(threshold)) ' \muV']);
        xlabel('Time (s)'); ylabel('EC VEOG amplitude');

        EEGtmp = NaN(NEOG,sum(chaneeg));
        for i = 1:NEOG
            if length(badEpoch2(i,1):badEpoch2(i,2)) == N
                EEGtmp(i,:) = mean(dataeeg(:,badEpoch2(i,1):badEpoch2(i,2)).^2,2);
            end
        end
        EEGtmp(isnan(EEGtmp(:,1)),:) = [];
        mEEGtmp = median(EEGtmp,1);

        myCmap = brewermap(128, 'Greys');
        % myCmap = brewermap(128, 'PuOr');
        nexttile;
        topoplot(mEEGtmp,chanlocs,'maplimits',max(abs(mEEGtmp))*[0 1],'headrad',0.5,'colormap',myCmap,'whitebk','on','electrodes','off','style','map','shading','interp');
        title('EC EEG timelocked to eye blinks');
        hcb = colorbar;
        hcb.Title.String = "uV^2";

        % Detected peaks
        nexttile(5, [1 2]); hold on;
        plot(times, dataDetect, 'Color', [0.4 0.4 0.4]);
        plot(locs, qrspeaks, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

        % Calculate block boundaries
        ec_blocks_only = EEG.ALSUTRECHT.blockinfo.rs_mask(~EEG.ALSUTRECHT.blockinfo.eo_mask, :);
        block_lengths = sum(ec_blocks_only, 2);
        boundary_indices = cumsum(block_lengths(1:end-1)); % Get indices for all boundaries except the very end

        % Draw vertical lines if boundaries exist
        if ~isempty(boundary_indices) && boundary_indices(1) > 0
            boundary_times = times(boundary_indices);
            xline(boundary_times, 'Color', [0.2 0.4 0.8], 'LineStyle', '--', 'LineWidth', 1.5);
        end

        title('Continuous EC VEOG with Detected Blinks');
        xlabel('Time (s)'); ylabel('Amplitude (uV)');
        xlim([times(1) times(end)]);

        % Report
        warning('The participant has eye blinks detected during EC.');
    else
        fprintf('Great! No eye blinks found during eyes closed resting-state recording.\n');
        mEEGtmp = NaN(1,128);
    end

    % Save
    % plotX=25; plotY=35;
    % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    % print(fh, fullfile(EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_detected_ecblinks']), '-dtiff', '-r200'); close(fh);
    save_figure(fh, EEG.ALSUTRECHT.subject.figures, [EEG.ALSUTRECHT.subject.id '_detected_ecblinks'], [25 35]);

else
    warning('Unexpected. The participant is missing EC data.');
    NEOG        = NaN;
    flagGoodFit = NaN;
    ec_metrics  = NaN;
    mEEGtmp     = NaN;
end

% Log
EEG.ALSUTRECHT.blockinfo.ec_NBlinks     = NEOG;
EEG.ALSUTRECHT.blockinfo.ec_flagGoodFit = flagGoodFit;
EEG.ALSUTRECHT.blockinfo.ec_shapeStats  = ec_metrics;
EEG.ALSUTRECHT.blockinfo.ec_blinksTopo  = mEEGtmp;

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Eyes closed resting-state eye blinks detection\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected: %d\n', NEOG);


end