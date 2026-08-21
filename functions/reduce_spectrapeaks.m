function DATA = reduce_spectrapeaks(DATA, cfg)
% Needs imporvements:
% 1. The minimum peak prominance
% 2. The check how well DSS is done, like ratios DSS1/DSS2 > X

% Use DSS to isolate the peaks
fprintf('\n================================\n');
fprintf('Removing additional peaks from the spectrum\n');
fprintf('================================\n');

% First components are most dominated by these peaks
NBLK = length(DATA);

% Power spectra
[psdspectra, freq, chaneeg] = estimate_power(DATA, 'speaks');

% % Special cases
% if strcmpi(EEG(1).ALSUTRECHT.subject.id,'ALS34280') && strcmpi(EEG(1).ALSUTRECHT.subject.visit,'T1') && ~strcmpi(EEG(1).ALSUTRECHT.subject.task,'SART')
%     % 1.65-1.7 Hz harmonics?
%     locs = [1.65 3.35 5.05 6.7 8.4 10.1 11.75 13.45 15.15];
%     KREMOVE = 1;
%
% elseif strcmpi(EEG(1).ALSUTRECHT.subject.id,'ALS34280') && strcmpi(EEG(1).ALSUTRECHT.subject.visit,'T1') && strcmpi(EEG(1).ALSUTRECHT.subject.task,'SART')
%     % 0.65-0.70 Hz and 1.65-1.7 Hz harmonics?
%     locs = [[0.65 1.35 2.00 2.65], [5.3 7.1 8.85 10.6 12.35 14.15]]';
%     KREMOVE = [1 5];
%
% elseif strcmpi(EEG(1).ALSUTRECHT.subject.id,'ALS36104') && strcmpi(EEG(1).ALSUTRECHT.subject.visit,'T1')
%     % 1.35 Hz harmonics?
%     % locs = [1.35 2.7 4.1 5.4 6.75 8 9.5 10.75 12.2 13.5 14.9]';
%     % locs = [1.35 2.7 4 4.15 5.45 6.75 6.9 8.15 12 12.5 13.5 14.9 16.3 17.6 19]';
%     % locs = [1.35 2.7 6.7]';
%     % locs = round((1:5)' * 1.35,2);
%     locs = [1.35 2.7 4.1 5.4 6.75];
%     KREMOVE = 1:4;
%
% elseif strcmpi(EEG(1).ALSUTRECHT.subject.id,'P111') && strcmpi(EEG(1).ALSUTRECHT.subject.visit,'T1')
%     % 1.5 Hz harmonics?
%     locs = [4.55 6.00 7.7];
%     KREMOVE = 1;
%
% elseif strcmpi(EEG(1).ALSUTRECHT.subject.id,'P117') && strcmpi(EEG(1).ALSUTRECHT.subject.visit,'T1')
%     % 1.5 Hz harmonics?
%     locs = [3.05 4.55 6.05 7.60];
%     KREMOVE = 1;
%
% else
%     % Identify the peaks automatically - does not work well
%     % locs = find_peaks(psdspectra,freq);
%     locs = [];
%     KREMOVE = NaN;
% end

% Find peaks
[locs, fh] = find_peaks(psdspectra, freq, cfg.figure.visible);

% plotX=30; plotY=30;
% set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
% set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
% print(fh, fullfile(DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_peaks_2']), '-dtiff', '-r200'); close(fh);
save_figure(fh, DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_peaks_1'], [30 30]);

% Rounding errors
locs = round(locs(:), 2);

% If any peaks found
if ~isempty(locs)
    fprintf('Peaks found at:\n');
    fprintf('f = %1.2f Hz\n',locs);

    % Covariance matrices of full band (c0) and filtered at the peaks (c1)
    NPTS = 2 * (length(freq)-1);
    data = cat(2, DATA(:).data)';
    [c0, c1] = nt_bias_fft(data, locs'/DATA(1).srate, NPTS);

    % DSS matrix
    [todss, pwr0, pwr1] = nt_dss0(c0, c1, [], []);
    p1 = pwr1 ./ pwr0;

    % figure;
    % for i = 1:NREMOVE
    %     mytopoplot(todss(:, i), [], '', nexttile); colorbar;
    % end

    % --- AUTOMATIC COMPONENT SELECTION ---
    % Calculate the baseline noise level using the tail of the DSS components
    % (assuming components 10 to end are mostly noise/baseline)
    if length(p1) > 10
        baseline_mean = mean(p1(10:end));
        baseline_std  = std(p1(10:end));
        score_threshold = baseline_mean + 5 * baseline_std;
    else
        % Fallback if very few components
        score_threshold = mean(p1) + std(p1);
    end

    KREMOVE = find(p1 > score_threshold);

    % Handle case where no components exceed the threshold
    if isempty(KREMOVE)
        fprintf('No knee detected. Better to avoid removing any data. Skipping...\n');
        % Log empty values before returning to keep EEG structure consistent
        for i = 1:NBLK
            DATA(i).ALSUTRECHT.otherPeakCleaning.peaks = locs;
            DATA(i).ALSUTRECHT.otherPeakCleaning.N     = 0;
            DATA(i).ALSUTRECHT.otherPeakCleaning.p1    = p1;
        end
        return;
    end

    % Safety check: do not remove more than 5 components
    if length(KREMOVE) > 5
        KREMOVE = KREMOVE(1:5);
    end

    % Check DSS1/DSS2 ratio safely
    if length(p1) >= 2
        dss_ratio = p1(1) / p1(2);
        if dss_ratio < 1.5
            fprintf('Warning: DSS separation ratio is low (%1.2f). Peak might be weak.\n', dss_ratio);
        end
    end

    % DSS components
    z = nt_mmat(data, todss);

    % Regress them out
    clean = nt_tsr(data, z(:, KREMOVE));

    % Split back the blocks
    EEGTMP = make_blockmasks(DATA);
    assert(size(clean,1) == size(EEGTMP(1).ALSUTRECHT.blockinfo.rs_mask,2));

    for i = 1:NBLK
        DATA(i).data = clean(EEGTMP(i).ALSUTRECHT.blockinfo.rs_mask(i,:),:)';
        assert(size(DATA(i).data,2) == DATA(i).pnts);
    end
    DATA = eeg_checkset(DATA);

    % Log / Report
    for i = 1:NBLK
        DATA(i).ALSUTRECHT.otherPeakCleaning.peaks = locs;
        DATA(i).ALSUTRECHT.otherPeakCleaning.N     = KREMOVE;
        DATA(i).ALSUTRECHT.otherPeakCleaning.p1    = p1;
    end
    fprintf('Done! Number of components removed: %d\n', length(KREMOVE));

    % =====================================================================
    psdspectra = log10(mean(psdspectra, 2));

    % tmp = nt_spect_plot2(nt_normcol(z(:,1:25)),NPTS,0,NPTS,EEG(1).srate);
    % for i = 1:5
    %      locs1{i} = find_peaks(tmp(:,i),freq);
    % end

    % Plot
    fh = figure('Visible', cfg.figure.visible);
    th = tiledlayout(2,3);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    freq = round(freq, 2);
    % plot the peaks
    nexttile; hold on;
    plot(freq,psdspectra);
    % plot([freqMin freqMin],[min(psdspectra) max(psdspectra)],'Color',0.5*ones(1,3)); axis tight
    scatter(locs, psdspectra(any(freq==locs',2)),'filled','MarkerFaceColor',0*ones(1,3)); axis tight;
    xlabel('Frequency (Hz)'); ylabel('log_{10}(Power)'); title('Spectra peak detection'); pbaspect([1.618 1 1]);

    % plot bias score
    nexttile;
    plot(p1,'.-'); xlabel('Component'); ylabel('Score'); title('DSS scores');
    axis tight; xlim([0 50]); pbaspect([1.618 1 1]);

    % plot spectra of DSS components
    nexttile;
    nt_spect_plot2(nt_normcol(z(:,1:25)),NPTS,0,NPTS,DATA(1).srate);
    title('Spectra of first 25 DSS components'); ylabel('Component'); xlabel('Frequency (Hz)');
    pbaspect([1.618 1 1]); colormap(brewermap(sum(chaneeg),'BrBG'));

    % plot the first DSS weights
    fromdss = pinv(todss);
    % mytopoplot(todss(:,1),[],['DSS1 score: ' num2str(round(p1(1)))],nexttile); % colorbar;
    mytopoplot(fromdss(1, 1:128),[],['DSS1 score: ' num2str(round(p1(1)))],nexttile); % colorbar;

    % plot spectra of data before and after removal of the peak component(s)
    nexttile; hold on;
    nt_spect_plot(data,NPTS,0,NPTS,DATA(1).srate);
    nt_spect_plot(clean,NPTS,0,NPTS,DATA(1).srate);
    nt_linecolors([],[3 1]);
    title('Power spectra, average over channels');
    legend('before','after'); legend boxoff;
    set(gca,'ygrid','on'); pbaspect([1.618 1 1]);

    nexttile;
    nt_spect_plot(data-clean,NPTS,0,NPTS,DATA(1).srate);
    title(['Noise power (removed), N = ' num2str(length(KREMOVE))]);
    set(gca,'ygrid','on'); pbaspect([1.618 1 1]);

    % plotX=30; plotY=15;
    % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    % print(fh, fullfile(DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_peaks_1']), '-dtiff', '-r200'); close(fh);
    save_figure(fh, DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_peaks_2'], [30 15]);

else
    fprintf('Great! No additional peaks found in the spectra.\n');

    % Log / Report
    for i = 1:NBLK
        DATA(i).ALSUTRECHT.otherPeakCleaning.peaks = NaN;
        DATA(i).ALSUTRECHT.otherPeakCleaning.N     = NaN;
        DATA(i).ALSUTRECHT.otherPeakCleaning.p1    = NaN;
    end
end

end