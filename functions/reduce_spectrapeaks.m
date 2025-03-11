function EEG = reduce_spectrapeaks(EEG)
% Needs imporvements:
% 1. The minimum peak prominance
% 2. The check how well DSS is done, like ratios DSS1/DSS2 > X

% Use DSS to isolate the peaks
fprintf('\n================================\n');
fprintf('Removing additional peaks from the spectrum\n');
fprintf('================================\n');

% First components are most dominated by these peaks
NREMOVE = 1;
NBLK = length(EEG);

[psdspectra, freq, chaneeg] = estimate_power(EEG,'speaks');
psdspectra = log10(mean(psdspectra,2));
NPTS = 2 * (length(freq)-1);

% Special cases
if strcmpi(EEG(1).ALSUTRECHT.subject.id,'ALS34280') && strcmpi(EEG(1).ALSUTRECHT.subject.visit,'T1')
    % 1.65-1.7 Hz harmonics?
    locs = [1.65 3.35 5.05 6.7 8.4 10.1 11.75 13.45 15.15]';

elseif strcmpi(EEG(1).ALSUTRECHT.subject.id,'ALS36104') && strcmpi(EEG(1).ALSUTRECHT.subject.visit,'T1')
    % 1.35 Hz harmonics?
    % locs = [1.35 2.7 4.1 5.4 6.75 8 9.5 10.75 12.2 13.5 14.9]';
    % locs = [1.35 2.7 4 4.15 5.45 6.75 6.9 8.15 12 12.5 13.5 14.9 16.3 17.6 19]';
    % locs = [1.35 2.7 6.7]';
    % locs = round((1:5)' * 1.35,2);
    locs = [1.35 2.7 4.1 5.4 6.75]';
    NREMOVE = 4;

else
    % % Identify the peaks which are always above 50 Hz
    % freqMask = freq>50 & freq<120;
    % 
    % % Initial step
    % [qrspeaks, locs, ~, proms] = findpeaks(psdspectra(freqMask),freq(freqMask));
    % % Guided detection
    % promFinal = quantile(proms, 0.98);
    % [qrspeaks, locs, ~, proms] = findpeaks(psdspectra(freqMask),freq(freqMask),'MinPeakProminence',promFinal);
    % 
    % % Prevent removing peaks that are not actual peaks/noise
    % if max(proms) < 0.25
    %     locs = [];
    % end
    locs = [];
end

% Rounding errors
locs = round(locs,2);

% If any peaks found
if ~isempty(locs)
    fprintf('Peaks found at:\n');
    fprintf('f = %1.2f Hz\n',locs);

    % Covariance matrices of full band (c0) and filtered at the peaks (c1)
    data = cat(2,EEG(:).data)';
    [c0, c1] = nt_bias_fft(data, locs'/EEG(1).srate, NPTS);

    % DSS matrix
    [todss, pwr0, pwr1] = nt_dss0(c0,c1,[],[]);
    p1 = pwr1 ./ pwr0;

    % Check quality
    % p1z = zscore(p1);
    % assert(p1z(1) > 5);
    % assert(p1(1)/p1(2) > 2);

    % figure;
    % for i = 1:NREMOVE
    %     mytopoplot(todss(:,i),[],'',nexttile); colorbar;
    % end

    % DSS components
    z = nt_mmat(data,todss);

    % Regress them out
    clean = nt_tsr(data,z(:,1:NREMOVE));

    % Split back the blocks
    EEG2 = make_rsmasks(EEG);
    assert(size(clean,1) == size(EEG2(1).ALSUTRECHT.blockinfo.rs_mask,2));

    for i = 1:NBLK
        EEG(i).data = clean(EEG2(i).ALSUTRECHT.blockinfo.rs_mask(i,:),:)';
        assert(size(EEG(i).data,2) == EEG(i).pnts);
    end
    EEG = eeg_checkset(EEG);

    % Log / Report
    for i = 1:NBLK
        EEG(i).ALSUTRECHT.otherPeakCleaning.peaks = locs;
        EEG(i).ALSUTRECHT.otherPeakCleaning.N     = NREMOVE;
        EEG(i).ALSUTRECHT.otherPeakCleaning.p1    = p1;
    end
    fprintf('Done! Number of components removed: %d\n',NREMOVE);

    % =====================================================================
    % Plot
    fh = figure;
    th = tiledlayout(2,3);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    freq = round(freq,2);
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
    nt_spect_plot2(nt_normcol(z(:,1:25)),NPTS,0,NPTS,EEG(1).srate);
    title('Spectra of first 25 DSS components'); ylabel('Component'); xlabel('Frequency (Hz)');
    pbaspect([1.618 1 1]); colormap(brewermap(sum(chaneeg),'BrBG'));

    % plot the first DSS weights
    fromdss = pinv(todss);
    % mytopoplot(todss(:,1),[],['DSS1 score: ' num2str(round(p1(1)))],nexttile); % colorbar;
    mytopoplot(fromdss(1,:),[],['DSS1 score: ' num2str(round(p1(1)))],nexttile); % colorbar;

    % plot spectra of data before and after removal of the peak component(s)
    nexttile; hold on;
    nt_spect_plot(data,NPTS,0,NPTS,EEG(1).srate);
    nt_spect_plot(clean,NPTS,0,NPTS,EEG(1).srate);
    nt_linecolors([],[3 1]);
    title('Power spectra, average over channels');
    legend('before','after'); legend boxoff;
    set(gca,'ygrid','on'); pbaspect([1.618 1 1]);

    nexttile;
    nt_spect_plot(data-clean,NPTS,0,NPTS,EEG(1).srate);
    title(['Noise power (removed), N = ' num2str(NREMOVE)]);
    set(gca,'ygrid','on'); pbaspect([1.618 1 1]);

    plotX=30; plotY=15;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_dss_peaks']),'-dtiff','-r300');
    close(fh);

else
    fprintf('Great! No additional peaks found in the spectra.\n');

    % Log / Report
    for i = 1:NBLK
        EEG(i).ALSUTRECHT.otherPeakCleaning.peaks = NaN;
        EEG(i).ALSUTRECHT.otherPeakCleaning.N     = NREMOVE;
        EEG(i).ALSUTRECHT.otherPeakCleaning.p1    = NaN;
    end
end

end