function EEG = reduce_spectrapeaks(EEG)
% Needs imporvements:
% 1. The minimum peak prominance
% 2. The check how well DSS is done, like ratios DSS1/DSS2 > X

% Use DSS to isolate the peaks
fprintf('\nRemoving possible additional peaks from the spectrum...\n');

% First components are most dominated by these peaks
NREMOVE = 2;
% winSizeCompleteSpectrum = 20; % [s]
%
% chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
% data = cat(2,EEG(:).data);
% NPTS = size(data,2);
NBLK = length(EEG);

% % we want at least 8 segments for proper usage of pwelch
% if winSizeCompleteSpectrum*EEG(1).srate > NPTS/8
%     winSizeCompleteSpectrum = floor(NPTS/8/EEG(1).srate);
%     warning('Dataset is short. Adjusted window size for whole data set spectrum calculation to be 1/8 of the length.')
% end

% % Make it 20s long -> 0.05 Hz freq resolution
% [NCHN,NPTSALL] = size(data);
% NPTS = winSizeCompleteSpectrum * EEG(1).srate;
% NTRL = floor(NPTSALL/NPTS);
% data = reshape(data(:,1:NTRL*NPTS),NCHN,NPTS,NTRL);
%
% % Compute power spectra
% [NCHN, NPTS, NTRL] = size(data);
% psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);
%
% for i = 1:NTRL
%     [psdspectra(:,:,i), freq] = pwelch(data(:,:,i)',NPTS,0,NPTS,EEG(1).srate);
% end
%
% % Average the spectra
% psdspectra = mean(psdspectra,3);
% psdspectra = psdspectra(:,chaneeg);
% psdspectra = log10(mean(psdspectra,2));

[psdspectra, freq, chaneeg] = estimate_power(EEG,'speaks');
psdspectra = log10(mean(psdspectra,2));
NPTS = 2 * (length(freq)-1);

% % Identify the peaks which are always above 50 Hz
% freqMin = 52;
% freqMask = freq > freqMin;
%
% % Initial step
% [qrspeaks, locs, ~, proms] = findpeaks(psdspectra(freqMask),freq(freqMask));
% % Guided detection
% promFinal = quantile(proms,0.98);
% [qrspeaks, locs] = findpeaks(psdspectra(freqMask),freq(freqMask),'MinPeakProminence',promFinal);

% ALS34280: 1.65-1.7 Hz harmonics?
locs = [1.65 3.35 5.05 6.7 8.4 10.1 11.75 13.45 15.15]';

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
    p1z = zscore(p1);
    assert(p1z(1) > 5);
    assert(p1(1)/p1(2) > 2);

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
        EEG(i).ALSUTRECHT.otherPeakCleaning.N = NREMOVE;
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
    mytopoplot(todss(:,1),[],['DSS1 score: ' num2str(round(p1(1),2))],nexttile); % colorbar;

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
    title('Noise power (removed)');
    set(gca,'ygrid','on'); pbaspect([1.618 1 1]);

    plotX=30; plotY=15;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_otherpeakremoval']),'-dtiff','-r300');
    close(fh);

else
    fprintf('Great! No additional peaks found in the spectra.\n');

    % Log / Report
    for i = 1:NBLK
        EEG(i).ALSUTRECHT.otherPeakCleaning.peaks = NaN;
        EEG(i).ALSUTRECHT.otherPeakCleaning.N = NREMOVE;
    end
end

end