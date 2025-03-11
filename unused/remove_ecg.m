function EEG = remove_ecg(EEG,EXT)

chanecg = find(strcmp({EXT.chanlocs.labels},'ECG'));
if ~any(EEG.ALSUTRECHT.ica.ICsMostLikelyHeart) && any(chanecg)
    % Detect ECG/QRS
    [ECGmask, ECGepoch, ~, ~, pulsEstimate] = detect_ecg(EXT,[-250 450],'Recorded');

    % % mask to emphasize eyeblink intervals
    % ecgchan = EXT.data(chanecg,:)';
    % ECGmask = mean(ecgchan.^2,2);
    % quantile = 0.8;
    % tmp = sort(ECGmask);
    % ECGmask = min(ECGmask, tmp(round(size(ECGmask,1)*quantile))); % avoid extreme weight

    % % Filter EEG
    % [blECG, alECG] = butter(4,20/(EEG.srate/2),'low');
    % [bhECG, ahECG] = butter(4,1/(EEG.srate/2),'high');
    % assert(isstable(blECG,alECG));
    % assert(isstable(bhECG,ahECG));

    EEGdata = EEG.data;
    % EEGdata = do_filteringcore(blECG,alECG,EEGdata,EEG.event,EEG.srate);
    % EEGdata = do_filteringcore(bhECG,ahECG,EEGdata,EEG.event,EEG.srate);

    EEGdata = nt_demean(EEGdata');
    c0 = nt_cov(EEGdata);
    c1 = nt_cov(nt_demean(EEGdata(ECGmask,:)));
    [todss, pwr0, pwr1] = nt_dss0(c0,c1);

    p1 = pwr1./pwr0;
    figure; plot(p1,'.-'); title('ECG DSS'); ylabel('Score'); xlabel('Component');

    z = nt_mmat(EEGdata, todss);

    figure; hold on;
    x = EXT.data(chanecg,1:5*EEG.srate);
    x = x ./ max(abs(x));
    plot(x);
    x = z(1:1:5*EEG.srate);
    x = x ./ max(abs(x));
    plot(x);
    x = double(ECGmask(1:1:5*EEG.srate));
    plot(x);

    NREMOVE = 1;
    clean = nt_tsr(EEGdata, z(:,1:NREMOVE));

    % Plot
    fh = figure;
    th = tiledlayout(2,3);
    th.TileSpacing = 'compact'; th.Padding = 'compact';
    % plot bias score
    nexttile;
    plot(p1,'.-'); xlabel('Component'); ylabel('Score'); title('DSS scores');
    axis tight; xlim([0 50]); pbaspect([1.618 1 1]);

    % plot spectra of DSS components
    NPTS = 1024;
    nexttile;
    nt_spect_plot2(nt_normcol(z(:,1:25)),NPTS,0,NPTS,EEG(1).srate);
    title('Spectra of first 25 DSS components'); ylabel('Component'); xlabel('Frequency (Hz)');
    pbaspect([1.618 1 1]); colormap(brewermap(128,'BrBG'));

    % plot the first DSS weights
    fromdss = pinv(todss);
    % mytopoplot(todss(:,1),[],['DSS1 score: ' num2str(round(p1(1)))],nexttile); % colorbar;
    mytopoplot(fromdss(1,:),[],['DSS1 score: ' num2str(round(p1(1)))],nexttile); % colorbar;

    % plot spectra of data before and after removal of the peak component(s)
    nexttile; hold on;
    nt_spect_plot(EEGdata,NPTS,0,NPTS,EEG(1).srate);
    nt_spect_plot(clean,NPTS,0,NPTS,EEG(1).srate);
    nt_linecolors([],[3 1]);
    title('Power spectra, average over channels');
    legend('before','after'); legend boxoff;
    set(gca,'ygrid','on'); pbaspect([1.618 1 1]);

    nexttile;
    nt_spect_plot(EEGdata-clean,NPTS,0,NPTS,EEG(1).srate);
    title(['Noise power (removed), N = ' num2str(NREMOVE)]);
    set(gca,'ygrid','on'); pbaspect([1.618 1 1]);

end