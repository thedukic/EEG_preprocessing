function EEG = denoise_emg(EEG)
% Denoising Source Separation
fprintf('Doing DSS to remove muscle activity...\n');
NREMOVE = 1;

% EMG signals
[b, a] = butter(10,35/(EEG.srate/2),'high');

% Covariance matrices of full band (c0) and filtered (c1)
chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
data = EEG.data(chaneeg,:)';
[c0,c1] = nt_bias_filter(data,b,a);

% DSS matrix
[todss, pwr0, pwr1] = nt_dss0(c0,c1,[],[]);
p1 = pwr1 ./ pwr0;

% % Finds kappa
% d = -diff(log10(p1));
% d = d./std(d);
% NREMOVE = find(d > 4,1,'first');
%
% fh = figure; hold on;
% plot(log10(p1),'o'); plot(NREMOVE*[1 1],ylim); axis tight;

% DSS components
z = nt_mmat(data,todss);

% Regress them out
clean = nt_tsr(data, z(:,1:NREMOVE));

% Return
EEG.data(chaneeg,:,:) = permute(reshape(clean,EEG.pnts,EEG.trials,sum(chaneeg)),[3 1 2]);

% =========================================================================
% Plot
NPTS = 512;

fh = figure;
th = tiledlayout(2,3);
th.TileSpacing = 'compact'; th.Padding = 'compact';

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

% Save
plotX=20; plotY=10;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_dss_muscle']),'-dtiff','-r300');
close(fh);

% figure;
% for i = 1:NREMOVE
%     mytopoplot(todss(:,i),[],'',nexttile); % colorbar;
% end

end