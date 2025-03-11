function [EEG, flagREDO] = report_leftovers(EEG,EXT,run,cfg)
%
% Based on the RELAX toolbox
% The code might not work very well
%

fprintf('\n================================\n');
fprintf('Detecting leftovers\n');
fprintf('================================\n');

%% ========================================================================
% [psdspectra, freq] = estimate_power(EEG,'freport');
%
% maskFreq  = freq>30 & freq<70;
% psdspectraGamma = mean(psdspectra(maskFreq,:),1);
% psdspectraGamma = zscore(psdspectraGamma);
% badElectrodes   = psdspectraGamma > 6;
%
% figure; tiledlayout(1,2);
% nexttile; plot(freq,log10(psdspectra));
% mytopoplot(psdspectraGamma,badElectrodes,'',nexttile); colorbar;
%
% chanLabelsEEG = {EEG.chanlocs(strcmp({EEG.chanlocs.type},'EEG')).labels};
% badElectrodes = chanLabelsEEG(badElectrodes);
% EEG.ALSUTRECHT.leftovers.badElectrodes = badElectrodes;
%
% % Interpolate bad electrodes
% if ~isempty(badElectrodes)
%     EEG = pop_select(EEG,'nochannel',badElectrodes);
%     EEG = pop_interp(EEG,chanlocs,'spherical');
% end

%% ========================================================================
fprintf('Checking muscle activity leftovers...\n');
muscleSlopeThreshold = cfg.bch.muscleSlopeThreshold;
muscleSlopeDuration  = cfg.bch.muscleSlopeTime;

% Estimate log-log power spectra
slopesChannelsxEpochs = detect_emg(EEG,cfg.bch);
[NCHANEEG, NTRL] = size(slopesChannelsxEpochs);

% Strong slow drifts are reflected as very steep negative slopes of the power spectrum
badchn = sum(slopesChannelsxEpochs > muscleSlopeThreshold,2);
badchn = badchn ./ NTRL;
badElectrodes = {EEG.chanlocs(find(badchn > muscleSlopeDuration)).labels};

% The following replaces all values that aren't above the muscle
% threshold with NaN, then sums the values that are above the threshold
% for each epoch to allow identification of the worst epochs:
slopesChannelsxEpochs(slopesChannelsxEpochs < muscleSlopeThreshold) = NaN;

% Shift the baseline of the values to the muscleSlopeThreshold, so that all
% muscle artifacts can be ranked in severity of EMG starting from 0 (least
% severe) and moving more positive as more severe:
slopesChannelsxEpochs = slopesChannelsxEpochs-muscleSlopeThreshold;

% Sum muscle slopes across all channels that show slopes above the
% threshold. This gives an indication of how badly each epoch is
% affected by muscle activity, with more affected electrodes within
% the epoch providing higher values:
slopesEpochs = sum(slopesChannelsxEpochs,1,'omitnan');

% Threshold = 0, because all slope values have had the threshold subtracted from them (so the threshold is now 0)
proportionOfDataShowingMuscleActivityTotal = mean(slopesEpochs > 0);

% Log
fprintf('Total amount of leftover muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Leftovers: muscle artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle log(7-75Hz) slope threshold: %1.2f\n',muscleSlopeThreshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of leftover muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);

EEG.ALSUTRECHT.leftovers.muscle1 = proportionOfDataShowingMuscleActivityTotal;

%% ========================================================================
fprintf('\nChecking eye blink letovers...\n');

fh = figure;
th = tiledlayout(2,2);
th.TileSpacing = 'compact'; th.Padding = 'compact';

chanlocs = readlocs('biosemi128_eeglab.ced');
myCmap1 = brewermap(128,'*RdBu');
myCmap2 = brewermap(128,'BuPu'); % BrBG

% Select only EEG + VEOG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
dataeeg  = EEG.data(chaneeg,:);

% =========================================================================
% Detect eye blinks
blinkLenght = 500;
[~, eyeBlinksEpochs, BlinkMaxLatency, dataeog, treshold] = detect_veog(EXT,blinkLenght,false);

% Find multi-blinks
multiBlink = detect_multiblinks(eyeBlinksEpochs,0);
eyeBlinksEpochs(multiBlink,:) = [];
NTRL1 = size(eyeBlinksEpochs,1);

if NTRL1 > 0
    L = mode(diff(eyeBlinksEpochs'))+1;
    timeBlink0 = (0:L-1)./EEG.srate*1000;
    timeBlink1  = timeBlink0-blinkLenght;

    dataeegepoched = NaN(NCHANEEG,L,NTRL1);
    dataeogepoched = NaN(L,NTRL1);
    for i = 1:NTRL1
        dataeegepoched(:,:,i) = dataeeg(:,eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
        dataeogepoched(:,i)   = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    end

    baselineTime = timeBlink0(end) * [0.05 0.95];
    timesel = timeBlink0<baselineTime(1) | timeBlink0>baselineTime(2);
    dataeegepoched = dataeegepoched - mean(dataeegepoched(:,timesel,:),2);
    % dataeegepoched = dataeegepoched - mean(dataeegepoched,2);
    dataeegepoched = dataeegepoched - mean(dataeegepoched,1);

    % 1. Plot
    % maskChanBlink = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);
    % dataeegplot = dataeegepoched(maskChanBlink,:,:);
    % dataeegplot = squeeze(mean(dataeegplot,1))';
    % [h,p,ci,stats] = ttest(dataeegplot);
    % dataeegplot = stats.tstat;
    % % for i = 1:size(dataeegplot,2)
    % %     [p(i),h,stats] = signrank(dataeegplot(:,i));
    % %     statszval(i) = stats.zval;
    % % end
    % % dataeegplot = statszval;
    %
    % th = nexttile; hold on;
    % plot(timeBlink1,zeros(1,length(timeBlink1)),'Color','k');
    % plot(timeBlink1,dataeegplot,'LineWidth',1.2);
    % maskTmp = p<0.01;
    % scatter(timeBlink1(maskTmp),dataeegplot(maskTmp));
    % axis tight; ylim([-10 25]); pbaspect([1.618 1 1]);
    % ylabel('t-test'); title(['Frontal electrodes blink leftovers, N = ' num2str(NTRL)]);

    % 2. Plot
    nexttile(1); hold on;
    plot(timeBlink1,dataeogepoched,'LineWidth',1.2);
    plot(timeBlink1,treshold*ones(size(timeBlink1)),'Color','k');
    dataCmapTmp = brewermap(NTRL1,'Spectral'); % BuGn
    set(gca,'ColorOrder',[0 0 0; dataCmapTmp]);
    title(['Detected blinks, N = ' num2str(NTRL1)]);
    pbaspect([1.618 1 1]); ylabel('EOG amplitude (\muV)');

    dataeegepoched = squeeze(mean(dataeegepoched,2));
    [h,p,ci,stats] = ttest(dataeegepoched');
    % dataeegplot = mean(dataeegepoched,2);

    tstatMean = mean(stats.tstat);

    nexttile(2); hold on;
    % plotThis = dataeegplot;
    plotThis = stats.tstat;
    % plotThis = -log10(p);
    topoplot(plotThis,chanlocs,'maplimits',max(abs(plotThis))*[-1 1],'headrad',0.5,'colormap',myCmap1,'whitebk','on','electrodes','off','style','map','shading','interp');
    title('EEG timelocked to the eye blinks');
    hcb = colorbar;
    hcb.Title.String = "T-value";

    % plot(timeBlink1,zeros(1,length(timeBlink1)),'Color','k');
    % plot(timeBlink1,dataeegplot,'LineWidth',1.2);
    % maskTmp = p<0.01;
    % scatter(timeBlink1(maskTmp),dataeegplot(maskTmp));
    % axis tight; ylim([-10 25]); pbaspect([1.618 1 1]);
    % ylabel('t-test'); title(['Frontal electrodes blink leftovers, N = ' num2str(NTRL)]);
else
    warning('No data to make an estimate of blink leftovers...');
    stats.tstat = NaN;
    tstatMean = NaN;
end

% =========================================================================
% Detect eye blinks
% Not ideal, the code does not care about boundary events
blinkLenght = 2000;
[~, eyeBlinksEpochs, BlinkMaxLatency, dataeog, treshold] = detect_veog(EXT,blinkLenght,false);

% Find multi-blinks
multiBlink = detect_multiblinks(eyeBlinksEpochs,0);

fprintf('Number of detected blinks is %d.\n',size(eyeBlinksEpochs,1));
fprintf('Number of detected multiple blinks within each evalulation window is %d.\n',sum(multiBlink));

eyeBlinksEpochs(multiBlink,:) = [];
NTRL2 = size(eyeBlinksEpochs,1);

L = mode(diff(eyeBlinksEpochs'))+1;
timeBlink0 = (0:L-1)./EEG.srate*1000;
timeBlink1  = timeBlink0 - blinkLenght;

if NTRL2 > 0
    dataeegepoched = NaN(NCHANEEG,L,NTRL2);
    dataeogepoched = NaN(L,NTRL2);
    for i = 1:NTRL2
        dataeegepoched(:,:,i) = dataeeg(:,eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
        dataeogepoched(:,i)   = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    end

    % yrange = 1.05*[min(dataeogepoched(:)), max(dataeogepoched(:))];
    % fh = figure;
    % th = tiledlayout(1,3);
    % th.TileSpacing = 'compact'; th.Padding = 'compact';

    % figure;
    % hold on;
    % for i = 1:NTRL
    %     cla;
    %     plot(timeBlink1,treshold*ones(size(timeBlink1)),'Color','k');
    %     plot(timeBlink1,dataeogepoched(:,i),'LineWidth',1.2);
    %     pause;
    % end

    % 2. Plot
    nexttile(3); hold on;
    plot(timeBlink1,treshold*ones(size(timeBlink1)),'Color','k');
    plot(timeBlink1,dataeogepoched,'LineWidth',1.2);
    dataCmapTmp = brewermap(NTRL2,'Spectral'); % BuGn
    set(gca,'ColorOrder',[0 0 0; dataCmapTmp]);
    title(['Detected blinks, N = ' num2str(NTRL2)]);
    pbaspect([1.618 1 1]); ylabel('EOG amplitude (\muV)');

    % plot([500 500],yrange,'Color',[0 0 0]);
    % plot([3500 3500],yrange,'Color',[0 0 0]);
    % plot([1500 1500],yrange,'Color',[0 0 0]);
    % plot([2500 2500],yrange,'Color',[0 0 0]);
    % axis tight;

    % RELAX code
    % Calculate the absolute difference between the EEG.times and the timepoints we need:
    absDiff_500ms = abs(timeBlink0 - 500);
    % Find the minimum absolute difference
    minDiff_500ms = min(absDiff_500ms(:));
    % Find the indices of the closest number
    [~, col_500ms] = find(absDiff_500ms == minDiff_500ms);
    % 3500 ms:
    absDiff_3500ms = abs(timeBlink0 - 3500);
    % Find the minimum absolute difference
    minDiff_3500ms = min(absDiff_3500ms(:));
    % Find the indices of the closest number
    [~, col_3500ms] = find(absDiff_3500ms == minDiff_3500ms);
    % 1500 ms:
    absDiff_1500ms = abs(timeBlink0 - 1500);
    % Find the minimum absolute difference
    minDiff_1500ms = min(absDiff_1500ms(:));
    % Find the indices of the closest number
    [~, col_1500ms] = find(absDiff_1500ms == minDiff_1500ms);
    % 2500 ms:
    absDiff_2500ms = abs(timeBlink0 - 2500);
    % Find the minimum absolute difference
    minDiff_2500ms = min(absDiff_2500ms(:));
    % Find the indices of the closest number
    [~, col_2500ms] = find(absDiff_2500ms == minDiff_2500ms);
    % 4000 ms:
    col_4000ms = L;

    % disp(([1 col_500ms col_1500ms col_2500ms col_3500ms col_4000ms]-1)/EEG.srate);

    % Baseline correct data
    dataeegepoched = dataeegepoched - mean(dataeegepoched(:,[1:col_500ms, col_3500ms:col_4000ms],:),2);

    % Convert to absolute values
    absolutevaluesblink = abs(dataeegepoched);

    % Calculate
    BlinkAmplitudeRatioAllEpochs = NaN(NCHANEEG,NTRL2);
    for i = 1:NTRL2
        BlinkAmplitudeRatioAllEpochs(:,i) = mean(absolutevaluesblink(:,col_1500ms:col_2500ms,i),2) ./ mean(absolutevaluesblink(:,[1:col_500ms, col_3500ms:col_4000ms],i),2);
    end
    BlinkAmplitudeRatio = mean(BlinkAmplitudeRatioAllEpochs,2)-1;
    BlinkAmplitudeRatioMean = mean(BlinkAmplitudeRatio)*100;

    % 3. Plot
    nexttile(4);
    % mask = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);
    % maxBlinkRatio = max(BlinkAmplitudeRatio);
    maxBlinkRatio = prctile(BlinkAmplitudeRatio,95);
    maxBlinkRatio = max(maxBlinkRatio, 0.2); % minimum is this value

    topoplot(BlinkAmplitudeRatio,chanlocs,'maplimits',[0 maxBlinkRatio],'headrad',0.5,'colormap',myCmap2,'whitebk','on','electrodes','off','style','map','shading','interp'); % ,'emarker2',{find(mask),'d','k',10,1}
    title({'Mean blink amplitude leftover', [num2str(round(BlinkAmplitudeRatioMean)) '%']});
    hcb = colorbar;
    hcb.Title.String = "%";

else
    warning('No data to make an estimate of blink leftovers...');
    BlinkAmplitudeRatio = NaN;
    BlinkAmplitudeRatioMean = NaN;
end

% Save
plotX=25; plotY=25;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_leftovers_' num2str(run)]),'-dtiff','-r300');
close(fh);

% Log
fprintf('Average amount of leftover eye blink artifact: %1.0f%%\n', BlinkAmplitudeRatioMean);
fprintf('Average amount of leftover eye blink artifact: %1.1f T-stat\n', tstatMean);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Leftovers: eye blink artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Average amount of leftover eye blink artifact: %1.0f%%\n', BlinkAmplitudeRatioMean);
fprintf(EEG.ALSUTRECHT.subject.fid,'Average amount of leftover eye blink artifact: %1.1f T-stat\n', tstatMean);

EEG.ALSUTRECHT.leftovers.blinksRatio = BlinkAmplitudeRatio;
EEG.ALSUTRECHT.leftovers.blinksTstat = stats.tstat;

%% ========================================================================
load(fullfile(EEG.ALSUTRECHT.subject.mycodes,'files','Blinkweights'),'Blinkweights');
maskChanBlink = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);

% blinkchans = {'C14','C15','C16','C17','C18','C19','C27','C28','C29'};
% maskChanBlink = ismember({EEG.chanlocs(:).labels},blinkchans);

blinkTresholdCorr   = 0.7;
blinkTreshold2Tstat = 4;
blinkTreshold2Perc  = 0.1;

if (~isnan(tstatMean) && NTRL1>5) && (~isnan(BlinkAmplitudeRatioMean) && NTRL2>5)
    % 1. Compare the leftover maps to the blink IC
    % But do not do it if there are not enough eyeblinks detected
    BlinkAmplitudeTstatNorm = stats.tstat' ./ norm(stats.tstat);
    BlinkAmplitudeRatioNorm = BlinkAmplitudeRatio ./ norm(BlinkAmplitudeRatio);
    BlinkweightsNorm        = Blinkweights ./ norm(Blinkweights);

    BlinkAmplitudeTstatNorm = BlinkAmplitudeTstatNorm - mean(BlinkAmplitudeTstatNorm);
    BlinkAmplitudeRatioNorm = BlinkAmplitudeRatioNorm - mean(BlinkAmplitudeRatioNorm);
    BlinkweightsNorm        = BlinkweightsNorm - mean(BlinkweightsNorm);

    % Smoothen to remove noise
    BlinkAmplitudeTstatNorm = estimate_invlaplacian(BlinkAmplitudeTstatNorm,EEG.chanlocs,1);
    BlinkAmplitudeRatioNorm = estimate_invlaplacian(BlinkAmplitudeRatioNorm,EEG.chanlocs,1);

    corrMatTstat = abs(corr(BlinkAmplitudeTstatNorm, BlinkweightsNorm));
    corrMatPerc  = abs(corr(BlinkAmplitudeRatioNorm, BlinkweightsNorm));

    % figure;
    % mytopoplot(BlinkAmplitudeRatioNorm,[],[],nexttile);
    % mytopoplot(BlinkweightsNorm,[],[],nexttile);

    % Lowered to capture imperfect leftovers
    flagCorrTstat = corrMatTstat > blinkTresholdCorr;
    flagCorrPercent = corrMatPerc  > blinkTresholdCorr;

    % 2.Frontal electrodes should have high blink leftover
    meanFrontalBlinkLeftoverTstat = mean(stats.tstat(maskChanBlink));
    flagMeanTstat = meanFrontalBlinkLeftoverTstat > blinkTreshold2Tstat;

    meanFrontalBlinkLeftoverPerc = mean(BlinkAmplitudeRatio(maskChanBlink));
    flagMeanPercent = meanFrontalBlinkLeftoverPerc>blinkTreshold2Perc;

    % Combine
    flagREDO = (flagCorrTstat | flagCorrPercent) & (flagMeanTstat | flagMeanPercent);

    % Report
    fprintf('Blink leftover map correlations are %1.2f and %1.2f.\n',corrMatTstat,corrMatPerc);
    fprintf('Average frontal leftovers are %1.2f (T-stat) and %1.2f (%%).\n',meanFrontalBlinkLeftoverTstat,meanFrontalBlinkLeftoverPerc);

    EEG.ALSUTRECHT.leftovers.flagCorrTstat   = flagCorrTstat;
    EEG.ALSUTRECHT.leftovers.flagCorrPercent = flagCorrPercent;
    EEG.ALSUTRECHT.leftovers.flagMeanTstat   = flagMeanTstat;
    EEG.ALSUTRECHT.leftovers.flagMeanPercent = flagMeanPercent;

elseif ~isnan(tstatMean) && NTRL1>5
    % 1. Compare the leftover map to the blink IC
    % But do not do it if there are not enough eyeblinks detected
    BlinkAmplitudeTstatNorm = stats.tstat' ./ norm(stats.tstat);
    BlinkweightsNorm        = Blinkweights ./ norm(Blinkweights);

    BlinkAmplitudeTstatNorm = BlinkAmplitudeTstatNorm - mean(BlinkAmplitudeTstatNorm);
    BlinkweightsNorm        = BlinkweightsNorm - mean(BlinkweightsNorm);

    % Smoothen to remove noise
    BlinkAmplitudeTstatNorm = estimate_invlaplacian(BlinkAmplitudeTstatNorm,EEG.chanlocs,1);

    corrMatTstat = abs(corr(BlinkAmplitudeTstatNorm, BlinkweightsNorm));
    % figure;
    % mytopoplot(BlinkAmplitudeTstatNorm,[],[],nexttile);
    % mytopoplot(icawinvSmooth,[],[],nexttile);
    % mytopoplot(BlinkweightsNorm,[],[],nexttile);

    % Lowered to capture imperfect leftovers
    flagCorrTstat = corrMatTstat > blinkTresholdCorr;

    % 2.Frontal electrodes should have high blink leftover
    meanFrontalBlinkLeftoverTstat = mean(stats.tstat(maskChanBlink));
    flagMeanTstat = meanFrontalBlinkLeftoverTstat > blinkTreshold2Tstat;

    % Combine
    flagREDO = flagCorrTstat & flagMeanTstat;

    % Report
    fprintf('Blink leftover map correlation is %1.2f.\n',corrMatTstat);
    fprintf('Average frontal leftover is %1.2f (T-stat).\n',meanFrontalBlinkLeftoverTstat);

    EEG.ALSUTRECHT.leftovers.flagCorrTstat   = flagCorrTstat;
    EEG.ALSUTRECHT.leftovers.flagCorrPercent = NaN;
    EEG.ALSUTRECHT.leftovers.flagMeanTstat   = flagMeanTstat;
    EEG.ALSUTRECHT.leftovers.flagMeanPercent = NaN;
end

% %%
% EEG.ALSUTRECHT.JD.flagREDO = flagREDO;
% if flagREDO
%     warning('Blink leftover is too big, joint decorrelation will be done...\n');
%     EEG = do_jointdecorrelation(EEG,EXT,cfg);
%     [y,z,mask] = nt_eyeblink(EEG.data',maskChanBlink,1,EEG.srate);
% end

% Remove (not needed)
EEG.icaact = [];

end

% =========================================================================
% Helper functions
% =========================================================================
function multiBlink = detect_multiblinks(eyeBlinksEpochs,overlap)
% Find multi-blinks
% Function to create the range specified by each row
% Apply the function to each row of X

% Narrow down the epoch to get more of them
% It is okay if they overlap 1s (256 samples)
eyeBlinksEpochs(:,1) = eyeBlinksEpochs(:,1) + overlap;
eyeBlinksEpochs(:,2) = eyeBlinksEpochs(:,2) - overlap;

createRange = @(row) row(1):row(2);
rangesCell = arrayfun(@(i) createRange(eyeBlinksEpochs(i,:)), 1:size(eyeBlinksEpochs,1), 'UniformOutput', false);

NTRL = length(rangesCell);
multiBlink = false(1,NTRL);
for i = 1:NTRL
    concatenatedRanges    = rangesCell;
    concatenatedRanges(i) = [];
    multiBlink(i) = any(ismember(rangesCell{i},[concatenatedRanges{:}]));
end

end