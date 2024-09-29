function EEG = report_leftovers(EEG,EXT,cfg)
%
% RELAX toolbox
% SDukic, July 2024
%

%% =========================================================================
fprintf('\nChecking muscle activity leftovers...\n');
muscleSlopeThreshold = cfg.bch.muscleSlopeThreshold;
muscleSlopeDuration  = 0.5;

% Estimate log-log power spectra
slopesChannelsxEpochs = dected_emg(EEG);
[NCHANEEG, NTRL] = size(slopesChannelsxEpochs);

% Strong slow drifts are reflected as very steep negative slopes of the power spectrum
badchn = sum(slopesChannelsxEpochs>muscleSlopeThreshold,2);
badchn = badchn./NTRL;
badElectrodes = {EEG.chanlocs(find(badchn>muscleSlopeDuration)).labels};

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
proportionOfDataShowingMuscleActivityTotal = mean(slopesEpochs>0);

% Log
fprintf('Total amount of leftover muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Leftovers: muscle artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle log(7-75Hz) slope threshold: %1.2f\n',muscleSlopeThreshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of leftover muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);

EEG.ALSUTRECHT.leftovers.muscle1 = proportionOfDataShowingMuscleActivityTotal;

%% =========================================================================
fprintf('\nChecking eye blink letovers...\n');

fh = figure;
th = tiledlayout(1,3);
th.TileSpacing = 'compact'; th.Padding = 'compact';

% Select only EEG + VEOG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
dataeeg  = EEG.data(chaneeg,:);

% =========================================================================
% Detect eye blinks
blinkLenght = 500;
[~, eyeBlinksEpochs] = detect_eog(EXT,blinkLenght,false);

% Find multi-blinks
multiBlink = detect_multiblinks(eyeBlinksEpochs);
eyeBlinksEpochs(multiBlink,:) = [];
NTRL = size(eyeBlinksEpochs,1);

L = mode(diff(eyeBlinksEpochs'))+1;
timeBlink0 = (0:L-1)./EEG.srate*1000;
timeBlink1  = timeBlink0-blinkLenght;

dataeegepoched = NaN(NCHANEEG,L,NTRL);
for i = 1:NTRL
    dataeegepoched(:,:,i) = dataeeg(:,eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    % dataeogepoched(:,i)   = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
end

timesel = timeBlink0<200 | timeBlink0>800;
dataeegepoched = dataeegepoched - mean(dataeegepoched(:,timesel,:),2);
dataeegepoched = dataeegepoched - mean(dataeegepoched,1);

% 1. Plot
maskChanBlink = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);
dataeegplot = dataeegepoched(maskChanBlink,:,:);
dataeegplot = squeeze(mean(dataeegplot,1))';
[h,p,ci,stats] = ttest(dataeegplot);
dataeegplot = stats.tstat;
% for i = 1:size(dataeegplot,2)
%     [p(i),h,stats] = signrank(dataeegplot(:,i));
%     statszval(i) = stats.zval;
% end
% dataeegplot = statszval;

th = nexttile; hold on;
plot(timeBlink1,zeros(1,length(timeBlink1)),'Color','k');
plot(timeBlink1,dataeegplot,'LineWidth',1.2);
maskTmp = p<0.01;
scatter(timeBlink1(maskTmp),dataeegplot(maskTmp));
axis tight; ylim([-10 25]); pbaspect([1.618 1 1]);
ylabel('t-test'); title(['Frontal electrodes blink leftovers, N = ' num2str(NTRL)]);

% =========================================================================
% Detect eye blinks
% Not ideal, the code does not care about boundary events
blinkLenght = 2000;
[~, eyeBlinksEpochs, BlinkMaxLatency, dataeog, treshold] = detect_eog(EXT,blinkLenght,false);

% Find multi-blinks
multiBlink = detect_multiblinks(eyeBlinksEpochs);

fprintf('Number of detected blinks is %d.\n',size(eyeBlinksEpochs,1));
fprintf('Number of detected multiple blinks within each evalulation window is %d.\n',sum(multiBlink));

eyeBlinksEpochs(multiBlink,:) = [];
NTRL = size(eyeBlinksEpochs,1);

L = mode(diff(eyeBlinksEpochs'))+1;
timeBlink0 = (0:L-1)./EEG.srate*1000;
timeBlink1  = timeBlink0-blinkLenght;

if NTRL>0
    dataeegepoched = NaN(NCHANEEG,L,NTRL);
    dataeogepoched = NaN(L,NTRL);
    for i = 1:NTRL
        dataeegepoched(:,:,i) = dataeeg(:,eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
        dataeogepoched(:,i)   = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    end

    % yrange = 1.05*[min(dataeogepoched(:)), max(dataeogepoched(:))];
    % fh = figure;
    % th = tiledlayout(1,3);
    % th.TileSpacing = 'compact'; th.Padding = 'compact';

    % 2. Plot
    nexttile; hold on;
    plot(timeBlink1,treshold*ones(size(timeBlink1)),'Color','k');
    plot(timeBlink1,dataeogepoched,'LineWidth',1.2);
    dataCmapTmp = brewermap(NTRL,'Spectral'); % BuGn
    set(gca,'ColorOrder',[0 0 0; dataCmapTmp]);
    title(['Detected blinks, N = ' num2str(NTRL)]);
    pbaspect([1.618 1 1]); ylabel('EOG amplitude (/muV)');

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
    BlinkAmplitudeRatioAllEpochs = NaN(NCHANEEG,NTRL);
    for i = 1:NTRL
        BlinkAmplitudeRatioAllEpochs(:,i) = mean(absolutevaluesblink(:,col_1500ms:col_2500ms,i),2) ./ mean(absolutevaluesblink(:,[1:col_500ms, col_3500ms:col_4000ms],i),2);
    end
    BlinkAmplitudeRatio = mean(BlinkAmplitudeRatioAllEpochs,2);
    BlinkAmplitudeRatioMean = (mean(BlinkAmplitudeRatio)-1)*100;

    % 3. Plot
    th = nexttile;
    % mask = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);
    chanlocs = readlocs('biosemi128_eeglab.ced'); myCmap = brewermap(128,'BuPu'); % BrBG
    topoplot(BlinkAmplitudeRatio,chanlocs,'headrad',0.5,'colormap',myCmap,'whitebk','on','electrodes','off','style','map','shading','interp'); % ,'emarker2',{find(mask),'d','k',10,1}

    % maxBlinkRatio = max(BlinkAmplitudeRatio);
    maxBlinkRatio = prctile(BlinkAmplitudeRatio,95);
    maxBlinkRatio = max(maxBlinkRatio,1.5); % minimum is this value
    clim(th,[1, maxBlinkRatio]); colorbar;
    title({'Mean blink amplitude leftover', [num2str(round(BlinkAmplitudeRatioMean)) '%']});

else
    warning('No data to make an estimate of blink leftovers...');
    BlinkAmplitudeRatio = NaN;
    BlinkAmplitudeRatioMean = NaN;
end

% Save
plotX=25; plotY=10;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_leftovers']),'-dtiff','-r300');
close(fh);

% Log
fprintf('Average amount of leftover eye blink artifact: %1.0f%%\n', BlinkAmplitudeRatioMean);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Leftovers: eye blink artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of leftover eye blink artifact: %1.0f%%\n', BlinkAmplitudeRatioMean);

EEG.ALSUTRECHT.leftovers.blinks = BlinkAmplitudeRatio;

end

function multiBlink = detect_multiblinks(eyeBlinksEpochs)
% Find multi-blinks
% Function to create the range specified by each row
% Apply the function to each row of X

% Narrow down the epoch to get more of them
% It is okay if they overlap 1s (256 samples)
eyeBlinksEpochs(:,1) = eyeBlinksEpochs(:,1)+256;
eyeBlinksEpochs(:,2) = eyeBlinksEpochs(:,2)-256;

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