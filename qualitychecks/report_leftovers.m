function EEG = report_leftovers(EEG,EXT,cfg)
%
% RELAX toolbox
% SDukic, July 2024
%

%% =========================================================================
fprintf('\nChecking muscle activity leftovers...\n');
muscleSlopeThreshold = cfg.bch.muscleSlopeThreshold;
muscleSlopeDuration  = 0.5;

% Select only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
% dataeeg = EEG.data(chaneeg,:,:);
dataeeg = EEG.data(chaneeg,:);

NCHANEEG = sum(chaneeg);

% Epoch into 1s (or maybe better into 1s with 0.5 overlap, but OK)
L = EEG.srate;
N = floor(size(dataeeg,2)/L);
dataeeg = reshape(dataeeg(:,1:N*L),NCHANEEG,L,N);

% Compute (log) power spectra (7-75 Hz)
[NCHN,NPTS,NTRL] = size(dataeeg);

psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);
for i = 1:NTRL
    [psdspectra(:,:,i),freq] = pwelch(dataeeg(:,:,i)',NPTS,0,NPTS,EEG.srate);
end

foi = [7 75];
frqmsk = freq>=foi(1) & freq<=foi(2);
psdspectra = permute(psdspectra,[1 3 2]);
logpow = log10(psdspectra(frqmsk,:,:));
logfoi = log10(freq(frqmsk));

% Fit linear regression to log-log data, and store the slope
slopesEpochsxChannels = NaN(NCHN,NTRL);
for i = 1:NCHN
    for j = 1:NTRL
        p = polyfit(logfoi,logpow(:,j,i),1);
        slopesEpochsxChannels(i,j) = p(1);
    end
end

% Strong slow drifts are reflected as very steep negative slopes of the power spectrum
badchn = sum(slopesEpochsxChannels>muscleSlopeThreshold,2);
badchn = badchn./NTRL;
badElectrodes = {EEG.chanlocs(find(badchn>muscleSlopeDuration)).labels};

% The following replaces all values that aren't above the muscle
% threshold with NaN, then sums the values that are above the threshold
% for each epoch to allow identification of the worst epochs:
sortingOutWorstMuscleEpochs = slopesEpochsxChannels;
sortingOutWorstMuscleEpochs(sortingOutWorstMuscleEpochs < muscleSlopeThreshold) = NaN;

% Shift the baseline of the values to the muscleSlopeThreshold, so that all
% muscle artifacts can be ranked in severity of EMG starting from 0 (least
% severe) and moving more positive as more severe:
sortingOutWorstMuscleEpochs = sortingOutWorstMuscleEpochs-muscleSlopeThreshold;

% Sum muscle slopes across all channels that show slopes above the
% threshold. This gives an indication of how badly each epoch is
% affected by muscle activity, with more affected electrodes within
% the epoch providing higher values:
sortingOutWorstMuscleEpochs = sum(sortingOutWorstMuscleEpochs,1,'omitnan');

% Threshold = 0, because all slope values have had the threshold subtracted from them (so the threshold is now 0)
templateMarkedForMuscleArtifacts(sortingOutWorstMuscleEpochs>0) = 1;
proportionOfDataShowingMuscleActivityTotal = mean(templateMarkedForMuscleArtifacts,'omitnan');

% Log
fprintf('Total amount of leftover muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Leftovers: muscle artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle log(7-75Hz) slope threshold: %1.2f\n',muscleSlopeThreshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of leftover muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);

EEG.ALSUTRECHT.leftovers.muscle = proportionOfDataShowingMuscleActivityTotal;

%% =========================================================================
fprintf('\nChecking eye blink letovers...\n');

% Select only EEG + VEOG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
dataeeg  = EEG.data(chaneeg,:);

% Detect eye blinks
[~, eyeBlinksEpochs, BlinkMaxLatency, dataeog, treshold] = detect_eog(EXT,2000);

% Find multi-blinks
% Function to create the range specified by each row
% Apply the function to each row of X
createRange = @(row) row(1):row(2);
rangesCell = arrayfun(@(i) createRange(eyeBlinksEpochs(i,:)), 1:size(eyeBlinksEpochs,1), 'UniformOutput', false);

NTRL = length(rangesCell);
multiBlink = false(1,NTRL);
for i = 1:NTRL
    concatenatedRanges    = rangesCell;
    concatenatedRanges(i) = [];
    multiBlink(i) = any(ismember(rangesCell{i},[concatenatedRanges{:}]));
end

fprintf('Number of detected blinks is %d.\n',NTRL);
fprintf('Number of detected multiple blinks within each evalulation window is %d.\n',sum(multiBlink));

eyeBlinksEpochs(multiBlink,:) = [];
NTRL = size(eyeBlinksEpochs,1);

L = mode(diff(eyeBlinksEpochs'))+1;
timeBlink = (0:L-1)./EEG.srate*1000-2000;

if NTRL>5
    dataeegepoched = NaN(NCHANEEG,L,NTRL);
    dataeogepoched = NaN(L,NTRL);
    for i = 1:NTRL
        dataeegepoched(:,:,i) = dataeeg(:,eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
        dataeogepoched(:,i)   = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    end

    % yrange = 1.05*[min(dataeogepoched(:)), max(dataeogepoched(:))];
    fh = figure;
    th = tiledlayout(1,2);
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    nexttile; hold on;
    plot(timeBlink,treshold*ones(size(timeBlink)),'Color','k');
    plot(timeBlink,dataeogepoched,'LineWidth',1.2);
    set(gca,'ColorOrder',[0 0 0; brewermap(NTRL,'BuGn')]);
    title(['N = ' num2str(NTRL)]);

    % plot([500 500],yrange,'Color',[0 0 0]);
    % plot([3500 3500],yrange,'Color',[0 0 0]);
    % plot([1500 1500],yrange,'Color',[0 0 0]);
    % plot([2500 2500],yrange,'Color',[0 0 0]);
    % axis tight;

    % RELAX code
    % Calculate the absolute difference between the EEG.times and the timepoints we need:
    absDiff_500ms = abs(timeBlink - 500);
    % Find the minimum absolute difference
    minDiff_500ms = min(absDiff_500ms(:));
    % Find the indices of the closest number
    [~, col_500ms] = find(absDiff_500ms == minDiff_500ms);
    % 3500ms:
    absDiff_3500ms = abs(timeBlink - 3500);
    % Find the minimum absolute difference
    minDiff_3500ms = min(absDiff_3500ms(:));
    % Find the indices of the closest number
    [~, col_3500ms] = find(absDiff_3500ms == minDiff_3500ms);
    % 1500ms:
    absDiff_1500ms = abs(timeBlink - 1500);
    % Find the minimum absolute difference
    minDiff_1500ms = min(absDiff_1500ms(:));
    % Find the indices of the closest number
    [~, col_1500ms] = find(absDiff_1500ms == minDiff_1500ms);
    % 2500ms:
    absDiff_2500ms = abs(timeBlink - 2500);
    % Find the minimum absolute difference
    minDiff_2500ms = min(absDiff_2500ms(:));
    % Find the indices of the closest number
    [~, col_2500ms] = find(absDiff_2500ms == minDiff_2500ms);

    % Baseline correct data
    dataeegepoched = dataeegepoched - mean(dataeegepoched(:,[1:col_500ms,col_3500ms:L],:),2);

    % Convert to absolute values
    absolutevaluesblink = abs(dataeegepoched);

    % Calculate
    BlinkAmplitudeRatioAllEpochs = NaN(NCHANEEG,NTRL);
    for i = 1:NTRL
        BlinkAmplitudeRatioAllEpochs(:,i) = mean(absolutevaluesblink(:,col_1500ms:col_2500ms,i),2) ./ mean(absolutevaluesblink(:,[1:col_500ms,col_3500ms:L],i),2);
    end
    BlinkAmplitudeRatio = mean(BlinkAmplitudeRatioAllEpochs,2);

    th = nexttile;
    % mask = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);
    chanlocs = readlocs('biosemi128_eeglab.ced'); myCmap = brewermap(128,'BuPu'); % BrBG
    topoplot(BlinkAmplitudeRatio,chanlocs,'headrad',0.5,'colormap',myCmap,'whitebk','on','electrodes','off','style','map','shading','interp'); % ,'emarker2',{find(mask),'d','k',10,1}

    maxBlinkRatio = max(BlinkAmplitudeRatio);
    maxBlinkRatio = max(maxBlinkRatio,1.5);
    clim(th,[1, maxBlinkRatio]); colorbar;

    % Save
    plotX=15; plotY=10;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_leftovers']),'-dtiff','-r300');
    close(fh);

else
    warning('Too little data for a robust estimate...');
    BlinkAmplitudeRatio = NaN;
end

% Log
fprintf('Total amount of leftover eye blink artifact: %1.1f%%\n', (mean(BlinkAmplitudeRatio)-1)*100);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Leftovers: eye blink artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of leftover eye blink artifact: %1.1f%%\n', (mean(BlinkAmplitudeRatio)-1)*100);

EEG.ALSUTRECHT.leftovers.blinks = BlinkAmplitudeRatio;

end