function EEG = report_leftovers(EEG)
%
% Inspired by:
% RELAX_metrics_muscle
% RELAX_metrics_blinks
%
% SDukic, March 2023
%

% =========================================================================
fprintf('\nChecking muscle activity letovers...\n');
muscleSlopeThreshold = - 0.59;
muscleSlopeDuration  = 0.05;

% Select only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
dataeeg = EEG.data(chaneeg,:);

% Epoch into 1s (or maybe better into 1s with 0.5 overlap, but OK)
L = EEG.srate;
N = floor(size(dataeeg,2)/L);
dataeeg = reshape(dataeeg(:,1:N*L),sum(chaneeg),L,N);

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
ProportionOfDataShowingMuscleActivityTotal = mean(templateMarkedForMuscleArtifacts,'omitnan');

% Log
fprintf('Total amount of leftover muscle artifact: %1.2f\n', ProportionOfDataShowingMuscleActivityTotal);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Leftovers: muscle artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle log(7-75Hz) slope threshold: %1.2f\n',muscleSlopeThreshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of leftover muscle artifact: %1.2f\n', ProportionOfDataShowingMuscleActivityTotal);

EEG.ALSUTRECHT.leftovers.muscle = ProportionOfDataShowingMuscleActivityTotal;

% =========================================================================
fprintf('\nChecking eye blink letovers...\n');

% Select only EEG + VEOG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
dataeeg  = EEG.data(chaneeg,:);
dataeog  = EEG.data(chaneog,:);

% VEOG
EOGIQR = iqr(dataeog);
EOG75P = prctile(dataeog,75);
treshold = EOG75P + 3*EOGIQR;

% Detect eye blinks
EOGfocus  = 2000; % ms
mspersamp = 1000/EEG.srate;
EOGfocussamples = round(EOGfocus/mspersamp);
fprintf('VEOG evaluation window is %ds around the detected peaks.\n',2*EOGfocus/1000);

% times = EEG.times/1000; % EEGLAB time is in [ms]
times = (0:prod(size(EEG.data,[2 3]))-1)./EEG.srate;
[qrspeaks,locs] = findpeaks(dataeog,times,'MinPeakHeight',treshold);

% figure; hold on;
% plot(times(1:20*256),dataeog(1:20*256));
% plot(locs(1:6),qrspeaks(1:6),'ro');

NEOG = length(locs);
goodEyeBlinks = false(NEOG,1);
for i = 2:NEOG-1
    d = [locs(i)-locs(i-1), locs(i+1)-locs(i)];
    if all(d>2.25)
        goodEyeBlinks(i) = true;
    end
end
locs = locs(goodEyeBlinks);
NEOG = length(locs);
% fprintf('Number of detected peaks is %d.\n',NEOG);

eyeBlinksEpochs = round(locs*EEG.srate); % THIS SHOULD BE INTEGER BY DEFAULT
eyeBlinksEpochs = [eyeBlinksEpochs-EOGfocussamples; eyeBlinksEpochs+EOGfocussamples]';

% N = length(badEpoch2(1,1):badEpoch2(1,2));
N = 2*EOGfocussamples+1;
T = (0:N-1)./EEG.srate*1000;

eyeBlinksEpochs(eyeBlinksEpochs<1) = 1;
eyeBlinksEpochs(eyeBlinksEpochs>EEG.pnts) = EEG.pnts;

goodEyeBlinks = false(NEOG,1);
for i = 1:NEOG
    if length(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2))==N
        goodEyeBlinks(i) = true;
    end
end
eyeBlinksEpochs = eyeBlinksEpochs(goodEyeBlinks,:);

NEOG = size(eyeBlinksEpochs,1);
fprintf('Number of detected peaks is %d.\n',NEOG);

if NEOG>4
    % % Check
    % EEGTMP = EEG;
    % mask = false(size(EEGTMP.times));
    % for i = 1:NEOG
    %     mask(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2)) = true;
    % end
    % EEGTMP.data(:,~mask) = 0;
    % vis_artifacts(EEG,EEGTMP);

    dataeegepoched = NaN(sum(chaneeg),N,NEOG);
    dataeogepoched = NaN(N,NEOG);
    for i = 1:NEOG
        dataeegepoched(:,:,i) = dataeeg(:,eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
        dataeogepoched(:,i) = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));
    end

    % yrange = 1.05*[min(dataeogepoched(:)), max(dataeogepoched(:))];
    % figure; hold on;
    % plot(T,dataeogepoched);
    % plot([500 500],yrange,'Color',[0 0 0]);
    % plot([3500 3500],yrange,'Color',[0 0 0]);
    % plot([1500 1500],yrange,'Color',[0 0 0]);
    % plot([2500 2500],yrange,'Color',[0 0 0]);
    % axis tight;

    % RELAX code
    % Calculate the absolute difference between the EEG.times and the timepoints we need:
    absDiff_500ms = abs(T - 500);
    % Find the minimum absolute difference
    minDiff_500ms = min(absDiff_500ms(:));
    % Find the indices of the closest number
    [~, col_500ms] = find(absDiff_500ms == minDiff_500ms);
    % 3500ms:
    absDiff_3500ms = abs(T - 3500);
    % Find the minimum absolute difference
    minDiff_3500ms = min(absDiff_3500ms(:));
    % Find the indices of the closest number
    [~, col_3500ms] = find(absDiff_3500ms == minDiff_3500ms);
    % 1500ms:
    absDiff_1500ms = abs(T - 1500);
    % Find the minimum absolute difference
    minDiff_1500ms = min(absDiff_1500ms(:));
    % Find the indices of the closest number
    [~, col_1500ms] = find(absDiff_1500ms == minDiff_1500ms);
    % 2500ms:
    absDiff_2500ms = abs(T - 2500);
    % Find the minimum absolute difference
    minDiff_2500ms = min(absDiff_2500ms(:));
    % Find the indices of the closest number
    [~, col_2500ms] = find(absDiff_2500ms == minDiff_2500ms);

    % Baseline correct data
    dataeegepoched = dataeegepoched - mean(dataeegepoched(:,[1:col_500ms,col_3500ms:N],:),2);

    % Convert to absolute values
    absolutevaluesblink = abs(dataeegepoched);

    % Calculate
    BlinkAmplitudeRatioAllEpochs = NaN(sum(chaneeg),NEOG);
    for i = 1:NEOG
        BlinkAmplitudeRatioAllEpochs(1:size(absolutevaluesblink,1),i) = mean(absolutevaluesblink(:,col_1500ms:col_2500ms,i),2) ./ mean(absolutevaluesblink(:,[1:col_500ms,col_3500ms:N],i),2);
    end
    BlinkAmplitudeRatio = mean(BlinkAmplitudeRatioAllEpochs,2);

    % figure;
    % chanlocs = readlocs('biosemi128_eeglab.ced'); myCmap = brewermap(128,'BrBG');
    % topoplot(BlinkAmplitudeRatio,chanlocs,'headrad','rim','colormap',myCmap,'whitebk','on','electrodes','off','style','map');
    % clim([1, max(BlinkAmplitudeRatio)]); colorbar;

else
    warning('Maybe too little data for a robust estimate...?');
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