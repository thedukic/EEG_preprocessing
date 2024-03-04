function [EEG, badElectrodes, noiseMask] = detect_channelemg(EEG,cfgbch)
%
% Test whether frequency slope indicating high freq oscillations (EMG) is present
% Based on: RELAX_excluding_channels_and_epoching / RELAX_excluding_extreme_values / RELAX_muscle
% SDukic, March 2023
%
fprintf('\nMWF (EMG) muscle artifacts...\n');

% Select only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
dataEEG = EEG.data(chaneeg,:);

% Epoch into 1s (or maybe better into 1s with 0.5 overlap, but OK)
L = EEG.srate;
N = floor(size(dataEEG,2)/L);
dataEEG = reshape(dataEEG(:,1:N*L),sum(chaneeg),L,N);

% Compute (log) power spectra (7-75 Hz)
[NCHN,NPTS,NTRL]= size(dataEEG);
NWIN = NPTS;
psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);
for i = 1:NTRL
    [psdspectra(:,:,i),freq] = pwelch(dataEEG(:,:,i)',NWIN,0,NWIN,EEG.srate);
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
badchn = sum(slopesEpochsxChannels>cfgbch.muscleSlopeThreshold,2);
badchn = badchn./NTRL;
badElectrodes = {EEG.chanlocs(find(badchn>cfgbch.slopeTime)).labels};

% =========================================================================
% OLD CODE
% Make a noise mask using very bad epochs
% badtrl_all = find(any(SlopesEpochsxChannels>cfgbch.muscleSlopeThreshold,1));
% badtrl_all = unique([badtrl_all badtrl_all-1 badtrl_all+1]);
% badtrl_all(badtrl_all==0 | badtrl_all>NTRL) = [];
%
% BadEpochs = sum(SlopesEpochsxChannels>cfgbch.muscleSlopeThreshold,1);
% K = ceil(max(6,median(BadEpochs(BadEpochs>0))));
% % K = max(8,mode(BadEpochs(BadEpochs>0)));
% badtrl_all = find(BadEpochs>K);
%
% if length(badtrl_all)>0.5*NTRL
%     badtrl_all = find(BadEpochs>K+1);
% end
% if length(badtrl_all)>0.5*NTRL
%     badtrl_all = find(BadEpochs>K+2);
% end
%
% % 1: [0 - 1]-0
% % 2: [1 - 2]-0.5
% % 3: [2 - 3]-1
% % 4: [3 - 4]-1.5
% % i: [i-1 i]-(i-1)*0.5
%
% % figure; plot(freq,squeeze(psdspectra(:,badtrl_all(4),:)));
% % figure; plot(freq,squeeze(psdspectra(:,badtrl_all(4),:)));

% =========================================================================
% The following replaces all values that aren't above the muscle
% threshold with NaN, then sums the values that are above the threshold
% for each epoch to allow identification of the worst epochs:

% Exclude extreme epochs
assert(size(slopesEpochsxChannels,2)==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2));
slopesEpochsxChannels(:,EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2) = NaN;

sortingOutWorstMuscleEpochs = slopesEpochsxChannels;
sortingOutWorstMuscleEpochs(sortingOutWorstMuscleEpochs < cfgbch.muscleSlopeThreshold) = NaN;

% Shift the baseline of the values to the cfgbch.muscleSlopeThreshold, so that all
% muscle artifacts can be ranked in severity of EMG starting from 0 (least
% severe) and moving more positive as more severe:
sortingOutWorstMuscleEpochs = sortingOutWorstMuscleEpochs-cfgbch.muscleSlopeThreshold;

% Sum muscle slopes across all channels that show slopes above the
% threshold. This gives an indication of how badly each epoch is
% affected by muscle activity, with more affected electrodes within
% the epoch providing higher values:
sortingOutWorstMuscleEpochs = sum(sortingOutWorstMuscleEpochs,1,'omitnan');

% Work out proportion of data marked as showing muscle activity
assert(NTRL==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2));
templateMarkedForMuscleArtifacts = zeros(1,NTRL);
templateMarkedForMuscleArtifacts(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2) = NaN;

% Threshold = 0, because all slope values have had the threshold subtracted from them (so the threshold is now 0)
muscleSlopeThresholdAfterAdjustment = 0;
templateMarkedForMuscleArtifacts(sortingOutWorstMuscleEpochs>muscleSlopeThresholdAfterAdjustment) = 1;
ProportionOfDataShowingMuscleActivityTotal = mean(templateMarkedForMuscleArtifacts,"omitnan");

fprintf('Total amount of muscle artifact: %1.2f\n', ProportionOfDataShowingMuscleActivityTotal);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'MWF (EMG) muscle artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle log(7-75Hz) slope threshold: %1.2f\n',cfgbch.muscleSlopeThreshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of muscle artifact: %1.2f\n', ProportionOfDataShowingMuscleActivityTotal);

% Checks the proportion of muscle artifact periods selected, and
% reduces the threshold for marking a period if the above section marked more than the allowed threshold
% => 50% of the data
MaxProportionOfDataCanBeMarkedAsMuscle = 0.5;
if ProportionOfDataShowingMuscleActivityTotal > MaxProportionOfDataCanBeMarkedAsMuscle
    fprintf('That is too much for MWF. Limiting it to the worst 50%% ...\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'That is too much for MWF. Limiting it to worst 50%%.\n');

    templateMarkedForMuscleArtifacts = zeros(1,NTRL);
    templateMarkedForMuscleArtifacts(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2) = NaN;
    muscleSlopeThresholdAfterAdjustment = prctile(sortingOutWorstMuscleEpochs,100-(MaxProportionOfDataCanBeMarkedAsMuscle*100));
    templateMarkedForMuscleArtifacts(sortingOutWorstMuscleEpochs>muscleSlopeThresholdAfterAdjustment) = 1;
end

% Find those trials
badtrl_all = templateMarkedForMuscleArtifacts;
badtrl_all(isnan(badtrl_all))=0;
badtrl_all = find(badtrl_all); % otherwise NaNs are "found" as well

%% Multi-channel Wiener Filter
% Noise mask:
% NaN - ignored samples (== very bad data)
% 0   - good samples
% 1   - bad samples (== bad data that will be corrected)

% Noise mask
noiseMask = zeros(1,EEG.pnts);
if ~isempty(badtrl_all)
    badEpoch1 = [badtrl_all-1; badtrl_all]'; % in [s]
    badEpoch2 = badEpoch1*EEG.srate;         % in [samples]
    badEpoch2(:,1) = badEpoch2(:,1)+1;

    for i = 1:length(badtrl_all)
        noiseMask(badEpoch2(i,1):badEpoch2(i,2)) = 1;
    end
end

assert(length(noiseMask)==size(EEG.data,2));
assert(length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1)==length(noiseMask));
noiseMask(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = NaN;

% Log info
EEG.ALSUTRECHT.MWF.EMG.badElectrodes          = badElectrodes;
EEG.ALSUTRECHT.MWF.EMG.noiseMask              = noiseMask;
EEG.ALSUTRECHT.MWF.EMG.proportionMarkedForMWF = mean(noiseMask,'omitnan');
EEG.ALSUTRECHT.MWF.EMG.ProportionOfDataShowingMuscleActivityTotal = ProportionOfDataShowingMuscleActivityTotal;

fprintf('Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.EMG.proportionMarkedForMWF);
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.EMG.proportionMarkedForMWF);

if EEG.ALSUTRECHT.MWF.EMG.proportionMarkedForMWF>0.05
    [cleanEEG, d, W, SER, ARR] = mwf_process(dataEEG(:,:),noiseMask,8);

    % % Check
    % EEG0 = EEG;
    % EEG0.data(1:NCHNEEG,:) = cleanEEG;
    % vis_artifacts(EEG0,EEG);

    EEG.data(chaneeg,:) = cleanEEG;

    % EEG.ALSUTRECHT.MWF.EMG.estimatedArtifactInEachChannel = d;
    % EEG.ALSUTRECHT.MWF.EMG.matrixUsedToEstimateArtifacts  = W;
    EEG.ALSUTRECHT.MWF.EMG.status                 = true;
    EEG.ALSUTRECHT.MWF.EMG.signalToErrorRatio     = SER;
    EEG.ALSUTRECHT.MWF.EMG.artifactToResidueRatio = ARR;

    fprintf(EEG.ALSUTRECHT.subject.fid,'Signal to error ratio:     %1.2f\n',SER);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Artifact to residue ratio: %1.2f\n',ARR);
else
    EEG.ALSUTRECHT.MWF.EMG.status                 = false;
    EEG.ALSUTRECHT.MWF.EMG.signalToErrorRatio     = NaN;
    EEG.ALSUTRECHT.MWF.EMG.artifactToResidueRatio = NaN;

    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF will not be done. Too little data.\n');
    fprintf('MWF will not be done. Too little data.\n');
end

end