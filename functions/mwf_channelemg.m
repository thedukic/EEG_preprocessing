function [EEG, badElectrodes, noiseMask] = mwf_channelemg(EEG,cfgbch)
%
% Test whether frequency slope indicating high freq oscillations (EMG) is present
% Based on: RELAX_excluding_channels_and_epoching / RELAX_excluding_extreme_values / RELAX_muscle
%
% MuscleSlopeThreshold: (-0.31 = data from paralysed participants showed no
% independent components with a slope value more positive than this (so
% excluding slopes above this threshold means only excluding data that we
% know must be influenced by EMG). Using -0.31 as the threshold means
% possibly leaving low level EMG data in, and only eliminating the data we
% know is definitely EMG)
% (-0.59 is where the histograms between paralysed ICs and EMG ICs cross,
% so above this value contains a very small amount of the brain data,
% and over 50% of the EMG data. Above this point, data is more likely to be
% EMG than brain)
% (-0.72 is the maximum of the histogram of the paralysed IC data, so
% excluding more positive values than this will exclude most of the EMG
% data, but also some brain data).
%
% Fitzgibbon, S. P., DeLosAngeles, D., Lewis, T. W., Powers, D. M. W., Grummett, T. S., Whitham, E. M., ... & Pope, K. J. (2016). Automatic determination of EMG-contaminated components and validation of independent component analysis using EEG during pharmacologic paralysis. Clinical Neurophysiology, 127(3), 1781-1793.
%
% SDukic, July 2024
%

fprintf('\nMWF (EMG) muscle artifacts...\n');

% =========================================================================
% Estimate log-log power spectra
[slopesChannelsxEpochs, other] = dected_emg(EEG);

% Account for very bad epochs affected across all channels
if any(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2)
    assert(size(slopesChannelsxEpochs,2)==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2));
    slopesChannelsxEpochs(:,EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2) = NaN;
end

% Detect noisy channels that are left in the data
badElectrodes = mean(slopesChannelsxEpochs>cfgbch.muscleSlopeThreshold,2,'omitnan');
badElectrodes = {EEG.chanlocs(badElectrodes>cfgbch.muscleSlopeTime).labels};

% =========================================================================
% The following replaces all values that aren't above the muscle
% threshold with NaN, then sums the values that are above the threshold
% for each epoch to allow identification of the worst epochs:

sortingOutWorstMuscleEpochs = slopesChannelsxEpochs;
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
NTRL = size(slopesChannelsxEpochs,2);
templateMarkedForMuscleArtifacts = zeros(1,NTRL);

% Threshold = 0, because all slope values have had the threshold subtracted from them (so the threshold is now 0)
muscleSlopeThresholdAfterAdjustment = 0;
templateMarkedForMuscleArtifacts(sortingOutWorstMuscleEpochs>muscleSlopeThresholdAfterAdjustment) = 1;

% Mark extremly bad epochs as they may now appear as EMG-free!
assert(NTRL==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2));
templateMarkedForMuscleArtifacts(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2) = NaN;

proportionOfDataShowingMuscleActivityTotal = mean(templateMarkedForMuscleArtifacts,"omitnan");

fprintf('Total amount of muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'MWF (EMG) muscle artifacts\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle log(7-75Hz) slope threshold: %1.2f\n',cfgbch.muscleSlopeThreshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Total amount of muscle artifact: %1.2f\n', proportionOfDataShowingMuscleActivityTotal);

% Checks the proportion of muscle artifact periods selected, and
% reduces the threshold for marking a period if the above section marked more than the allowed threshold
% i.e. >50% of the data
maxProportionOfDataCanBeMarkedAsMuscle = 0.5;
if proportionOfDataShowingMuscleActivityTotal > maxProportionOfDataCanBeMarkedAsMuscle
    fprintf('That is too much for MWF. Limiting it to the worst 50%% ...\n');
    fprintf(EEG.ALSUTRECHT.subject.fid,'That is too much for MWF. Limiting it to worst 50%%.\n');

    templateMarkedForMuscleArtifacts = zeros(1,NTRL);
    % 1. Original
    % muscleSlopeThresholdAfterAdjustment = prctile(sortingOutWorstMuscleEpochs,100-(MaxProportionOfDataCanBeMarkedAsMuscle*100));
    % 2. Use only those epochs that are not extremly bad, as they can obscure the following estimation
    extremelyBadEpochsExcluded = sortingOutWorstMuscleEpochs(~EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2);
    muscleSlopeThresholdAfterAdjustment = prctile(extremelyBadEpochsExcluded,100-(maxProportionOfDataCanBeMarkedAsMuscle*100));

    templateMarkedForMuscleArtifacts(sortingOutWorstMuscleEpochs>=muscleSlopeThresholdAfterAdjustment) = 1;
    templateMarkedForMuscleArtifacts(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2) = NaN;
end

% Find those trials
badTrialsAll = templateMarkedForMuscleArtifacts;
badTrialsAll = find(badTrialsAll==1);

%% Multi-channel Wiener Filter
% Noise mask:
% NaN - ignored samples (== very bad data)
% 0   - good samples
% 1   - bad samples (== bad data that will be corrected)

% Noise mask
noiseMask = zeros(1,EEG.pnts);
if ~isempty(badTrialsAll)
    badEpochs1 = [badTrialsAll-1; badTrialsAll]'; % in [s]
    badEpochs2 = badEpochs1*EEG.srate;            % in [samples]
    badEpochs2(:,1) = badEpochs2(:,1)+1;

    for i = 1:length(badTrialsAll)
        noiseMask(badEpochs2(i,1):badEpochs2(i,2)) = 1;
    end
end

% Mark extremly bad epochs
assert(length(noiseMask)==size(EEG.data,2));
assert(length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1)==length(noiseMask));
noiseMask(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = NaN;

% Ignore the very last samples
% because we do not know if they are good or not
if other.modulus>0
    noiseMask(other.N*other.L+1:end) = NaN;
end

% % Visual inspection
% EEG1 = EEG;
% EEG1.data(:,noiseMask==0 | isnan(noiseMask)) = 0;
% vis_artifacts(EEG,EEG1);

% Log info
EEG.ALSUTRECHT.MWF.R1.badElectrodes          = badElectrodes;
EEG.ALSUTRECHT.MWF.R1.noiseMask              = noiseMask;
EEG.ALSUTRECHT.MWF.R1.proportionMarkedForMWF = mean(noiseMask,'omitnan');
EEG.ALSUTRECHT.MWF.R1.ProportionOfDataShowingMuscleActivityTotal = proportionOfDataShowingMuscleActivityTotal;

fprintf('Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R1.proportionMarkedForMWF);
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R1.proportionMarkedForMWF);

if EEG.ALSUTRECHT.MWF.R1.proportionMarkedForMWF>0.05
    % Use only EEG channels
    chaneeg = strcmp({EEG.chanlocs.type},'EEG');
    params  = mwf_params('delay',12,'delay_spacing',2);
    [cleanEEG, d, W, SER, ARR] = mwf_process(EEG.data(chaneeg,:),noiseMask,params);

    % Check if there were any problems
    if contains(lastwarn,"eigenvectors")
        warning('The MWF delay is too long?'); SER = Inf; ARR = Inf;
    end
    if isnan(SER) || isnan(ARR)
        warning('MWF did not fail but the MWF quality measures (SER/ARR) are NaN. The bad data might be too short.');
    end

    % % Visual inspection
    % EEG0 = EEG;
    % EEG0.data(chaneeg,:) = cleanEEG;
    % vis_artifacts(EEG0,EEG);

    % Return the clean data
    EEG.data(chaneeg,:) = cleanEEG;

    % Log
    % EEG.ALSUTRECHT.MWF.R1.estimatedArtifactInEachChannel = d;
    % EEG.ALSUTRECHT.MWF.R1.matrixUsedToEstimateArtifacts  = W;
    EEG.ALSUTRECHT.MWF.R1.status                 = 1;
    EEG.ALSUTRECHT.MWF.R1.signalToErrorRatio     = SER;
    EEG.ALSUTRECHT.MWF.R1.artifactToResidueRatio = ARR;

    fprintf('Signal-to-error ratio:     %1.2f\n',SER);
    fprintf('Artifact-to-residue ratio: %1.2f\n',ARR);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Signal-to-error ratio:     %1.2f\n',SER);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Artifact-to-residue ratio: %1.2f\n',ARR);
else
    EEG.ALSUTRECHT.MWF.R1.status                 = 0;
    EEG.ALSUTRECHT.MWF.R1.signalToErrorRatio     = NaN;
    EEG.ALSUTRECHT.MWF.R1.artifactToResidueRatio = NaN;

    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF is not done. Too little data.\n');
    fprintf('MWF is not done. Too little data.\n');
end

end