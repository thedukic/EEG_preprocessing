function EEG = detect_extremelybadepochs(EEG)

cfg = [];
cfg.extremeAbsoluteVoltageThreshold               = 500;   % microvolts max or min above which will be excluded from cleaning and deleted from data
cfg.extremeSingleChannelKurtosisThreshold         = 10;    % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis
cfg.extremeAllChannelKurtosisThreshold            = 10;    % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis
cfg.extremeImprobableVoltageDistributionThreshold = 10;    % SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data (default: 8)
cfg.extremeDriftSlopeThreshold                    = -4;    % slope of log frequency log power below which to reject as drift without neural activity (default: -4)
cfg.extremeMuscleSlopeThreshold                   = -0.25; % Less stringent = -0.31, Middle Stringency = -0.59, More stringent = -0.72

%%
fprintf('\nDetecting extremly bad epochs...\n');

% Epoch into 1s
L = EEG.srate;
NTRL = floor(size(EEG.data,2)/L);

EEGTMP = pop_select(EEG,'channel',{EEG.chanlocs(strcmp({EEG.chanlocs.type},'EEG')).labels});
EEGTMP = eeg_regepochs(EEGTMP,'recurrence',1,'limits',[0 1-1/EEGTMP.srate]);
EEGTMP = eeg_checkset(EEGTMP);
assert(EEGTMP.trials==NTRL);

%% Absolute threshold to identify absolute amplitude extreme values:
maxInEpoch = squeeze(max(EEGTMP.data,[],2));
minInEpoch = squeeze(min(EEGTMP.data,[],2));

% If maximum voltage within the epoch exceeds the maximum threshold, mark as artifact:
EEG.ALSUTRECHT.extremeNoise.absoluteAmplitudeExceededThreshold = zeros(1,NTRL);
EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections   = zeros(1,NTRL);

for i = 1:NTRL
    if sum(maxInEpoch(:,i)>cfg.extremeAbsoluteVoltageThreshold,1)>0
        EEG.ALSUTRECHT.extremeNoise.absoluteAmplitudeExceededThreshold(i) = 1;
        EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections(i)   = 1;
    end
end

% If minimum voltage within the epoch exceeds the minimum threshold, mark as artifact:
for i = 1:NTRL
    if sum(minInEpoch(:,i)<-cfg.extremeAbsoluteVoltageThreshold,1)>0
        EEG.ALSUTRECHT.extremeNoise.absoluteAmplitudeExceededThreshold(i) = 1;
        EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections(i)   = 1;
    end
end

%% Kurtosis to identify epochs with abnormally peaky or flat distributions of voltage values and mark those epochs in the extreme outlier template:
EEGTMP = pop_rejkurt(EEGTMP,1,1:EEGTMP.nbchan,cfg.extremeSingleChannelKurtosisThreshold,cfg.extremeAllChannelKurtosisThreshold,0,0);

% Obtain the highest single electrode kurtosis value within each epoch across all electrodes:
% If the highest value in the epoch is below the threshold, mark as 0, if above, mark as 1:
EEG.ALSUTRECHT.extremeNoise.marksForSingleChannelKurt = max(EEGTMP.stats.kurtE);
EEG.ALSUTRECHT.extremeNoise.marksForSingleChannelKurt(EEG.ALSUTRECHT.extremeNoise.marksForSingleChannelKurt <= cfg.extremeSingleChannelKurtosisThreshold) = 0;
EEG.ALSUTRECHT.extremeNoise.marksForSingleChannelKurt(EEG.ALSUTRECHT.extremeNoise.marksForSingleChannelKurt > cfg.extremeSingleChannelKurtosisThreshold)  = 1;

% For all channel kurtosis, if the highest value in the epoch is below the threshold, mark as 0, if above, mark as 1:
EEG.ALSUTRECHT.extremeNoise.marksForAllChannelKurt = (EEGTMP.stats.kurt)';
EEG.ALSUTRECHT.extremeNoise.marksForAllChannelKurt(EEG.ALSUTRECHT.extremeNoise.marksForAllChannelKurt <= cfg.extremeAllChannelKurtosisThreshold) = 0;
EEG.ALSUTRECHT.extremeNoise.marksForAllChannelKurt(EEG.ALSUTRECHT.extremeNoise.marksForAllChannelKurt > cfg.extremeAllChannelKurtosisThreshold)  = 1;

% Add markings from single and all channel kurtosis together in single matrix:
EEG.ALSUTRECHT.extremeNoise.marksForBothSingleAndAllChannelKurt = EEG.ALSUTRECHT.extremeNoise.marksForSingleChannelKurt;
for i = 1:NTRL
    if EEG.ALSUTRECHT.extremeNoise.marksForAllChannelKurt(1,i) == 1
        EEG.ALSUTRECHT.extremeNoise.marksForBothSingleAndAllChannelKurt(1,i) = 1;
        EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections(1,i)    = 1;
    end
end

%% Marking bad epochs for improbable distributions of voltages within the epoch:
% Entering into the mask any epochs that show any channel with more than the threshold SD deviation in absolute voltage:
EEGTMP = pop_jointprob(EEGTMP,1,1:EEGTMP.nbchan,cfg.extremeImprobableVoltageDistributionThreshold, 1000, 0, 0, 0);

% Obtain the highest single electrode improbable voltage distribution values within each epoch across all electrodes:
EEG.ALSUTRECHT.extremeNoise.improbableVoltageDistributionExceededThreshold = max(EEGTMP.stats.jpE);

% If the highest value in the epoch is below the threshold, mark as 0, if above, mark as 1:
EEG.ALSUTRECHT.extremeNoise.improbableVoltageDistributionExceededThreshold(EEG.ALSUTRECHT.extremeNoise.improbableVoltageDistributionExceededThreshold <= cfg.extremeImprobableVoltageDistributionThreshold) = 0;
EEG.ALSUTRECHT.extremeNoise.improbableVoltageDistributionExceededThreshold(EEG.ALSUTRECHT.extremeNoise.improbableVoltageDistributionExceededThreshold > cfg.extremeImprobableVoltageDistributionThreshold)  = 1;
for i = 1:NTRL
    if EEG.ALSUTRECHT.extremeNoise.improbableVoltageDistributionExceededThreshold(1,i)==1
        EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections(1,i) = 1;
    end
end

%% Test drifts based on frequency slopes
% Compute (log) power spectra (1-75 Hz)
[NCHN,NPTS,NTRL]= size(EEGTMP.data);
psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);
for i = 1:NTRL
    [psdspectra(:,:,i),freq] = pwelch(EEGTMP.data(:,:,i)',NPTS,0,NPTS,EEGTMP.srate);
end
psdspectra = permute(psdspectra,[1 3 2]);

foi = [1 75];
frqmsk = freq>=foi(1) & freq<=foi(2);
logpow = log10(psdspectra(frqmsk,:,:));
logfoi = log10(freq(frqmsk));

% Fit linear regression to log-log data, and store the slope
driftRatioEpochsxChannels = NaN(NCHN,NTRL);
for i = 1:NCHN
    for j = 1:NTRL
        p = polyfit(logfoi,logpow(:,j,i),1);
        driftRatioEpochsxChannels(i,j) = p(1);
    end
end

% If log-frequency log-power slope exceeds threshold
EEG.ALSUTRECHT.extremeNoise.driftSlopeMaskEpochs = zeros(1,NTRL);
for i = 1:NTRL
    if any(driftRatioEpochsxChannels(:,i)<cfg.extremeDriftSlopeThreshold,1)
        EEG.ALSUTRECHT.extremeNoise.driftSlopeMaskEpochs(1,i)             = 1;
        EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections(1,i) = 1;
    end
end

%% Test EMG based on frequency slopes
% Compute (log) power spectra (7-75 Hz)
foi = [7 75];
frqmsk = freq>=foi(1) & freq<=foi(2);
logpow = log10(psdspectra(frqmsk,:,:));
logfoi = log10(freq(frqmsk));

% Fit linear regression to log-log data, and store the slope
driftRatioEpochsxChannels = NaN(NCHN,NTRL);
for i = 1:NCHN
    for j = 1:NTRL
        p = polyfit(logfoi,logpow(:,j,i),1);
        driftRatioEpochsxChannels(i,j) = p(1);
    end
end

% If log-frequency log-power slope exceeds threshold
EEG.ALSUTRECHT.extremeNoise.muscleSlopeMaskEpochs = zeros(1,NTRL);
for i = 1:NTRL
    if mean(driftRatioEpochsxChannels(:,i)>cfg.extremeMuscleSlopeThreshold) > 0.5
        EEG.ALSUTRECHT.extremeNoise.muscleSlopeMaskEpochs(1,i)            = 1;
        EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections(1,i) = 1;
    end
end

%% Combine the pop masks into the full noise mask:
extremeEpochs = find(EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections);
% extremeEpochs = find(EEG.ALSUTRECHT.extremeNoise.muscleSlopeMaskEpochs);

extremeNoiseMask = false(1,EEG.pnts);
if ~isempty(extremeEpochs)
    badepoch1 = [extremeEpochs-1; extremeEpochs]'; % in [s]
    badepoch2 = badepoch1*EEGTMP.srate;            % in [samples]
    badepoch2(:,1) = badepoch2(:,1)+1;

    for i = 1:length(extremeEpochs)
        extremeNoiseMask(badepoch2(i,1):badepoch2(i,2)) = true;
    end
else
    badepoch2 = [];
end

% Log
EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier = mean(extremeNoiseMask);
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1                 = extremeNoiseMask;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2                 = logical(EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections);
% EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3                 = extremeEpochs;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3                 = badepoch2;

fprintf('Percentage of extremly bad EEG: %1.2f\n', EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);
fprintf('These data will be excluded from MWF and ICA.\n');

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Extremly bad epochs\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'These data will be excluded from MWF and ICA.\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected: %1.2f\n', EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);

% % Visual check
% EEG1 = EEG;
% EEG1.data(:,~logical(extremeNoiseMask)) = 0;
% vis_artifacts(EEG,EEG1);

end