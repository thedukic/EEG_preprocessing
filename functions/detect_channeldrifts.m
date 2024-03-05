function [EEG, badElectrodes, noiseMask] = detect_channeldrifts(EEG)
%
% Test whether only drift is present for each channel in each epoch based on frequency slope:
% SDukic, March 2023
% Based on: RELAX_excluding_channels_and_epoching / RELAX_excluding_extreme_values
% Performance improvement
%
fprintf('\nMWF (HEOG) horizontal eye movements...\n');

% Select only EEG + HEOG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chaneog  = strcmp({EEG.chanlocs.labels},'HEOG');
dataEEG0 = EEG.data(chaneeg,:);
dataeog  = EEG.data(chaneog,:);

% =========================================================================
% HEOG
% Temporarily filter  HEOG, bandpass 1-25 Hz
[bl, al] = butter(2,5/(EEG.srate/2),'low');
assert(isstable(bl,al));
dataeog = filtfilt(bl,al,dataeog);
dataeog = abs(dataeog);

% See: RELAX_blinks_IQR_method
EOGIQR = iqr(dataeog);
EOG75P = prctile(dataeog,75);
treshold = EOG75P + 2*EOGIQR;
eyeBlinksEpochs = dataeog>=treshold;

noiseMask1 = zeros(1,EEG.pnts);
if any(eyeBlinksEpochs)
    jump = find(diff([false, eyeBlinksEpochs, false])~=0);
    durall = jump(2:2:end)-jump(1:2:end);

    % Minimum of 25 ms of HEOG
    mspersamp = 1000/EEG.srate;
    mindrftdur = round(25/mspersamp);

    % Horizontal eye movements should last at least this long [ms]
    eyeBlinksEpochs = find(durall>mindrftdur);
    NHEOG = length(eyeBlinksEpochs);

    if NHEOG>0
        jumpStart = jump(1:2:end);
        jumpStop  = jump(2:2:end);

        jumpStart = jumpStart-round(25/mspersamp);
        jumpStart(jumpStart<1) = 1;
        jumpStop = jumpStop+round(25/mspersamp);
        jumpStop(jumpStop>EEG.pnts) = EEG.pnts;

        fh = figure; hold on;
        durheog = NaN(NHEOG,1);
        for i = 1:NHEOG
            y = dataeog(jumpStart(eyeBlinksEpochs(i)):jumpStop(eyeBlinksEpochs(i)));
            durheog(i) = length(y);
            t = linspace(0,1,durheog(i));
            plot(t,y,'LineWidth',1.2);
        end
        set(gca,'ColorOrder',brewermap(NHEOG,'BuGn'));
        title(['N = ' num2str(NHEOG)]);
        xlabel('Time a.u.'); ylabel('HEOG amplitude');

        % Save
        plotX=35; plotY=20;
        set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
        set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
        print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_HEOG']),'-dtiff','-r400');
        close(fh);

        TEOG = mean(durheog)/EEG.srate*1000;
        fprintf('HEOG average duration is %2.0fms wrt the detected peaks.\n',TEOG);

        for i = 1:NHEOG
            noiseMask1(jumpStart(eyeBlinksEpochs(i)):jumpStop(eyeBlinksEpochs(i))) = 1;
        end
    end
end

% =========================================================================
% % EEG
% % Temporarily filter EEG, lowpass 5 Hz
% dataEEG1 = dataEEG0;
% [bl, al] = butter(4,5/(EEG.srate/2),'low');
% assert(isstable(bl,al));
% dataEEG1 = filtfilt(bl,al,dataEEG1')';
% dataEEG1 = dataEEG1 - trimmean(dataEEG1,10);
%
% % Epoch EEG into 1s
% L = EEG.srate;
% N = floor(size(dataEEG1,2)/L);
% dataEEG1 = reshape(dataEEG1(:,1:N*L),sum(chaneeg),L,N);
%
% % Calculate the median across all channels, and the MAD, and the
% % threshold based on these from the number of MAD from the median that
% % you have chosen as the threshold:
% medianAmplitude = median(dataEEG1(:,:,:),1);
% MADAmplitude    = mad(dataEEG1(:,:,:),1);
% upperBound      = medianAmplitude+(cfgbch.DriftSeverityThreshold*MADAmplitude);
% lowerBound      = medianAmplitude-(cfgbch.DriftSeverityThreshold*MADAmplitude);
%
% % If voltage exceeds MAD * your drift severity threshold for a total of
% % 100ms then mark the epoch as showing drift in the MWF mask
% mindrftdur = round(100/mspersamp);
%
% NTRL = size(dataEEG1,3);
% driftEEGEpochs = zeros(1,NTRL);
% for i = 1:NTRL
%     for j = 1:size(dataEEG1,1)
%         maxdur = find_maxduration(dataEEG1(j,:,i)>=upperBound(1,:,i));
%         if maxdur>mindrftdur
%             driftEEGEpochs(i) = 1;
%         end
%         maxdur = find_maxduration(dataEEG1(j,:,i)<=lowerBound(1,:,i));
%         if maxdur>mindrftdur
%             driftEEGEpochs(i) = 1;
%         end
%     end
% end
% driftEEGEpochs = find(driftEEGEpochs);

% The above apporach would be probably good if
% slow drifts would be frequent and nonrandom
driftEEGEpochs = [];

% Noise mask
noiseMask2 = zeros(1,EEG.pnts);
if ~isempty(driftEEGEpochs)
    badEpoch1 = [driftEEGEpochs-1; driftEEGEpochs]'; % in [s]
    badEpoch2 = badEpoch1*EEG.srate;                 % in [samples]
    badEpoch2(:,1) = badEpoch2(:,1)+1;

    for i = 1:length(driftEEGEpochs)
        noiseMask2(badEpoch2(i,1):badEpoch2(i,2)) = 1;
    end
end

% % Visual check
% EEG1 = EEG;
% EEG1.data(:,~logical(noiseMask2)) = 0;
% vis_artifacts(EEG,EEG1);

% =========================================================================
% Cobine the two apporaches for detecting slow drifts
noiseMask = noiseMask1+noiseMask2;
badElectrodes = {};

%% Multi-channel Wiener filter
% Noise mask:
% NaN - ignored samples (== very bad data)
% 0   - good samples
% 1   - bad samples (== bad data that will be corrected)

assert(length(noiseMask)==size(EEG.data,2));
assert(length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1)==length(noiseMask));
noiseMask(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = NaN;

% Log info
EEG.ALSUTRECHT.MWF.drift.badElectrodes          = badElectrodes;
EEG.ALSUTRECHT.MWF.drift.noiseMask              = noiseMask;
EEG.ALSUTRECHT.MWF.drift.proportionMarkedForMWF = mean(noiseMask,'omitnan');

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'MWF (HEOG) horizontal eye movements\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Horizontal eye movement amplitude threshold: %1.2f\n',treshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of horizontal eye movements detected: %d\n',NHEOG);
fprintf(EEG.ALSUTRECHT.subject.fid,'Average horizontal eye movement duration detected: %1.2f\n',TEOG);
% fprintf(EEG.ALSUTRECHT.subject.fid,'Number of slow drifts detected: %d\n',length(driftEEGEpochs));
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.drift.proportionMarkedForMWF);

fprintf('Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.drift.proportionMarkedForMWF);

if EEG.ALSUTRECHT.MWF.drift.proportionMarkedForMWF>0.05
    [cleanEEG, d, W, SER, ARR] = mwf_process(dataEEG0,noiseMask,8);

    % EEG0 = EEG;
    % EEG0.data(chaneeg,:) = cleanEEG;
    % vis_artifacts(EEG0,EEG);

    EEG.data(chaneeg,:) = cleanEEG;

    % EEG.ALSUTRECHT.MWF.drift.estimatedArtifactInEachChannel = d;
    % EEG.ALSUTRECHT.MWF.drift.matrixUsedToEstimateArtifacts  = W;
    EEG.ALSUTRECHT.MWF.drift.status                 = 1;
    EEG.ALSUTRECHT.MWF.drift.signalToErrorRatio     = SER;
    EEG.ALSUTRECHT.MWF.drift.artifactToResidueRatio = ARR;

    fprintf(EEG.ALSUTRECHT.subject.fid,'Signal to error ratio:     %1.2f\n',SER);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Artifact to residue ratio: %1.2f\n',ARR);
else
    EEG.ALSUTRECHT.MWF.drift.status                 = 0;
    EEG.ALSUTRECHT.MWF.drift.signalToErrorRatio     = NaN;
    EEG.ALSUTRECHT.MWF.drift.artifactToResidueRatio = NaN;

    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF will not be done. Too little data.\n');
    fprintf('MWF will not be done. Too little data.\n');
end

end