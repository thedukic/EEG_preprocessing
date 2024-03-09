function [EEG, badElectrodes, noiseMask] = mwf_channeldrifts(EEG)
%
% Test whether only drift is present for each channel in each epoch based on frequency slope:
% SDukic, March 2023
% Based on: RELAX_excluding_channels_and_epoching / RELAX_excluding_extreme_values
% Performance improvement
%
fprintf('\nMWF (HEOG) horizontal eye movements...\n');

% Select only EEG + HEOG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chanheog = strcmp({EEG.chanlocs.labels},'HEOG');
dataEEG0 = EEG.data(chaneeg,:);
dataheog = EEG.data(chanheog,:);

% =========================================================================
% HEOG
% Temporarily filter  HEOG, bandpass 1-25 Hz
[bl, al] = butter(4,20/(EEG.srate/2),'low');
assert(isstable(bl,al));
dataheog = filtfilt(bl,al,dataheog);

% Demean and take absolute value so that L/R HEOG can be both detected (?)
dataheog = dataheog - trimmean(dataheog,10);
dataheogabs = abs(dataheog);

% HEOG treshold
EOGIQR = iqr(dataheogabs);
EOG75P = prctile(dataheogabs,75);
tresholdh = EOG75P + 1.5*EOGIQR;
eyeBlinksEpochs = dataheogabs>=tresholdh;

% VEOG treshold
[dataveog, ~, tresholdv] = detect_veog(EEG);

noiseMask1 = zeros(1,EEG.pnts);
if any(eyeBlinksEpochs)
    jump = find(diff([false, eyeBlinksEpochs, false])~=0);
    durall = jump(2:2:end)-jump(1:2:end);

    % Minimum of 25 ms of HEOG
    mspersamp = 1000/EEG.srate;
    mindrftdur = round(50/mspersamp);

    % Horizontal eye movements should last at least this long [ms]
    eyeBlinksEpochs = find(durall>mindrftdur);
    NHEOG = length(eyeBlinksEpochs);

    if NHEOG>0
        jumpStart = jump(1:2:end);
        jumpStop  = jump(2:2:end);

        EOGfocus = 300; % ms
        EOGfocussamples = round(EOGfocus/mspersamp);

        jumpStart = jumpStart-EOGfocussamples;
        jumpStart(jumpStart<1) = 1;
        jumpStop = jumpStop+EOGfocussamples;
        jumpStop(jumpStop>EEG.pnts) = EEG.pnts;

        actuallyBlinks = false(NHEOG,1);
        for i = 1:NHEOG
            actuallyBlinks(i) = any(dataveog(jumpStart(eyeBlinksEpochs(i)):jumpStop(eyeBlinksEpochs(i)))>=tresholdv);
        end

        eyeBlinksEpochs(actuallyBlinks) = [];
        % jumpStart(actuallyBlinks)       = [];
        % jumpStop(actuallyBlinks)        = [];
        NHEOG = length(eyeBlinksEpochs);

        % % Check: A lot of times HEOG captures blinks too !!!
        % EEGTMP = EEG;
        % mask = false(size(EEGTMP.times));
        % for i = 1:NHEOG
        %     mask(jumpStart(eyeBlinksEpochs(i)):jumpStop(eyeBlinksEpochs(i))) = true;
        % end
        % EEGTMP.data(:,~mask) = 0;
        % vis_artifacts(EEG,EEGTMP);

        fh = figure; hold on;
        durheog = NaN(NHEOG,1);
        for i = 1:NHEOG
            y = dataheog(jumpStart(eyeBlinksEpochs(i)):jumpStop(eyeBlinksEpochs(i)));
            durheog(i) = length(y);
            T = linspace(-0.5,0.5,durheog(i));

            absDiff_150ms = abs(T - (-0.4));
            minDiff_150ms = min(absDiff_150ms(:));
            [~, col_m150ms] = find(absDiff_150ms == minDiff_150ms);
            absDiff_150ms = abs(T - 0.4);
            minDiff_150ms = min(absDiff_150ms(:));
            [~, col_p150ms] = find(absDiff_150ms == minDiff_150ms);
            y = y - mean(y([1:col_m150ms, col_p150ms:end]));
            % y = y - mean(y(1:col_m150ms));

            plot(T,y,'LineWidth',1.2);
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
        fprintf('HEOG average duration is %2.0fms.\n',TEOG);

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
EEG.ALSUTRECHT.MWF.R3.badElectrodes          = badElectrodes;
EEG.ALSUTRECHT.MWF.R3.noiseMask              = noiseMask;
EEG.ALSUTRECHT.MWF.R3.proportionMarkedForMWF = mean(noiseMask,'omitnan');

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'MWF (HEOG) horizontal eye movements\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Horizontal eye movement amplitude threshold: %1.2f\n',tresholdh);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of horizontal eye movements detected: %d\n',NHEOG);
fprintf(EEG.ALSUTRECHT.subject.fid,'Average horizontal eye movement duration detected: %1.2f\n',TEOG);
% fprintf(EEG.ALSUTRECHT.subject.fid,'Number of slow drifts detected: %d\n',length(driftEEGEpochs));
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R3.proportionMarkedForMWF);

fprintf('Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R3.proportionMarkedForMWF);

if EEG.ALSUTRECHT.MWF.R3.proportionMarkedForMWF>0.05
    [cleanEEG, d, W, SER, ARR] = mwf_process(dataEEG0,noiseMask,8);

    % EEG0 = EEG;
    % EEG0.data(chaneeg,:) = cleanEEG;
    % vis_artifacts(EEG0,EEG);

    EEG.data(chaneeg,:) = cleanEEG;

    % EEG.ALSUTRECHT.MWF.R3.estimatedArtifactInEachChannel = d;
    % EEG.ALSUTRECHT.MWF.R3.matrixUsedToEstimateArtifacts  = W;
    EEG.ALSUTRECHT.MWF.R3.status                 = 1;
    EEG.ALSUTRECHT.MWF.R3.signalToErrorRatio     = SER;
    EEG.ALSUTRECHT.MWF.R3.artifactToResidueRatio = ARR;

    fprintf(EEG.ALSUTRECHT.subject.fid,'Signal to error ratio:     %1.2f\n',SER);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Artifact to residue ratio: %1.2f\n',ARR);
else
    EEG.ALSUTRECHT.MWF.R3.status                 = 0;
    EEG.ALSUTRECHT.MWF.R3.signalToErrorRatio     = NaN;
    EEG.ALSUTRECHT.MWF.R3.artifactToResidueRatio = NaN;

    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF will not be done. Too little data.\n');
    fprintf('MWF will not be done. Too little data.\n');
end

end