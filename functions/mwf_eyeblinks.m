function [EEG, noiseMask] = mwf_eyeblinks(EEG)
%
% Test whether only eye blinks are present using bipolar VEOG
% SDukic, March 2023
%
fprintf('\nMWF (VEOG) eye blinks...\n');

% Select EEG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
dataeeg  = EEG.data(chaneeg,:);

% Detect eye blinks
[dataeog, locs, treshold] = detect_veog(EEG);
NEOG = length(locs);

% figure; hold on;
% plot(times(1:100*256),dataeog(1:100*256));
% plot(locs(1:11),qrspeaks(1:11),'ro');

EOGfocus  = 250; % ms
mspersamp = 1000/EEG.srate;
EOGfocussamples = round(EOGfocus/mspersamp);
fprintf('VEOG duration is +-%dms wrt the detected peaks.\n',EOGfocus);

badEpoch2 = locs*EEG.srate;
badEpoch2 = [badEpoch2-EOGfocussamples; badEpoch2+EOGfocussamples]';
N = length(badEpoch2(1,1):badEpoch2(1,2));
% TEOG = (N-1)/2/EEG.srate*1000;

badEpoch2(badEpoch2<1) = 1;
badEpoch2(badEpoch2>EEG.pnts) = EEG.pnts;

% Calculate the absolute difference between the EEG.times and the timepoints we need:
T = (0:N-1)./EEG.srate;
T = T-mean(T);

absDiff_150ms = abs(T - (-0.1));
minDiff_150ms = min(absDiff_150ms(:));
[~, col_m150ms] = find(absDiff_150ms == minDiff_150ms);
absDiff_150ms = abs(T - 0.1);
minDiff_150ms = min(absDiff_150ms(:));
[~, col_p150ms] = find(absDiff_150ms == minDiff_150ms);

EOG = NaN(NEOG,N);
for i = 1:NEOG
    if length(badEpoch2(i,1):badEpoch2(i,2))==N
        EOG(i,:) = dataeog(badEpoch2(i,1):badEpoch2(i,2));
        EOG(i,:) = EOG(i,:) - mean(EOG(i,[1:col_m150ms, col_p150ms:end]));
        % EOG(i,:) = EOG(i,:) - mean(EOG(i,1:col_m150ms));
    end
end
EOG(isnan(EOG(:,1)),:) = [];
mEOG = mean(EOG,1);

fh = figure; hold on;
plot(T,EOG','LineWidth',1.2);
set(gca, 'ColorOrder',brewermap(size(EOG,1),'BuGn'));
plot(T,mEOG,'Color',[0.8 0.1 0.1],'LineWidth',3);
title(['N = ' num2str(NEOG)]);
xlabel('Time (s)'); ylabel('VEOG amplitude');

% Save
plotX=35; plotY=20;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_VEOG']),'-dtiff','-r400');
close(fh);

%% Multi-channel Wiener filter
% Noise mask:
% NaN - ignored samples (== very bad data)
% 0   - good samples
% 1   - bad samples (== bad data that will be corrected)

% % Apporach 1: Noise mask
% noiseMask = zeros(1,EEG.pnts);
% if ~isempty(badtrl_all)
%     badEpoch1 = [badtrl_all-1; badtrl_all]'; % in [s]
%     badEpoch2 = badEpoch1*EEG.srate;         % in [samples]
%     badEpoch2(:,1) = badEpoch2(:,1)+1;
%
%     for i = 1:length(badtrl_all)
%         noiseMask(badEpoch2(i,1):badEpoch2(i,2)) = 1;
%     end
% end

% Apporach 2: Noise mask
noiseMask = zeros(1,EEG.pnts);
if ~isempty(badEpoch2)
    for i = 1:size(badEpoch2,1)
        noiseMask(badEpoch2(i,1):badEpoch2(i,2)) = 1;
    end
end

% Mask very bad periods of data
assert(length(noiseMask)==size(EEG.data,2));
assert(length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1)==length(noiseMask));
noiseMask(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = NaN;

% Log info
EEG.ALSUTRECHT.MWF.R2.badElectrodes          = NaN;
EEG.ALSUTRECHT.MWF.R2.noiseMask              = noiseMask;
EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF = mean(noiseMask,'omitnan');

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'MWF (VEOG) eye blinks\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Eye blink amplitude threshold: %1.2f\n',treshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of eye blinks detected: %d\n',NEOG);
fprintf(EEG.ALSUTRECHT.subject.fid,'Eye blink duration for MWF: %1.2f\n',2*EOGfocus);
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF);

fprintf('Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF);

if EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF>0.05
    [cleanEEG, d, W, SER, ARR] = mwf_process(dataeeg,noiseMask,8);

    % % Check
    % EEG0 = EEG;
    % EEG0.data(chaneeg,:) = cleanEEG;
    % vis_artifacts(EEG0,EEG);

    EEG.data(chaneeg,:) = cleanEEG;

    % EEG.ALSUTRECHT.MWF.R2.estimatedArtifactInEachChannel = d;
    % EEG.ALSUTRECHT.MWF.R2.matrixUsedToEstimateArtifacts  = W;
    EEG.ALSUTRECHT.MWF.R2.status                 = 1;
    EEG.ALSUTRECHT.MWF.R2.signalToErrorRatio     = SER;
    EEG.ALSUTRECHT.MWF.R2.artifactToResidueRatio = ARR;

    fprintf(EEG.ALSUTRECHT.subject.fid,'Signal to error ratio:     %1.2f\n',SER);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Artifact to residue ratio: %1.2f\n',ARR);
else
    EEG.ALSUTRECHT.MWF.R2.status                 = 0;
    EEG.ALSUTRECHT.MWF.R2.signalToErrorRatio     = NaN;
    EEG.ALSUTRECHT.MWF.R2.artifactToResidueRatio = NaN;

    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF will not be done. Too little data.\n');
    fprintf('MWF will not be done. Too little data.\n');
end

end