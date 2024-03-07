function [EEG, noiseMask] = detect_eyeblinks(EEG)
%
% Test whether only eye blinks are present using bipolar VEOG
% SDukic, March 2023
%
fprintf('\nMWF (VEOG) eye blinks...\n');

% Select only EEG + VEOG
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
dataeeg  = EEG.data(chaneeg,:);
dataeog  = EEG.data(chaneog,:);

% Temporarily filter VEOG, bandpass 1-25 Hz
% FNYQ = EEG.srate/2;
% [bl, al] = butter(2,25/FNYQ,'low');
% [bh, ah] = butter(2,1/FNYQ,'high');
% assert(isstable(bl,al));
% assert(isstable(bh,ah));
% dataeog = filtfilt(bh,ah,filtfilt(bl,al,dataeog'))';
% assert(isstable(bh,ah));

% % Temporarily filter VEOG, lowpass 6 Hz
% [bl, al] = butter(2,6/(EEG.srate/2),'low');
% assert(isstable(bl,al));
% dataeog = filtfilt(bl,al,dataeog);
% % dataeog = abs(dataeog);

% See: RELAX_blinks_IQR_method
EOGIQR = iqr(dataeog);
EOG75P = prctile(dataeog,75);
treshold = EOG75P + 3*EOGIQR;

% Approach 1: Detect eye blinks
% eyeBlinksEpochs = dataeog>=treshold;
%
% % Epoch into 1s
% L = EEG.srate;
% N = floor(length(eyeBlinksEpochs)/L);
% dataeog = reshape(dataeog(1:N*L),L,N);
% % eyeBlinksEpochs = any(reshape(eyeBlinksEpochs(1:N*L),L,N));
% eyeBlinksEpochs = reshape(eyeBlinksEpochs(1:N*L),L,N);
%
% % Binks should last at least this long [ms]
% mspersamp = 1000/EEG.srate;
% minpopdur = round(50/mspersamp);
% eyeBlinksEpochs = (find_maxduration(eyeBlinksEpochs')>minpopdur)';
%
% badtrl_all = find(eyeBlinksEpochs);
%
% fh = figure; hold on;
% F = (0:size(dataeog,1)-1)./EEG.srate;
% plot(F,dataeog(:,badtrl_all));

% Approach 2: Detect eye blinks
times = EEG.times/1000; % EEGLAB time is in [ms]
[qrspeaks,locs] = findpeaks(dataeog,times,'MinPeakHeight',treshold);
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

EOG = NaN(NEOG,N);
for i = 1:NEOG
    if length(badEpoch2(i,1):badEpoch2(i,2))==N
        EOG(i,:) = dataeog(badEpoch2(i,1):badEpoch2(i,2));
    end
end
EOG(isnan(EOG(:,1)),:) = [];
mEOG = median(EOG,1);

fh = figure; hold on;
F = (0:size(EOG,2)-1)./EEG.srate;
F = F-mean(F);
plot(F,EOG','LineWidth',1.2);
set(gca, 'ColorOrder',brewermap(size(EOG,1),'BuGn'));
plot(F,mEOG,'Color',[0.8 0.1 0.1],'LineWidth',3);
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