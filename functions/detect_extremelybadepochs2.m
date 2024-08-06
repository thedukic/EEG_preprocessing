function EEG = detect_extremelybadepochs2(EEG)
%
% Find extreme periods based on:
% 1. Very high amplitudes (e.g. >400-500 uV)
% 2. A lot of EMG across many electrodes (e.g. 1/4 of the dataset)
%
% Very brief periods of "good" data (e.g. <1s) that are sandwitched
% between "bad" data are also marked as "bad"
%

fprintf('\nDetecting extremely bad segments of the data...\n');
fprintf('These data will be excluded from MWF, ASR and ICA and from the final dataset.\n');

% Select these
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
dataeeg  = EEG.data(chaneeg,:);
dataeog  = EEG.data(chaneog,:);

%% Absolute threshold to identify absolute amplitude extreme values:
% More than +- 500 uV
extremeMask = any(abs(dataeeg)>500);

if any(extremeMask)
    mspersamp = 1000/EEG.srate;
    popshift  = round(25/mspersamp);
    extremeMask = movmean(extremeMask,popshift)>0;

    % jump = find(diff([false, extremeMask, false])~=0);
    % mspersamp = 1000/EEG.srate;
    % popshift  = round(50/mspersamp);
    % extremeNoiseEpochs3 = [jump(1:2:end)-popshift; jump(2:2:end)+popshift]';
    % extremeNoiseEpochs3(extremeNoiseEpochs3<1) = 1;
    % extremeNoiseEpochs3(extremeNoiseEpochs3>EEG.pnts) = EEG.pnts;
else
    extremeMask = false(1,size(dataeeg,2));
end

% Log
EEG.ALSUTRECHT.extremeNoise.absoluteAmplitudeExceededThreshold = extremeMask;

%% Strong EMG

% 1. EOG: Temporarily filter, lowpass 6 Hz
[bl, al] = butter(2,15/(EEG.srate/2),'low');
assert(isstable(bl,al));
dataeog = filtfilt(bl,al,dataeog')';

treshold = prctile(dataeog,75,2) + 3*iqr(dataeog,2);
mask = dataeog>treshold;

data2 = dataeog;
data2(~mask)  = 0;
data2(mask) = 20;

data2 = conv(data2,ones(1,8),'same');
data2 = movmean(data2',16)';

% % EOG mask looks like it is lagging wrt EEG mask by ~d samples
% d = 20;
% data2(1:d) = [];
% data2 = [data2, data2(end)*ones(1,d)];
% assert(length(data2)==length(dataeog));

% 2. EEG: Temporarily filter, highpass 25 Hz
[bl, al] = butter(5,70/(EEG.srate/2),'high');
assert(isstable(bl,al));
dataeeg = filtfilt(bl,al,dataeeg')';
% dataeeg = [zeros(sum(chaneeg),1), diff(dataeeg')'];
dataeeg = abs(dataeeg);

% Zero out parts around EOG
dataeeg(:,data2>0) = 0;
% dataeeg = movmean(dataeeg',64)';

% % Check
% EEG0 = EEG;
% EEG0.data(chaneeg,:) = dataeeg;
% vis_artifacts(EEG,EEG0);

dataTmp = dataeeg(:,data2==0);
dataTmp(dataTmp==0) = [];
treshold = prctile(dataTmp(:),50);
dataeeg(dataeeg<treshold) = 0;

dataTmp = dataeeg(:,data2==0);
dataTmp(dataTmp==0) = [];
treshold = prctile(dataTmp(:),75) + 3*iqr(dataTmp(:));
mask = dataeeg>treshold;

data1 = dataeeg;
data1(~mask) = 0;
data1(mask)  = 1;

% for i = 1:size(data1,1)
%     data1(i,:) = conv(data1(i,:),ones(1,10),'same');
% end
data1 = movmean(data1',64)';

% Mask
% extremeMask = data1>0;
% extremeMask = 100*mean(extremeMask,1);
% extremeMask(extremeMask<25) = 0;

% At least 25% of electrodes must be affected
dataTmp = data1(:,data2==0);
dataTmp(dataTmp==0) = [];
treshold = prctile(dataTmp,50);

% h = histogram(dataTmp);
% % Retrieve some properties from the histogram
% V = h.Values;
% E = h.BinEdges;
% % Use islocalmax
% L = islocalmax(V);
% % Find the centers of the bins that islocalmax identified as peaks
% left = E(L);
% right = E([false L]);
% center = (left + right)/2;
% % Plot markers on those bins
% figure; plot(V);
% figure; plot(center, V(L), 'o');
% % figure; histogram(dataTmp);

extremeMask = data1>treshold;
% extremeMask = data1>0.1*max(data1(:,~data2),[],"all");
extremeMask = 100*mean(extremeMask,1);
extremeMask(extremeMask<20) = 0;
extremeMask = movmean(extremeMask,64)>0;

% extremeMask = sum(data1,1);
% treshold = prctile(extremeMask,75,2) + 2*iqr(extremeMask,2);
% extremeMask(extremeMask<treshold) = 0;
% extremeMask(extremeMask>=treshold) = 1;

% % Check
% EEG0 = EEG;
% EEG0.data = EEG0.data*0;
% EEG0.data(chaneeg,:) = 10*data1;
% EEG0.data(chaneog,:) = 10*data2;
% EEG0.data(end,:)     = 500*extremeMask;
% vis_artifacts(EEG,EEG0);

% Log
EEG.ALSUTRECHT.extremeNoise.muscleExceededThreshold = extremeMask;

%% Finish and log
% Cobine masks
EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections = ...
    EEG.ALSUTRECHT.extremeNoise.absoluteAmplitudeExceededThreshold | ...
    EEG.ALSUTRECHT.extremeNoise.muscleExceededThreshold;

% Epoch into 1s, needed only for EMG MWF
extremeMask0 = EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections;

if any(extremeMask0)
    % Make sure good epochs have minimal length
    % too short good data is likely a miss in between two bad epochs
    jump = find(diff([false, ~extremeMask0, false])~=0);
    goodEpochs = [jump(1:2:end); jump(2:2:end)-1]';
    goodEpochsDur = diff(goodEpochs')+1;
    mindur = EEG.srate; % 1s
    goodEpochs(goodEpochsDur<mindur,:) = [];

    % New noise mask
    extremeMask = true(size(extremeMask0));
    for i = 1:size(goodEpochs,1)
        extremeMask(goodEpochs(i,1):goodEpochs(i,2)) = false;
    end

    % % Check
    % EEG0 = EEG;
    % EEG0.data = EEG0.data*0;
    % EEG0.data(chaneeg,:) = 100*data1;
    % EEG0.data(chaneog,:) = 10*data2;
    % EEG0.data(end,:)     = 500*extremeMask;
    % vis_artifacts(EEG,EEG0);

    % Get epoch start/stop samples
    jump = find(diff([false, extremeMask, false])~=0);
    extremeNoiseEpochs3 = [jump(1:2:end); jump(2:2:end)-1]';

    % Makes sure that the two masks are equivalent
    assert(sum(diff(extremeNoiseEpochs3'))+size(extremeNoiseEpochs3,1) == sum(extremeMask));
else
    % No extreme noise found...
    extremeMask = extremeMask0;
    extremeNoiseEpochs3 = [];
end

L = EEG.srate;
N = floor(length(extremeMask)/L);
extremeNoiseEpochs2 = any(reshape(extremeMask(1:N*L),L,N));

% Log
EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier = mean(extremeMask);
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs0                 = extremeMask0;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1                 = extremeMask;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2                 = extremeNoiseEpochs2;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3                 = extremeNoiseEpochs3;

fprintf('Percentage of extremly bad EEG: %1.2f\n', EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Extremly bad epochs\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'These data will be excluded from MWF and ICA.\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected: %1.2f\n', EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);

end