function [EEG, EMG] = remove_extremeperiods(EEG,EMG)
%
% Find extreme periods based on:
% 1. Very high amplitudes (e.g. >400-500 uV)
% 2. A lot of EMG across many electrodes (e.g. 1/4 of the dataset)
%
% Very brief periods of "good" data (e.g. <1s)
% that are sandwitched between "bad" data are also marked as "bad"
%
% TODO:
% 1. It would be better if this was done right after CMD/DRL drop-out detection
% 2.
%

fprintf('\nDetecting extremely bad segments of the data...\n');
fprintf('These data will be excluded from MWF, ASR and ICA and from the final dataset.\n');

% Select these
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
dataall  = EEG.data;
dataeeg  = EEG.data(chaneeg,:);
dataeog  = EEG.data(chaneog,:);

%% Absolute threshold to identify absolute amplitude extreme values:
% More than +- 350-500 uV
extremeMask = any(abs(dataeeg) > 350);

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

% 1. EOG: Temporarily lowpass filter
[bl, al] = butter(2,15/(EEG.srate/2),'low');
% dataeog = filtfilt(bl,al,dataeog')';
dataeog = do_filteringcore(bl,al,dataeog,EEG.event,EEG.srate);

treshold = prctile(dataeog,75,2) + 3*iqr(dataeog,2);
maskEOG = dataeog > treshold;

data2 = dataeog;
data2(~maskEOG) = 0;
data2(maskEOG)  = 20;

data2 = conv(data2,ones(1,8),'same');
data2 = movmean(data2',16)';

% 2. EEG: Temporarily highpass filter
[bh, ah] = butter(4,70/(EEG.srate/2),'high');
% dataeeg = filtfilt(bl,al,dataeeg')';
dataall = do_filteringcore(bh,ah,dataall,EEG.event,EEG.srate);

% Traces of EMG power
dataall = 10 * abs(dataall);

% Zero out parts around EOG
data2(data2 < 0) = 0;
maskEOG = data2 > 0;

dataeeg = dataall(chaneeg,:);
dataeeg(:,maskEOG) = 0;

% % Check
% EEG0 = EEG;
% EEG0.data(chaneeg,:) = dataeeg;
% vis_artifacts(EEG,EEG0);

% =========================================================================
% Try to find clusters of higher numbers
dataTmp  = dataeeg(:,~maskEOG);
treshold = prctile(dataTmp(:),99);
maskEMG = dataeeg > treshold;

% % Check
% EEG0 = EEG;
% EEG0.data(chaneeg,:) = data1;
% vis_artifacts(EEG,EEG0);

% dataTmp = dataeeg(:,~maskEOG);
% dataTmp(dataTmp==0) = [];
% treshold = prctile(dataTmp(:),75) + 3*iqr(dataTmp(:));
% mask = dataeeg>treshold;

data1 = dataeeg;
data1(~maskEMG) = 0;
data1(maskEMG)  = 1;

data1 = movmean(data1',64)';

% =========================================================================
% At least 20-25% of electrodes must be affected
Tperc = 20;
dataTmp = data1(:,~maskEOG);
dataTmp(dataTmp==0) = [];
treshold = prctile(dataTmp,90);

extremeMask = data1 > treshold;
extremeMask = 100 * mean(extremeMask,1);

% Mask1
% -> Many channels affected together with the EOG
% -> This happens due in large EMG/movement artifacts
extremeMask1 = extremeMask>=2*Tperc & maskEOG;

% Mask
extremeMask(extremeMask < Tperc) = 0;
extremeMask(extremeMask1) = 2*Tperc;
extremeMask = movmean(extremeMask, 256) > 0;

% % Check
% EEG0 = EEG;
% EEG0.data = EEG0.data*0;
% EEG0.data(chaneeg,:) = 10*data1;
% EEG0.data(chaneog,:) = 10*data2;
% EEG0.data(end,:)     = 500*extremeMask;
% vis_artifacts(EEG,EEG0);

% Log
EEG.ALSUTRECHT.extremeNoise.muscleExceededThreshold = extremeMask;

%% Finalise the mask
% Cobine masks
EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections = ...
    EEG.ALSUTRECHT.extremeNoise.absoluteAmplitudeExceededThreshold | ...
    EEG.ALSUTRECHT.extremeNoise.muscleExceededThreshold;

% Epoch into 1s, needed only for EMG MWF
extremeMask0 = EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections;

if any(extremeMask0)
    % Find start/stop of good periods
    jumpsGood = find(diff([false, ~extremeMask0, false])~=0);
    jumpsGood = reshape(jumpsGood,2,[]);

    nMax = length(extremeMask0);
    jumpsGood(jumpsGood>nMax) = nMax;

    maskDur = diff(jumpsGood) < 2*EEG.srate;
    jumpsGood(:,maskDur) = [];

    extremeMask = true(size(extremeMask0));
    for i = 1:size(jumpsGood,2)
        extremeMask(jumpsGood(1,i):jumpsGood(2,i)) = false;
    end

    % % Check
    % EEGTMP = EEG;
    % EEGTMP.data = 0 * EEGTMP.data;
    % EEGTMP.data(chaneeg,:) = 100 * data1;
    % EEGTMP.data(chaneog,:) = 10  * data2;
    % EEGTMP.data(end,:)     = 500 * extremeMask;
    % vis_artifacts(EEG,EEGTMP);

    % Get epoch start/stop samples
    jump = find(diff([false, extremeMask, false])~=0);
    extremeNoiseEpochs3 = [jump(1:2:end); jump(2:2:end)-1]';

    % Makes sure that the two masks are equivalent
    assert(sum(diff(extremeNoiseEpochs3')) + size(extremeNoiseEpochs3,1) == sum(extremeMask));
else
    % No extreme noise found
    extremeMask = extremeMask0;
    extremeNoiseEpochs3 = [];
end

L = EEG.srate;
N = floor(length(extremeMask)/L);
extremeNoiseEpochs2 = any(reshape(extremeMask(1:N*L),L,N));

%% Log
EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier = mean(extremeMask);
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs0                 = extremeMask0;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1                 = extremeMask;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs2                 = extremeNoiseEpochs2;
EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3                 = extremeNoiseEpochs3;

fprintf('Proportion of extremly bad EEG: %1.2f\n', EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Extremly bad epochs\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'These data will be excluded from MWF and ICA.\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected: %1.2f\n', EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);

% Report visually
if ~isempty(extremeNoiseEpochs3)
    report_badsegments(EEG,{extremeNoiseEpochs3},'extremeperiods');
end

%% Remove the chunks
if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)
    fprintf('\nRemoving extremely bad chunks of data...\n');
    EEG = eeg_eegrej(EEG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);

    if nargin == 2
        EMG = eeg_eegrej(EMG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    end
end

%% Fix the masks as the extreme epochs have been removed already!
if isfield(EEG.ALSUTRECHT,'blockinfo')
    maskGood = ~EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1;
    maskRS = EEG.ALSUTRECHT.blockinfo.rs_mask;

    % Find start/stop of good periods
    jump = find(diff([false, maskGood, false])~=0);
    jumpStop = jump(2:2:end)-1; jumpStop(end) = [];

    % Mark one sample just before the bad segment
    % 1   1   1   1   1   1   0   0   1   1 ---> 6
    maskGood(jumpStop) = false;

    maskGood(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = [];
    maskRS(:,EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = [];

    % Double-check
    assert(size(EEG.data,2) == length(maskGood));
    assert(size(EEG.data,2) == length(maskRS));
    assert(length(jumpStop) == sum(~maskGood));

    % Return
    EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1 = ~maskGood;
    EEG.ALSUTRECHT.blockinfo.rs_mask = maskRS;
end

end