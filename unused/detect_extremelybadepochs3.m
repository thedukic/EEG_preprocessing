function EEG = detect_extremelybadepochs3(EEG,cfgbch)
%
% Find extreme periods
fprintf('\nDetecting extremely bad segments of the data...\n');
fprintf('These data will be excluded from further analysis.\n');

% Select these
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
chaneog  = strcmp({EEG.chanlocs.labels},'VEOG');
dataeeg  = EEG.data(chaneeg,:);
dataeog  = EEG.data(chaneog,:);

%% Absolute threshold to identify absolute amplitude extreme values:
% More than +- 350-500 uV
extremeMask = any(abs(dataeeg)>350);

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

%% ASR
chaneeg  = strcmp({EEG.chanlocs.type},'EEG');
extchan = {EEG.chanlocs(~chaneeg).labels};
EXT = pop_select(EEG,'channel',extchan);

originalEEG = EEG;
EEG = pop_clean_rawdata(pop_select(EEG,'nochannel',extchan),'FlatlineCriterion','off','ChannelCriterion','off','LineNoiseCriterion','off','Highpass','off', ...
    'BurstCriterion',cfgbch.asr,'WindowCriterion','off','BurstRejection','on','Distance','Euclidian');

survivedDataIdx = EEG.etc.clean_sample_mask;

originalEEG = pop_select(originalEEG,'nochannel',extchan);
vis_artifacts(EEG,originalEEG);

fprintf('Done using ASR.\n');

% =========================================================================
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