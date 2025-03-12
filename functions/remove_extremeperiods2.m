function [EEG, flagExclude] = remove_extremeperiods2(EEG)
%
% Find extreme periods based on:
% 1. Very high amplitudes (e.g. >350 uV)
% 2. A lot of EMG across many electrodes
%
% Very brief periods of "good" data (e.g. <1s)
% that are sandwitched between "bad" data are also marked as "bad"
%
% Prevent data being excluded just because one electrode is going rogue
% An electrode that is taken out often also has very extreme values (unstable)
% It is wise to remove these flat electrodes before running this function
%
% TODO:
% 1. Maybe very brief bad data should not be removed? If false positive
% 2.
%

fprintf('\n================================\n');
fprintf('Detecting extremely bad segments\n');
fprintf('================================\n');
fprintf('These data will be excluded from further analysis.\n');

% =========================================================================
% Detect extremly bad chunks of data
% =========================================================================
% EEG channels
chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');

% Monopolar VEOG channels
eleclabels = {EEG(1).chanlocs.labels};
m1 = ismember(eleclabels,'VEOGS');
m2 = ismember(eleclabels,'VEOGI');

fs = EEG(1).srate;
NBLK = length(EEG);
extremeMaskFinalVisual = cell(NBLK,1);
badChans = false(sum(chaneeg),NBLK);

for i_blk = 1:NBLK
    fprintf('Block %d: Checking for extreme noise...\n',i_blk);

    % dataall = EEG(i_blk).data;
    dataeeg = EEG(i_blk).data(chaneeg,:);
    dataeog = EEG(i_blk).data(m1,:) - EEG(i_blk).data(m2,:);
    % dataeog = EEG(i_blk).data(chaneog,:);

    % Demean
    dataeeg = dataeeg - mean(dataeeg,2);
    dataeog = dataeog - mean(dataeog);

    % Temporarily highpass filter
    [bh, ah] = butter(2,0.5/(fs/2),'high');
    dataeeg = do_filteringcore(bh,ah,dataeeg,EEG(i_blk).event,fs);
    dataeog = do_filteringcore(bh,ah,dataeog,EEG(i_blk).event,fs);

    % =====================================================================
    % Absolute threshold to identify absolute amplitude extreme values:
    % More than +- 350-500 uV
    % =====================================================================
    % T = 350 uV can also capture very strong EMG here, when fs = 512
    extremeMask1 = abs(dataeeg) > 350;

    % Avearge extreme per channel
    extremeMask1ChanAvg = mean(extremeMask1,2);

    % Remove them from further checks to avoid remove unnecessary data
    badChans(:,i_blk) = extremeMask1ChanAvg > 0.15;
    extremeMask1(badChans(:,i_blk),:) = false;
    dataeeg(badChans(:,i_blk),:) = [];

    % Make the mask now
    extremeMask1 = any(extremeMask1);

    % Report
    if any(badChans(:,i_blk))
        fprintf('Channels with extreme values found (N = %d).\n', sum(badChans(:,i_blk)));
    end

    if any(extremeMask1)
        mspersamp = 1000/fs;
        popshift  = round(25/mspersamp);
        extremeMask1 = movmean(extremeMask1,popshift) > 0;

        % jump = find(diff([false, extremeMask, false])~=0);
        % mspersamp = 1000/fs;
        % popshift  = round(50/mspersamp);
        % extremeNoiseEpochs3 = [jump(1:2:end)-popshift; jump(2:2:end)+popshift]';
        % extremeNoiseEpochs3(extremeNoiseEpochs3<1) = 1;
        % extremeNoiseEpochs3(extremeNoiseEpochs3>EEG(i_blk).pnts) = EEG(i_blk).pnts;
    else
        extremeMask1 = false(1,size(dataeeg,2));
    end

    % =====================================================================
    % Strong EMG
    % =====================================================================

    % 1. EOG: Temporarily lowpass filter
    [bl, al] = butter(2,15/(fs/2),'low');
    dataeog = do_filteringcore(bl,al,dataeog,EEG(i_blk).event,fs);

    % Estimate the treshold
    treshold = prctile(dataeog,75,2) + 3*iqr(dataeog,2);

    % Minimum treshold set
    % -> Prevents from noise being detected as blinks
    % -> Especially good for EC blocks!
    treshold = max(treshold, 150);
    maskEOG = dataeog > treshold;

    data2 = dataeog;
    data2(~maskEOG) = 0;
    data2(maskEOG)  = 20;

    data2 = conv(data2,ones(1,8),'same');
    data2 = movmean(data2',16)';

    % 2. EEG: Temporarily highpass filter
    [bh, ah] = butter(4,70/(fs/2),'high');
    dataeeg = do_filteringcore(bh,ah,dataeeg,EEG(i_blk).event,fs);

    % Traces of EMG power
    dataeeg = 10 * abs(dataeeg);

    % Zero out parts around EOG
    data2(data2 < 0) = 0;
    maskEOG = data2 > 0;

    dataeegPlot = dataeeg;
    dataeeg(:,maskEOG) = 0;

    % % Check
    % EEG0 = EEG;
    % EEG0.data(chaneeg,:) = dataeeg;
    % vis_artifacts(EEG,EEG0);

    % ======================================
    % Try to find clusters of higher numbers
    % ======================================
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

    data1 = movmean(data1', fs/4)';

    % ======================================
    % At least 20-25% of electrodes must be affected
    % ======================================
    Tperc = 20;
    dataTmp = data1(:,~maskEOG);
    dataTmp(dataTmp==0) = [];
    treshold = prctile(dataTmp, 90);

    extremeMask2 = data1 > treshold;
    extremeMask2 = 100 * mean(extremeMask2,1);

    % Mask1
    % -> Many channels affected together with the EOG
    % -> This happens due in large EMG/movement artifacts
    extremeMaskTmp2 = extremeMask2>=2*Tperc & maskEOG;

    % Mask
    extremeMask2(extremeMask2 < Tperc) = 0;
    extremeMask2(extremeMaskTmp2) = 2*Tperc;
    extremeMask2 = movmean(extremeMask2, fs) > 0;

    % % Check
    % EEG00 = EEG(i_blk);
    % EEG00.data(chaneeg,:) = dataeeg;
    %
    % EEG0 = EEG(i_blk);
    % EEG0.data = EEG0.data*0;
    % EEG0.data(chaneeg,:) = 10*data1;
    % EEG0.data(sum(chaneeg)+1,:) = 10*data2;
    % EEG0.data(sum(chaneeg)+2,:) = 500*extremeMask2;
    % vis_artifacts(EEG00,EEG0);

    % =====================================================================
    % Combine the masks
    % =====================================================================
    % Make sure that  we dont exclude data that is actually large EOG
    extremeMask1(maskEOG) = false;

    % Combine
    extremeMaskCombined = extremeMask1 | extremeMask2;

    % =====================================================================
    % Epoch into 1s, needed for EMG MWF only
    % =====================================================================
    if any(extremeMaskCombined)
        % % Find start/stop of good periods
        % jumpsGood = find(diff([false, ~extremeMaskCombined, false])~=0);
        % jumpsGood = reshape(jumpsGood,2,[]);
        %
        % nMax = length(extremeMaskCombined);
        % jumpsGood(jumpsGood > nMax) = nMax;
        %
        % % If the "good" period is too short -> mark it as bad
        % maskDur = diff(jumpsGood) < 2*fs;
        % jumpsGood(:,maskDur) = [];
        %
        % extremeMaskFinal = true(size(extremeMaskCombined));
        % for i_chunk = 1:size(jumpsGood,2)
        %     extremeMaskFinal(jumpsGood(1,i_chunk):jumpsGood(2,i_chunk)) = false;
        % end

        % Fix too short good periods and then the same for bad periods
        extremeMaskFinal = fix_mask(extremeMaskCombined,fs,1);
        extremeMaskFinal = fix_mask(extremeMaskFinal,fs,2);

        % % Check
        % chaneegTmp = 1:sum(chaneeg);
        % chaneegTmp(badChans(:,i_blk)) = [];
        %
        % EEGTMP1 = EEG(i_blk);
        % EEGTMP1.data(chaneegTmp,:) = dataeegPlot;
        % EEGTMP1.data(129,:)        = dataeog;
        % EEGTMP1.data(130,:)        = randn(size(dataeog));
        % EEGTMP1.data(131:end,:)    = 0 * EEGTMP1.data(131:end,:);
        %
        % EEGTMP2 = EEG(i_blk);
        % EEGTMP2.data               = 0 * EEGTMP2.data;
        % EEGTMP2.data(chaneegTmp,:) = 500 * data1; % EOG not visible due to filtering
        % EEGTMP2.data(129,:)        = 10  * data2;
        % EEGTMP2.data(130,:)        = 500 * extremeMaskFinal;
        % EEGTMP1.data(131:end,:)    = 0 * EEGTMP1.data(131:end,:);
        % vis_artifacts(EEGTMP1,EEGTMP2);

        % Get epoch start/stop samples
        jump = find(diff([false, extremeMaskFinal, false])~=0);
        extremeNoiseEpochs3 = [jump(1:2:end); jump(2:2:end)-1]';

        % Makes sure that the two masks are equivalent
        assert(sum(diff(extremeNoiseEpochs3')) + size(extremeNoiseEpochs3,1) == sum(extremeMaskFinal));
    else
        % No extreme noise found
        extremeMaskFinal = extremeMaskCombined;
        extremeNoiseEpochs3 = [];
    end

    L = fs;
    N = floor(length(extremeMaskFinal) / L);
    extremeNoiseEpochs2 = any(reshape(extremeMaskFinal(1:N*L),L,N));

    % =====================================================================
    % Make a mask for resting-state data
    % =====================================================================
    maskGoodRS = ~extremeMaskFinal;

    % Find start/stop of good periods
    jump = find(diff([false, maskGoodRS, false])~=0);
    if ~isempty(jump)
        jumpStop = jump(2:2:end)-1;
        jumpStop(end) = [];

        % Mark one sample just before the bad segment
        % 1   1   1   1   1   1   0   0   1   1 ---> 6
        maskGoodRS(jumpStop) = false;
        maskGoodRS(extremeMaskFinal) = [];

        % Double-check
        assert(length(jumpStop) == sum(~maskGoodRS));
    else
        % The whole data block is bad!
        maskGoodRS(extremeMaskFinal) = [];
    end

    % =====================================================================
    % Log
    % =====================================================================
    EEG(i_blk).ALSUTRECHT.extremeNoise.absoluteAmplitudeExceededThreshold = extremeMask1;
    EEG(i_blk).ALSUTRECHT.extremeNoise.muscleExceededThreshold            = extremeMask2;
    EEG(i_blk).ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections   = extremeMaskCombined;

    EEG(i_blk).ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier = mean(extremeMaskFinal);
    EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochs0                 = extremeMaskCombined;
    EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochs1                 = extremeMaskFinal;
    EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochs2                 = extremeNoiseEpochs2;
    EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochs3                 = extremeNoiseEpochs3;
    EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochsRS                = ~maskGoodRS;

    extremeMaskFinalVisual{i_blk} = extremeNoiseEpochs3;

    fprintf('Proportion of extremly bad EEG: %1.2f\n', EEG(i_blk).ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);

    fprintf(EEG(i_blk).ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
    fprintf(EEG(i_blk).ALSUTRECHT.subject.fid,'Block %d: Extremly bad epochs\n',i_blk);
    fprintf(EEG(i_blk).ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
    fprintf(EEG(i_blk).ALSUTRECHT.subject.fid,'These data will be excluded from further analysis.\n');
    fprintf(EEG(i_blk).ALSUTRECHT.subject.fid,'Detected: %1.2f\n', EEG(i_blk).ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier);

end

% =========================================================================
% Report visually
% =========================================================================
report_badsegments(EEG,extremeMaskFinalVisual,'extremeperiods');

% =========================================================================
% Remove bad chunks
% =========================================================================
% If any found, then remove them
maskRemoveblock = false(NBLK,1);

if any(~cellfun(@isempty, extremeMaskFinalVisual))
    fprintf('Extremely bad chunks detected. Removing them now...\n');

    for i_blk = 1:NBLK
        maskBlockTmp1 = EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochs3;
        maskBlockTmp2 = EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochsRS;

        % Remove if not empty
        if ~isempty(maskBlockTmp1)
            if isempty(maskBlockTmp2)
                assert(diff(maskBlockTmp1)+1 == EEG(i_blk).pnts);
                maskRemoveblock(i_blk) = true;
            else
                EEG(i_blk) = eeg_eegrej(EEG(i_blk), maskBlockTmp1);
            end
        end

        % Make sure that the block masks are the right size
        % -> Only needed for RS data
        if ~maskRemoveblock(i_blk)
            assert(size(EEG(i_blk).data,2) == length(EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochsRS));
        end
    end

    % Deal with very bad blocks of data
    if any(maskRemoveblock)
        warning('Removing %d complete recording blocks!', sum(maskRemoveblock));
        EEG(maskRemoveblock) = [];
        NBLK = length(EEG);

    end
else
    fprintf('Great, no extreme noise found at all!\n');
end

% =========================================================================
% RS mask
% =========================================================================
% If any, then remove them
maskRS = [];
L = NaN(NBLK,1);
for i_blk = 1:NBLK
    maskTmp    = EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochsRS;
    maskTmp(1) = true;
    maskRS     = [maskRS, maskTmp];
    L(i_blk)   = EEG(i_blk).pnts;
end

% Make sure that the total mask is the same length as the total data
assert(length(maskRS) == sum(L));

% Store
for i_blk = 1:NBLK
    EEG(i_blk).ALSUTRECHT.extremeNoise.extremeNoiseEpochsRSFinal = maskRS;
    EEG(i_blk).ALSUTRECHT.extremeNoise.extremeBadChannels        = badChans;
    EEG(i_blk).ALSUTRECHT.extremeNoise.maskRemoveblock           = maskRemoveblock;
end

% =========================================================================
% Check if there is enough data for further analysis
% =========================================================================
% This code tries to prevents the pipeline from breaking if data is very bad
% This should never be the case but it can happen
Ndata = 0;
for i_blk = 1:NBLK
    Ndata = Ndata + EEG(i_blk).pnts;
end

% Minimum of 90s of data
if Ndata < 90 * EEG(1).srate
    flagExclude = true;
else
    flagExclude = false;
end

end

% =========================================================================
% Helper functions
% =========================================================================
function extremeMaskFinal = fix_mask(extremeMaskCombined,fs,flagOpt)
% Removes too short bad/good periods

if flagOpt == 1
    % Fix good periods
    jumpsGood = find(diff([false, ~extremeMaskCombined, false]) ~= 0);
    minL = 2*fs;
elseif flagOpt == 2
    % Fix bad periods
    jumpsGood = find(diff([false, extremeMaskCombined, false]) ~= 0);
    % jumpsGood = find(diff([false, ~extremeMaskCombined, false]) ~= 0);
    minL = fs/4;
end

% Find start/stop of
% 1. good periods
% 2. bad periods
jumpsGood = reshape(jumpsGood,2,[]);

nMax = length(extremeMaskCombined);
jumpsGood(jumpsGood > nMax) = nMax;

% If the good/bad period is too short -> remove it
maskDur = diff(jumpsGood) < minL;
jumpsGood(:,maskDur) = [];

if flagOpt == 1
    % 1. mark good periods
    extremeMaskFinal = true(size(extremeMaskCombined));
    for i_chunk = 1:size(jumpsGood,2)
        extremeMaskFinal(jumpsGood(1,i_chunk):jumpsGood(2,i_chunk)) = false;
    end

elseif flagOpt == 2
    % 2. mark bad periods
    extremeMaskFinal = false(size(extremeMaskCombined));
    for i_chunk = 1:size(jumpsGood,2)
        extremeMaskFinal(jumpsGood(1,i_chunk):jumpsGood(2,i_chunk)) = true;
    end

end

end