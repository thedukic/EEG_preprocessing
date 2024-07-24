function EEG = epoch_rsdata2(EEG,epochLength,epochOverlap)
% epoch_rsdata finds the useful Epochs with length 'EpochLength' in the UV mask vector
%
% [IDX, IDXsimplevector] = bepocher(UV, EpochLength)
%
%   UV          : is a 1xn vector that represents bad data points with 0 and good data points with 1.
%   EpochLength : is the length of the data chunks (in samples) that need to be extracted from the good periods
%   Epochs      : includes the index number where good epochs with length 'EpochLength' start.
%
%   IDX         : indices of the data for extraction (Nepoch x EpochLength)
%
% SDukic, March 2024
% Can be done better so that epoching is possible after bad epoch rejection
%
% =========================================================================

fprintf('\nEpoching (L = %ds) resting-state data...\n',epochLength);

% Check
% assert(size(EEG.data,2)==length(~maskGood));
nsmp    = round(epochLength*EEG.srate);
nshift  = round((1-epochOverlap)*nsmp);
if nshift<=0, error('the overlap is too large'); end

% Chunk good data only, per block and per good period
NBLK     = length(EEG.ALSUTRECHT.blockinfo.eo_mask);
goodIndx = cell(NBLK,1);
leftover = cell(NBLK,1);
NTRL     = NaN(NBLK,1);

% Fix the masks as the extreme epochs have been removed already!
maskGood = ~EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1;
maskRS = EEG.ALSUTRECHT.blockinfo.rs_mask;

% Find start/stop of good periods
jump = find(diff([false, maskGood, false])~=0);
jumpStop  = jump(2:2:end)-1; jumpStop(end) = [];

% Mark one sample as the break point at each jump
maskGood(jumpStop) = false;

maskGood(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = [];
maskRS(:,EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = [];

assert(size(EEG.data,2)==length(maskGood));
assert(size(EEG.data,2)==length(maskRS));

for i = 1:NBLK
    % Good data in this block
    maskGoodTmp = maskGood(maskRS(i,:));

    % Find start/stop of good periods
    jump = find(diff([false, maskGoodTmp, false])~=0);
    jumpStart = jump(1:2:end);
    jumpStop  = jump(2:2:end)-1;
    durall    = jumpStop-jumpStart+1;

    % Epoch each period separately
    newtrl = [];
    newleftover = [];
    for j = 1:length(jumpStart)
        thistrl = (jumpStart(j):nshift:(jumpStop(j)+1-nsmp))';

        if isempty(thistrl)
            % Too short period to fit any epochs
            newleftover = cat(1,newleftover,durall(j));
        else
            newtrl = cat(1,newtrl,thistrl);
            newleftover = cat(1,newleftover,jumpStop(j)-(thistrl(end)+nsmp-1));
        end
    end
    goodIndx{i} = newtrl;
    leftover{i} = newleftover;
    NTRL(i)     = length(goodIndx{i});
end

% Build a minimal but valid EEG.event from scratch
% We ignore previous events (boundry only?)
% as we account for them in the code above!
EEG.event = [];

cnt = 0;
for i = 1:NBLK
    K = find(maskRS(i,:),1,"first")-1;
    for j = 1:NTRL(i)
        cnt = cnt+1;
        if EEG.ALSUTRECHT.blockinfo.eo_mask(i)
            EEG.event(cnt).type = ['EO' num2str(i)];
        else
            EEG.event(cnt).type = ['EC' num2str(i-sum(EEG.ALSUTRECHT.blockinfo.eo_mask))];
        end

        % EEG.event(cnt).latency = (K+goodIndx{i}(j,1)-1)./EEG.srate;
        EEG.event(cnt).latency = K+goodIndx{i}(j); % MUST BE IN SAMPLES!
    end
end
EEG = eeg_checkset(EEG);

% Epoch now
allLabels = unique({EEG.event(:).type});
EEG = pop_epoch(EEG,allLabels,[0 epochLength],'epochinfo','yes');

% EEG0 = pop_epoch(EEG,allLabels,[0 epochLength]./EEG.srate,'epochinfo','yes');
% figure; hold on;
% plot(EEG.data(1,1:6*256));
% plot(reshape(EEG0.data(1,:,1:3),1,[]));
% axis tight; legend('Continuous','Epoched');

% Check
EEG = eeg_checkset(EEG);
assert(sum(NTRL)==size(EEG.data,3));
assert(nsmp==size(EEG.data,2));
% vis_artifacts(EEG,EEG);

% % This is correct only if not overlapping is done!
% % % -> extreme noise + epoching
% % oldL = length(~maskGood);
% % % -> epoching only
% oldL = sum(maskGood);
% newL = prod(size(EEG.data,[2 3]));
% R = 1-newL./oldL;
% fprintf('By epoching (L = %ds), %1.2f of the data is lost.\n',epochLength./EEG.srate,R);

% Takes into account overlap
oldL = sum(maskGood);
lostL = sum(cell2mat(leftover));
R = lostL./oldL;
fprintf('Epoching resulted in %1.3f data loss.\n',R);

% Log
EEG.ALSUTRECHT.blockinfo.goodIndx   = goodIndx;
EEG.ALSUTRECHT.blockinfo.goodTrls   = NTRL;
EEG.ALSUTRECHT.blockinfo.trlLength  = epochLength;
EEG.ALSUTRECHT.blockinfo.trlOverlap = epochOverlap;
EEG.ALSUTRECHT.blockinfo.dataLostbyEpoching = R;

end