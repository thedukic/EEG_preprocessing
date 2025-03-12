function EEG = epoch_rsdata2(EEG,epochLength,epochOverlap)
% epoch_rsdata finds the useful Epochs with length 'EpochLength' and overlap 'epochOverlap'
%
% SDukic, November 2024
% Can be done better so that epoching is possible after bad epoch rejection
%
% =========================================================================

fprintf('Epoching resting-state data (L = %d s, overlap %1.2f)...\n',epochLength,epochOverlap);

% Check
% assert(size(EEG.data,2)==length(~maskGood));
nsmp    = round(epochLength*EEG.srate);
nshift  = round((1-epochOverlap)*nsmp);
if nshift<=0, error('the overlap is too large'); end

% Extract the RS masks
maskGood = ~EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochsRSFinal;
maskRS = EEG.ALSUTRECHT.blockinfo.rs_mask;

% % Not needed, this is now done immediately when the extreme epochs are removed
% % Fix the masks as the extreme epochs have been removed already!
% % Find start/stop of good periods
% jump = find(diff([false, maskGood, false])~=0);
% jumpStop = jump(2:2:end)-1; jumpStop(end) = [];
%
% % Mark one sample just before the bad segment
% % 1   1   1   1   1   1   0   0   1   1 ---> 6
% maskGood(jumpStop) = false;
%
% maskGood(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochsRSFinal) = [];
% maskRS(:,EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochsRSFinal) = [];

% assert(length(jumpStop) == sum(~maskGood));
assert(size(EEG.data,2) == length(maskGood));
assert(size(EEG.data,2) == length(maskRS));

% Add boundaries which could be:
% 1. Between appended blocks
% 2. Due to CMS dropouts or extrme outliers being removed
isBoundary = ismember({EEG.event.type},'boundary');
latencies  = round([EEG.event.latency]);
latencies  = latencies(isBoundary);
latencies(latencies==1 | latencies==EEG.pnts-1) = [];
maskGood(latencies) = false;

% Chunk good data only, per block and per good period
NBLK = length(EEG.ALSUTRECHT.blockinfo.eo_mask);
NTRL = NaN(NBLK,1);
validEpochs  = cell(NBLK,1);
leftoverData = cell(NBLK,1);

for i = 1:NBLK
    % Good data in this block
    maskGoodTmp = maskGood(maskRS(i,:));

    % The first sample of a block should be false
    % Correct that as it can be used
    if i > 1
        assert(~maskGoodTmp(1));
        maskGoodTmp(1) = true;
    end

    % Find start/stop of good periods
    jump = find(diff([false, maskGoodTmp, false])~=0);
    jumpStart = jump(1:2:end);
    jumpStop  = jump(2:2:end)-1;

    % Find epoch starts
    epochStarts = arrayfun(@(start, stop) start:nshift:(stop + 1 - nsmp), jumpStart, jumpStop, 'UniformOutput', false);
    validEpochs{i} = cell2mat(epochStarts);

    % Calculate leftovers
    % Chunks that are very short would lead to empty cells
    mask = ~cellfun(@(x) isempty(x),epochStarts);
    epochStarts = epochStarts(mask);
    jumpStopTmp = jumpStop(mask);
    leftoverData{i} = sum(cell2mat(cellfun(@(e, stop) stop - (e(end) + nsmp - 1), epochStarts, num2cell(jumpStopTmp), 'UniformOutput', false)));

    % Add those very short epochs to the leftovers
    if any(~mask)
        epochChunks = jumpStop - jumpStart + 1;
        leftoverData{i} = leftoverData{i} + sum(epochChunks(~mask));
    end

    NTRL(i) = length(validEpochs{i});
    fprintf('Block %d: %d epochs, %d leftover samples\n', i, NTRL(i), sum(leftoverData{i}));

    % % Epoch each period separately
    % newtrl = [];
    % newleftover = [];
    % for j = 1:length(jumpStart)
    %     thistrl = (jumpStart(j):nshift:(jumpStop(j)+1-nsmp))';
    %
    %     if isempty(thistrl)
    %         % Too short period to fit any epochs
    %         newleftover = cat(1,newleftover,durall(j));
    %     else
    %         newtrl = cat(1,newtrl,thistrl);
    %         newleftover = cat(1,newleftover,jumpStop(j)-(thistrl(end)+nsmp-1));
    %     end
    % end
    %
    % goodIndx{i} = newtrl;
    % leftover{i} = newleftover;
    % NTRL(i)     = length(goodIndx{i});
end

% Build a minimal but valid EEG.event from scratch
% We ignore previous events (boundry only?)
% as we account for them in the code above!
EEG.event = [];
EEG.event = struct('type', cell(1, sum(NTRL)), 'latency', cell(1, sum(NTRL)));

% Generate EO and EC labels
numEO = sum(EEG.ALSUTRECHT.blockinfo.eo_mask);

% Double-checks
maskEO = contains(EEG.ALSUTRECHT.subject.datablocks,'_EO');

numEO0 = sum(maskEO);
numEOlabel = 1:numEO0;
numEC0 = sum(~maskEO);
numEClabel = 1:numEC0;

% Remove those that were removed completely due to very high noise
maskRemoveblock = EEG.ALSUTRECHT.extremeNoise.maskRemoveblock;
assert(length(maskEO) == length(maskRemoveblock));

maskEO(maskRemoveblock) = [];
assert(sum(maskEO) == numEO);
assert(all(maskEO(1:numEO)));

if ~isempty(numEOlabel)
    numEOlabel(maskRemoveblock) = [];
end
if ~isempty(numEClabel)
    numEClabel(maskRemoveblock) = [];
end

% Not correct if some blocks were completely removed
% EO2 will be EO1 if EO1 is removed
EO_labels = arrayfun(@(x) ['EO' num2str(x)], numEOlabel, 'UniformOutput', false);
EC_labels = arrayfun(@(x) ['EC' num2str(x)], numEClabel, 'UniformOutput', false);

% Combine labels in block order
RS_labels = cell(1, NBLK);
RS_labels(EEG.ALSUTRECHT.blockinfo.eo_mask)  = EO_labels;
RS_labels(~EEG.ALSUTRECHT.blockinfo.eo_mask) = EC_labels;

cnt = 0;
for i = 1:NBLK
    % Calculate the starting index for latencies
    K = find(maskRS(i,:), 1, "first") - 1;

    % Prepare event types and latencies
    latencies = num2cell(K + validEpochs{i});

    % Create new events and use deal
    [EEG.event(cnt + (1:NTRL(i))).type] = deal(RS_labels{i});
    [EEG.event(cnt + (1:NTRL(i))).latency] = deal(latencies{:});

    % Update counter
    cnt = cnt + NTRL(i);
end

% cnt = 0;
% for i = 1:NBLK
%     K = find(maskRS(i,:),1,"first") - 1;
%     for j = 1:NTRL(i)
%         cnt = cnt+1;
%         if EEG.ALSUTRECHT.blockinfo.eo_mask(i)
%             EEG.event(cnt).type = ['EO' num2str(i)];
%         else
%             EEG.event(cnt).type = ['EC' num2str(i-sum(EEG.ALSUTRECHT.blockinfo.eo_mask))];
%         end
%
%         % MUST BE IN SAMPLES!
%         % EEG.event(cnt).latency = (K + goodIndx{i}(j,1)-1) ./ EEG.srate;
%         EEG.event(cnt).latency = K + validEpochs{i}(j);
%     end
% end

% Check
EEG = eeg_checkset(EEG);

% Epoch now
% allLabels = unique({EEG.event(:).type});
EEG = pop_epoch(EEG,RS_labels,[0 epochLength],'epochinfo','yes');

% EEG0 = pop_epoch(EEG,allLabels,[0 epochLength]./EEG.srate,'epochinfo','yes');
% figure; hold on;
% plot(EEG.data(1,1:6*256));
% plot(reshape(EEG0.data(1,:,1:3),1,[]));
% axis tight; legend('Continuous','Epoched');

% Check
EEG = eeg_checkset(EEG);
assert(sum(NTRL) == size(EEG.data,3));
assert(nsmp == size(EEG.data,2));
% vis_artifacts(EEG,EEG);

% Takes into account overlap
oldL = sum(maskGood);
lostL = sum(cell2mat(leftoverData));
R = lostL./oldL;
fprintf('Total epochs: %d\n', sum(NTRL));
fprintf('Data loss due epoching: %1.1f%%\n', 100 * R);

% Log
EEG.ALSUTRECHT.blockinfo.goodIndx   = validEpochs;
EEG.ALSUTRECHT.blockinfo.goodTrls   = NTRL;
EEG.ALSUTRECHT.blockinfo.trlLength  = epochLength;
EEG.ALSUTRECHT.blockinfo.trlOverlap = epochOverlap;
EEG.ALSUTRECHT.blockinfo.dataLostbyEpoching = R;

end