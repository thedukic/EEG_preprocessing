function EEG = epoch_rsdata2(EEG,epochLength,epochOverlap)
% epoch_rsdata finds the useful Epochs with length 'EpochLength' and overlap 'epochOverlap'
%
% SDukic, November 2024
% Can be done better so that epoching is possible after bad epoch rejection
%
% =========================================================================

fprintf('Epoching resting-state data (L = %d s, overlap %1.2f)...\n', epochLength, epochOverlap);

% Check
% assert(size(EEG.data, 2) == length(~mask_good));
nsmp    = round(epochLength * EEG.srate);
nshift  = round((1-epochOverlap) * nsmp);
if nshift <= 0, error('the overlap is too large'); end

% % Extract the RS masks
% if isfield(EEG.ALSUTRECHT, "bad_segments")
%     mask_good = ~EEG.ALSUTRECHT.bad_segments.mask_bad;
% else
%     mask_good = true(1, EEG.pnts);
% end
mask_good  = true(1, EEG.pnts);
mask_block = EEG.ALSUTRECHT.blocks.mask_block;

% % Not needed, this is now done immediately when the extreme epochs are removed
% % Fix the masks as the extreme epochs have been removed already!
% % Find start/stop of good periods
% jump = find(diff([false, mask_good, false])~=0);
% jumpStop = jump(2:2:end)-1; jumpStop(end) = [];
%
% % Mark one sample just before the bad segment
% % 1   1   1   1   1   1   0   0   1   1 ---> 6
% mask_good(jumpStop) = false;
%
% mask_good(EEG.ALSUTRECHT.bad_segments.bad_segmentsEpochsRSFinal) = [];
% maskRS(:,EEG.ALSUTRECHT.bad_segments.bad_segmentsEpochsRSFinal) = [];

% assert(length(jumpStop) == sum(~mask_good));
assert(size(EEG.data, 2) == length(mask_good));
assert(size(EEG.data, 2) == size(mask_block, 2));

% Add boundaries which could be:
% 1. Between appended blocks
% 2. Due to CMS dropouts or extrme outliers being removed
is_boundary = ismember({EEG.event.type}, 'boundary');
latencies  = round([EEG.event.latency]);
latencies  = latencies(is_boundary);
latencies(latencies==1 | latencies==EEG.pnts-1) = [];
mask_good(latencies) = false;

% Chunk good data only, per block and per good period
num_block    = size(mask_block, 1);
num_trial    = NaN(num_block, 1);
validEpochs  = cell(num_block, 1);
leftoverData = cell(num_block, 1);

for i_block = 1:num_block
    % Good data in this block
    mask_good_tmp = mask_good(mask_block(i_block, :));

    % % The first sample of a block should be false
    % % Correct that as it can be used
    % if i > 1
    %     assert(~mask_good_tmp(1));
    %     mask_good_tmp(1) = true;
    % end

    % Find start/stop of good periods
    jump = find(diff([false, mask_good_tmp, false])~=0);
    jumpStart = jump(1:2:end);
    jumpStop  = jump(2:2:end)-1;

    % Find epoch starts
    epochStarts = arrayfun(@(start, stop) start:nshift:(stop + 1 - nsmp), jumpStart, jumpStop, 'UniformOutput', false);
    validEpochs{i_block} = cell2mat(epochStarts);

    % Calculate leftovers
    % Chunks that are very short would lead to empty cells
    mask = ~cellfun(@(x) isempty(x), epochStarts);
    epochStarts = epochStarts(mask);
    jumpStopTmp = jumpStop(mask);
    leftoverData{i_block} = sum(cell2mat(cellfun(@(e, stop) stop - (e(end) + nsmp - 1), epochStarts, num2cell(jumpStopTmp), 'UniformOutput', false)));

    % Add those very short epochs to the leftovers
    if any(~mask)
        epochChunks = jumpStop - jumpStart + 1;
        leftoverData{i_block} = leftoverData{i_block} + sum(epochChunks(~mask));
    end

    num_trial(i_block) = length(validEpochs{i_block});
    fprintf('Block %d: %d epochs, %d leftover samples\n', i_block, num_trial(i_block), sum(leftoverData{i_block}));

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
    % num_trial(i)     = length(goodIndx{i});
end

% Build a minimal but valid EEG.event from scratch
% We ignore previous events (boundry only?)
% but we account for them in the code above!
EEG.event = [];
EEG.event = struct('type', cell(1, sum(num_trial)), 'latency', cell(1, sum(num_trial)));

% Generate EO and EC labels
num_eo = sum(EEG.ALSUTRECHT.blocks.mask_eo);

% Double-checks
mask_eo = contains(EEG.ALSUTRECHT.subject.datablocks, '_EO');

numEO0     = sum(mask_eo);
numEOlabel = 1:numEO0;
numEC0     = sum(~mask_eo);
numEClabel = 1:numEC0;

% % Remove those that were removed completely due to very high noise
% if isfield(EEG.ALSUTRECHT, "bad_segments")
%     mask_block_removed = EEG.ALSUTRECHT.bad_segments.maskRemoveblock;
% else
%     mask_block_removed = false(size(mask_eo));
% end
mask_block_removed = false(size(mask_eo));

assert(length(mask_eo) == length(mask_block_removed));
mask_eo(mask_block_removed) = [];
assert(sum(mask_eo) == num_eo);
assert(all(mask_eo(1:num_eo)));

if ~isempty(numEOlabel)
    numEOlabel(mask_block_removed) = [];
end
if ~isempty(numEClabel)
    numEClabel(mask_block_removed) = [];
end

% Not correct if some blocks were completely removed
% EO2 will be EO1 if EO1 is removed
EO_labels = arrayfun(@(x) ['EO' num2str(x)], numEOlabel, 'UniformOutput', false);
EC_labels = arrayfun(@(x) ['EC' num2str(x)], numEClabel, 'UniformOutput', false);

% Combine labels in block order
RS_labels = cell(1, num_block);
RS_labels(EEG.ALSUTRECHT.blocks.mask_eo)  = EO_labels;
RS_labels(~EEG.ALSUTRECHT.blocks.mask_eo) = EC_labels;

cnt = 0;
for i_block = 1:num_block
    % Calculate the starting index for latencies
    K = find(mask_block(i_block,:), 1, "first") - 1;

    % Prepare event types and latencies
    latencies = num2cell(K + validEpochs{i_block});

    % Create new events and use deal
    [EEG.event(cnt + (1:num_trial(i_block))).type] = deal(RS_labels{i_block});
    [EEG.event(cnt + (1:num_trial(i_block))).latency] = deal(latencies{:});

    % Update counter
    cnt = cnt + num_trial(i_block);
end

% Check
EEG = eeg_checkset(EEG);

% Epoch now
% allLabels = unique({EEG.event(:).type});
EEG = pop_epoch(EEG, RS_labels, [0 epochLength], 'epochinfo', 'yes');

% EEG0 = pop_epoch(EEG,allLabels,[0 epochLength]./EEG.srate,'epochinfo','yes');
% figure; hold on;
% plot(EEG.data(1,1:6*256));
% plot(reshape(EEG0.data(1,:,1:3),1,[]));
% axis tight; legend('Continuous','Epoched');

% Check
EEG = eeg_checkset(EEG);
assert(sum(num_trial) == size(EEG.data, 3));
assert(nsmp == size(EEG.data,2));
% vis_artifacts(EEG,EEG);

% Takes into account overlap
oldL = sum(mask_good);
lostL = sum(cell2mat(leftoverData));
R = lostL./oldL;
fprintf('Total epochs: %d\n', sum(num_trial));
fprintf('Data loss due epoching: %1.1f%%\n', 100 * R);

% Log
EEG.ALSUTRECHT.trials.goodIndx   = validEpochs;
EEG.ALSUTRECHT.trials.goodTrls   = num_trial;
EEG.ALSUTRECHT.trials.trlLength  = epochLength;
EEG.ALSUTRECHT.trials.trlOverlap = epochOverlap;
EEG.ALSUTRECHT.trials.dataLostbyEpoching = R;

end