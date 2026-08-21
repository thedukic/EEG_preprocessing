function EEG = detect_dropouts(EEG, cfg)
% This function detects CMS/DRL dropouts.
% As a side effect, it also excludes data with large deviations,
% such as large movements or coughing.

fprintf('\n================================\n');
fprintf('Detecting CMS/DRL drop-outs\n');
fprintf('================================\n');

fprintf('Temporarily filtering the EEG electrodes for better detection of CMS out of range.\n');
eegchan = strcmp({EEG(1).chanlocs.type}, 'EEG');
NCHN = sum(eegchan);

% Remove slow drifts
EEGTMP = filter_signal(EEG, [], [1 2], 1:length(eegchan), 'eeglab');
% vis_artifacts(EEGTMP(2),EEGTMP(2));
% EEGTMP = EEG;

NBLK = length(EEG);
flagDropout = false(NBLK, 1);
jumpsBadAll = cell(NBLK, 1);

for i_block = 1:NBLK
    tmpData1 = EEGTMP(i_block).data(1:NCHN,:);
    smoothFactor = round(2*EEG(i_block).srate);
    % figure; histogram(tmp(:));

    % 1. Extreme voltage: CMS out of range
    maskBadVoltage = abs(tmpData1) > 350;
    maskBadVoltage = sum(maskBadVoltage, 1);
    maskBadVoltage = movmean(maskBadVoltage, smoothFactor);
    maskBadVoltage = maskBadVoltage > 10;

    % 2. No signal: unsuccessful data recovery
    dataDiff = abs(diff(tmpData1,1,2));

    % This accounts for long dropouts, which leads to the threshold being 0!
    % dataDiffMedian = median(dataDiff,2);
    dataDiffMedian = NaN(NCHN, 1);
    for i_channel = 1:NCHN
        tmpData2 = dataDiff(i_channel, :);
        dataDiffMedian(i_channel) = median(tmpData2(tmpData2 ~= 0));
    end

    DiffMedianthreshold = prctile(dataDiffMedian, 2);
    maskBadZero = dataDiff <= DiffMedianthreshold;
    maskBadZero = sum(maskBadZero, 1);
    maskBadZero = movmean(maskBadZero, smoothFactor);
    % figure; plot(maskBadZero);

    % The same as above, solves the problem of long dropouts, when trashold becomes 128
    maskBadZeroTmp = maskBadZero(maskBadZero < NCHN);

    % threshold
    ZIQR = iqr(maskBadZeroTmp);
    Z75P = prctile(maskBadZeroTmp, 75);
    threshold = Z75P + 10 * ZIQR;
    maskBadZero = maskBadZero > threshold;

    % tmp(maskBad) = 0;
    % EEG2 = EEG;
    % EEG2(2).data(1:NCHN,:) = tmp;
    % vis_artifacts(EEG2(2),EEG(2));

    % Combine
    maskBad = maskBadVoltage | [false maskBadZero];
    % maskBad = movmean(maskBad,smoothFactor);
    % maskBad = maskBad>0;
    % figure; plot(maskBad);

    % tmp(:,maskBad) = 0;
    % EEG2 = EEGTMP;
    % EEG2(2).data(1:NCHN,:) = tmp;
    % EEGTMP(2).data(1:NCHN,:) = tmp;
    % vis_artifacts(EEG2(2),EEGTMP(2));

    % Find start/stop of bad periods
    jumpsBad = find(diff([false, maskBad, false]) ~= 0);
    jumpsBad = reshape(jumpsBad, 2, [])';

    if ~isempty(jumpsBad)
        fprintf('Block %d: Large deviations found in the data (N = %d chunks)!\n', i_block, size(jumpsBad, 1));

        jumpsBad(:, 1) = jumpsBad(:,1) - smoothFactor;
        jumpsBad(:, 2) = jumpsBad(:,2) + smoothFactor;
        jumpsBad(jumpsBad < 1) = 1;
        nMax = size(tmpData1, 2);
        jumpsBad(jumpsBad > nMax) = nMax;

        % maskNew = false(1,nMax);
        % for j = 1:size(jumpsBad,1)
        %     maskNew(jumpsBad(j,1):jumpsBad(j,2)) = true;
        % end
        %
        % tmp(:,maskNew) = 0;
        % EEG2 = EEGTMP;
        % EEG2(i).data(1:NCHN,:) = tmp;
        % EEG2(i).data(NCHN:end,maskNew) = 0;
        % vis_artifacts(EEG2(i),EEGTMP(i));

        % Samples -> s
        jumpsBadSec = (jumpsBad-1) ./ EEG(i_block).srate;
        disp(jumpsBadSec);

        % EEG(i) = eeg_eegrej(EEG(i), jumpsBad);
        % % EEG(i) = pop_select(EEG(i),'rmpoint',jumpsBad);

        jumpsBadAll{i_block} = jumpsBad;
        flagDropout(i_block) = true;
    else
        fprintf('Block %d: Nice, no large deviations found!\n',i_block);
    end
end

% EEGTMP = filter_signal(EEG,[],[1 8],1:length(eegchan),'eeglab');
% vis_artifacts(EEGTMP(1),EEGTMP(1));

% =====================================================================
% Log
NDropouts = cellfun(@(x) size(x, 1), jumpsBadAll);

for i_block = 1:NBLK
    EEG(i_block).ALSUTRECHT.cmsDropouts.flagDropout = flagDropout;
    EEG(i_block).ALSUTRECHT.cmsDropouts.jumpsBadAll = jumpsBadAll;
    EEG(i_block).ALSUTRECHT.cmsDropouts.NDropouts   = NDropouts;
end

fprintf(EEG(i_block).ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG(i_block).ALSUTRECHT.subject.fid,'CMS dropout detected (blue light flashing) \n');
fprintf(EEG(i_block).ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');

if any(flagDropout)
    fprintf(EEG(i_block).ALSUTRECHT.subject.fid,'Drop out detected:    Yes\n');
    fprintf(EEG(i_block).ALSUTRECHT.subject.fid,'Number of detections: %d\n', sum(NDropouts));
else
    fprintf(EEG(i_block).ALSUTRECHT.subject.fid,'Drop out detected: No\n');
end

% Report visually
report_badsegments(EEG, jumpsBadAll, 'cmsdropouts', cfg.figure.visible);

% Remove the chunks
if any(flagDropout)
    fprintf('Removing the drop-outs now...\n');
    for i_block = 1:NBLK
        if ~isempty(jumpsBadAll{i_block})
            EEG(i_block) = eeg_eegrej(EEG(i_block), jumpsBadAll{i_block});
            % EEG(i) = pop_select(EEG(i), 'rmpoint', jumpsBadAll{i});
        end
    end
end
end