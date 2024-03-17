function EEG = epoch_rsdata(EEG, epochLength)
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

NBLK = length(EEG.ALSUTRECHT.blockinfo.eo_mask);
goodIndx = cell(NBLK,1);
NTRL = NaN(NBLK,1);

% Check
assert(size(EEG.data,2)==length(EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections));

% Chunk "good" data only, block per block
for i = 1:NBLK
    assert(length(EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections)==length(EEG.ALSUTRECHT.blockinfo.rs_mask(i,:)));
    maskGood = ~EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections(EEG.ALSUTRECHT.blockinfo.rs_mask(i,:));

    jump   = [0 find(diff(maskGood)~=0) length(maskGood)];
    jindex = find(diff(jump)>epochLength);

    longjumpS = jump(jindex)+1;
    longjumpE = jump(jindex+1);

    longjumpSgood = longjumpS(maskGood(longjumpS)==1);
    longjumpEgood = longjumpE(maskGood(longjumpE)==1);

    % nEp = floor((longjumpEgood-longjumpSgood+1)/EpochLength);

    EpochsS = cell2mat(arrayfun(@(x,y) (x:epochLength:(max(x,y-epochLength))),longjumpSgood.',longjumpEgood.','UniformOutput',0).');
    goodIndx{i} = cell2mat(arrayfun(@(x) (x:(x+epochLength-1)).',EpochsS,'UniformOutput',0)).';
    NTRL(i) = size(goodIndx{i},1);
end

% Build a minimal but valid EEG.event from scratch
% We ignore previous events (boundry only?)
% as we account for them in the code above!
EEG.event = [];

cnt = 0;
for i = 1:NBLK
    K = find(EEG.ALSUTRECHT.blockinfo.rs_mask(i,:),1,"first")-1;
    for j = 1:NTRL(i)
        cnt = cnt+1;
        if EEG.ALSUTRECHT.blockinfo.eo_mask(i)
            EEG.event(cnt).type = ['EO' num2str(i)];
        else
            EEG.event(cnt).type = ['EC' num2str(i-sum(EEG.ALSUTRECHT.blockinfo.eo_mask))];
        end

        % EEG.event(cnt).latency = (K+goodIndx{i}(j,1)-1)./EEG.srate;
        EEG.event(cnt).latency = K+goodIndx{i}(j,1); % MUST BE IN SAMPLES!
    end
end
EEG = eeg_checkset(EEG);

% Epoch now
allLabels = unique({EEG.event(:).type});
EEG = pop_epoch(EEG,allLabels,[0 epochLength]./EEG.srate,'epochinfo','yes');

% EEG0 = pop_epoch(EEG,allLabels,[0 epochLength]./EEG.srate,'epochinfo','yes');
% figure; hold on;
% plot(EEG.data(1,1:6*256));
% plot(reshape(EEG0.data(1,:,1:3),1,[]));
% axis tight; legend('Continuous','Epoched');

% Check
EEG = eeg_checkset(EEG);
assert(sum(NTRL)==size(EEG.data,3));
assert(epochLength==size(EEG.data,2));
% vis_artifacts(EEG,EEG);

% % -> extreme noise + epoching
% oldL = length(EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections);
% % -> epoching only
oldL = sum(~EEG.ALSUTRECHT.extremeNoise.allMethodsExtremeEpochRejections);
newL = prod(size(EEG.data,[2 3]));
R = 1-newL./oldL;
fprintf('By epoching (L = %ds), %1.2f of the data is lost.\n',epochLength./EEG.srate,R);

% Log
EEG.ALSUTRECHT.blockinfo.goodIndx = goodIndx;
EEG.ALSUTRECHT.blockinfo.goodTrls = NTRL;
EEG.ALSUTRECHT.blockinfo.dataLostbyEpoching = R;

end