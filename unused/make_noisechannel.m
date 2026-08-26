function EEG = make_noisechannel(EEG,noisemask0)
% Make a mask for each trial
% eg
% 111 000 000
% 000 111 000
% 000 000 111
N = cellfun(@(x,y) size(x,2),{EEG(:).data});
assert(size(EEG(1).data,1)<=136);

NBLK = length(EEG);
mask_trial = false(NBLK,sum(N));
for j = 1:NBLK
    if j == 1
        mask_trial(j,1:N(1)) = true;
    else
        mask_trial(j,sum(N(1:j-1))+1:sum(N(1:j))) = true;
    end
end

NCHN = EEG(1).nbchan+1;
% noisemask = cell(1,NBLK);
for j = 1:NBLK
    % noisemask{j} = [noisemask0(mask_trial(j,:))];
    EEG(j).data(NCHN,:) = [noisemask0(mask_trial(j,:))];
    EEG(j).nbchan = NCHN;
    EEG(j).chanlocs(NCHN).labels = 'NOISE';
    EEG(j).chanlocs(NCHN).type   = 'MSK';
end

% Check
EEG = eeg_checkset(EEG);

end