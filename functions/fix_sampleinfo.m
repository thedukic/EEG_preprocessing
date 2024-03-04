function eeg = fix_sampleinfo(eeg)

if isfield(eeg,'trial')
    if length(eeg.trial)==1
        eeg.sampleinfo = [1 size(eeg.trial{1},2)];
    else
        % nsmp = size(eeg.trial{1},2)*ones(length(eeg.trial),1);
        nsmp = cellfun(@(x) size(x,2),eeg.trial)';
        begsample = cat(1,0,cumsum(nsmp(1:end-1))) + 1;
        endsample = begsample + nsmp - 1;
        eeg.sampleinfo = [begsample endsample];
    end

else
    eeg.sampleinfo = [1 size(eeg.avg,2)];
end

end