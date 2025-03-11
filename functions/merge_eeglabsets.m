function output = merge_eeglabsets(varargin)
% Order:
% EEG (+ EMG) + EXT

fprintf('\n================================\n');
fprintf('Merging signals (EEG/EMG/EXT)\n');
fprintf('================================\n');

output = varargin{1};
for i = 2:numel(varargin)
    output.data = [output.data; varargin{i}.data];

    tmp  = varargin{i}.chanlocs;
    flds = fields(output.chanlocs);
    NFLD = length(flds);

    cnt = length(output.chanlocs);
    for j = 1:length(tmp)
        cnt = cnt+1;
        for k = 1:NFLD
            output.chanlocs(cnt).(flds{k}) = tmp(j).(flds{k});
        end
    end
end

% output.chanlocs = [output.chanlocs; varargin{i}.chanlocs];

output.nbchan = size(output.data,1);
output.chaninfo.removedchans = [];

% Check
output = eeg_checkset(output);

output.icaact = [];
fprintf('Done!\n');

end