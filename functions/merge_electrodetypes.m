function output = merge_electrodetypes(varargin)
% Order: EEG (+ EMG) + EXT

fprintf('\n================================\n');
fprintf('Merging electrode types (EEG / EMG / EXT)\n');
fprintf('================================\n');

output = varargin{1};
for i_data = 2:numel(varargin)
    if ~isempty(varargin{i_data})
        % Append data
        output.data = [output.data; varargin{i_data}.data];

        % Append chanlocs info
        chanlocs_tmp = varargin{i_data}.chanlocs;
        chanlocs_fields = fields(output.chanlocs);

        cnt = length(output.chanlocs);
        for i_chan = 1:length(chanlocs_tmp)
            cnt = cnt+1;
            for i_field = 1:length(chanlocs_fields)
                output.chanlocs(cnt).(chanlocs_fields{i_field}) = chanlocs_tmp(i_chan).(chanlocs_fields{i_field});
            end
        end
    end
end

% output.chanlocs = [output.chanlocs; varargin{i}.chanlocs];

% Update
output.nbchan = size(output.data,1);
output.chaninfo.removedchans = [];

% Check
output = eeg_checkset(output);
output.icaact = [];

fprintf('Done!\n');

end