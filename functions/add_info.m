function EEG = add_info(EEG,subject,cfg)

fprintf('\n================================\n');
fprintf('Adding participant/channel info\n');
fprintf('================================\n');

% Channel locations
chanlocs = readlocs('biosemi128_eeglab.ced');
EEG = fix_chanlocs(EEG,chanlocs);

% Add subject info
NBLK = length(EEG);
for i = 1:NBLK
    EEG(i).ALSUTRECHT.subject = subject;
    EEG(i).ALSUTRECHT.cfg = cfg;
end

end