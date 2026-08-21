function EEG = fix_events(EEG)
% Added cos of DUB EO data
% -> was not needed for Utrecht data

fprintf('\n================================\n');
fprintf('Fixing events (numeric -> char)\n');
fprintf('================================\n');

NBLK = length(EEG);
for i_block = 1:NBLK
    for i_event = 1:numel(EEG(i_block).event)
        if isnumeric(EEG(i_block).event(i_event).type)
            fprintf('Fixed: %d -> ', EEG(i_block).event(i_event).type);

            EEG(i_block).event(i_event).type   = num2str(EEG(i_block).event(i_event).type);
            EEG(i_block).urevent(i_event).type = EEG(i_block).event(i_event).type;

            fprintf('%s\n', EEG(i_block).event(i_event).type);
        end
    end
end

% Double-check
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG);
fprintf('Done!\n');

end