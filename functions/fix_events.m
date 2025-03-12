function EEG = fix_events(EEG)
% Added beacse of DUB EO data
% -> was not needed for Utrecht data

fprintf('\n================================\n');
fprintf('Fixing events (numeric -> char)\n');
fprintf('================================\n');

NBLK = length(EEG);
for i_blk = 1:NBLK
    for i_evt = 1:numel(EEG(i_blk).event)
        if isnumeric(EEG(i_blk).event(i_evt).type)
            fprintf('Fixed: %d -> ', EEG(i_blk).event(i_evt).type);
            EEG(i_blk).event(i_evt).type   = num2str(EEG(i_blk).event(i_evt).type);
            EEG(i_blk).urevent(i_evt).type = EEG(i_blk).event(i_evt).type;
            fprintf('%s\n', EEG(i_blk).event(i_evt).type);
        end
    end
end

fprintf('Done!\n');
end