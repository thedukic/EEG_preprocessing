function EEG = do_channelinterp(EEG,thisType)

fprintf('\n================================\n');
fprintf('Channel interpolation (%s)\n',thisType);
fprintf('================================\n');

if ~isempty(EEG.ALSUTRECHT.badchaninfo.badElectrodes)
    chanlocs = readlocs('biosemi128_eeglab.ced');
    EEG = pop_interp(EEG,chanlocs,'spherical');
end

end