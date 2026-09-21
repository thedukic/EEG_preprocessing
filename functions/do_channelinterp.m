function EEG = do_channelinterp(EEG,thisType)

fprintf('\n================================\n');
fprintf('Interpolating EEG electrodes (%s)\n', thisType);
fprintf('================================\n');

if EEG.nbchan < 128
    chanlocs = readlocs('biosemi128_eeglab.ced');
    EEG = pop_interp(EEG, chanlocs, 'spherical');
else
    fprintf('Skipping! All electrodes are already in the data.\n');
end

end