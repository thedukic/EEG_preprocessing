function EEG = do_channelinterp(EEG,thisType)

fprintf('\n================================\n');
fprintf('Interpolating EEG electrodes (%s)\n', thisType);
fprintf('================================\n');

% % if ~isempty(EEG.ALSUTRECHT.badchaninfo.badElectrodes)
% NBLK = length(EEG);
% for i = 1:NBLK
%     chanlocs = readlocs('biosemi128_eeglab.ced');
%     EEG(i) = pop_interp(EEG(i), chanlocs, 'spherical');
% end
% % end

if EEG.nbchan < 128
    chanlocs = readlocs('biosemi128_eeglab.ced');
    EEG = pop_interp(EEG, chanlocs, 'spherical');
else
    fprintf('Skipping! All electrodes are already in the data.\n');
end

end