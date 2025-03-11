function EEG = remove_trends(EEG)

ORDER = 10;
eegchans = strcmp({EEG(1).chanlocs.type},'EEG');

for i = 1:length(EEG)
    % Detrend
    [x, w] = nt_detrend(EEG(i).data(eegchans,:)', ORDER);

    % Place back
    EEG(i).data(eegchans,:) = x';

    % figure(3); clf;
    % plot(x); title('detrended');
    % figure(4); clf;
    % nt_imagescc(w'); title('weights from nt_detrend','interpreter','none'); ylabel('channels');
end

end