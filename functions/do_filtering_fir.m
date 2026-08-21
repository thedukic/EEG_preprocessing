function DATA = do_filtering_fir(DATA)

thisFiltering = 'highpass';
fprintf('\n================================\n');
fprintf('Filtering data (%s)\n', thisFiltering);
fprintf('================================\n');

% % Separate
% [EEG, EMG, EXT] = separate_electrodetypes(DATA);
% 
% % Filter
% EEG = pop_eegfiltnew(EEG, 1, 0, 1650, 0, [], 0);
% 
% % Merge
% DATA = merge_electrodetypes(EEG, EMG, EXT);

% Filter
DATA = pop_eegfiltnew(DATA, 1, 0, 1650, 0, [], 0);

end