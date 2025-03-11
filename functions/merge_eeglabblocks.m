function EEG = merge_eeglabblocks(EEG)
% Merge all blocks of data

fprintf('\n================================\n');
fprintf('Merging blocks of data\n');
fprintf('================================\n');

EEG = pop_mergeset(EEG,1:length(EEG));
fprintf('Done!\n');

end