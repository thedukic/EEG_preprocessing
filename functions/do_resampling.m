function EEG = do_resampling(EEG,fsnew)

fprintf('\n================================\n');
fprintf('Resampling (fs = %d Hz)\n',fsnew);
fprintf('================================\n');

% Demean
EEG = pop_rmbase(EEG,[]);

% Resample
EEG = pop_resample(EEG,fsnew);

end