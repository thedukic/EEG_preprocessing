function EEG = do_resampling(EEG, fs_new)

fprintf('\n================================\n');
fprintf('Resampling\n');
fprintf('================================\n');

% Check if resampling is actually necessary
if EEG(1).srate == fs_new
    fprintf('Data is already sampled at %d Hz. Skipping.\n', fs_new); return;
end

% Report
fprintf('Resampling %d channels from %d Hz to %d Hz.\n', EEG(1).nbchan, EEG(1).srate, fs_new);

% Resample using EEGLAB (handles anti-aliasing low-pass filter and event latency scaling)
EEG = pop_resample(EEG, fs_new);

% Validate internal EEGLAB structure consistency
EEG = eeg_checkset(EEG);

fprintf('Done!\n');

end