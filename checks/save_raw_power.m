function save_raw_power(EEG, cfg)

% Merge blocks
EEG = merge_eeglabblocks(EEG);

% Epoch
EEG = epoch_data(EEG, cfg);

% Extract always the first way of epoching
EEG = EEG{1};

% Re-reference to Common Average Reference
EEG = do_reref(EEG, 'aRegular');

% Identify EEG channels robustly
mask_eeg = strcmpi({EEG.chanlocs.type}, 'EEG');
assert(sum(mask_eeg) == 128, 'Expected 128 EEG channels, found %d.', sum(mask_eeg));

% Extract EEG data: [128 x n_pnts x n_trials]
eeg_data = double(EEG.data(mask_eeg, :, :));
[n_chans, n_pnts, n_trials] = size(eeg_data);

% Demean each epoch across the time dimension (identical to estimate_power)
eeg_data = eeg_data - mean(eeg_data, 2);

% Window parameters matching estimate_power (scalar window defaults to Hamming)
win      = n_pnts;
noverlap = 0;
nfft     = n_pnts;

% Compute trial 1 to obtain frequency vector and preallocate
[pxx_init, freq_pre] = pwelch(eeg_data(:, :, 1)', win, noverlap, nfft, EEG.srate);
n_freqs = length(freq_pre);

psd_pre_trials = zeros(n_chans, n_freqs, n_trials);
psd_pre_trials(:, :, 1) = pxx_init';

for i_epoch = 2:n_trials
    pxx = pwelch(eeg_data(:, :, i_epoch)', win, noverlap, nfft, EEG.srate);
    psd_pre_trials(:, :, i_epoch) = pxx';
end

% Grand average across all trials (pre-cleaning baseline)
psd_pre = mean(psd_pre_trials, 3, 'omitnan');

% Save output with metadata
subject   = EEG.ALSUTRECHT.subject;
data_name = fullfile(subject.data, [subject.filename '_rawpower.mat']);
save(data_name, 'psd_pre', 'freq_pre', 'n_trials', '-v7.3');

fprintf('Raw PSD saved successfully to:\n  %s\n', data_name);

end