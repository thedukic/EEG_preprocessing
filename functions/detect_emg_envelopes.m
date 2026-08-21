function mask_final = detect_emg_envelopes(dataeeg, srate, filter_emg)

% =====================================================================
% Strong EMG
% =====================================================================
% 2. EEG: Temporarily highpass filter
[bh, ah] = butter(4, filter_emg/(srate/2), 'high');
dataeeg = filtfilt(bh, ah, dataeeg')';

% Traces of EMG power
dataeeg = 15 * abs(dataeeg);

% ======================================
% Try to find clusters of higher numbers
% ======================================
P1 = 99;
treshold = prctile(dataeeg(:), P1);
mask_emg_tmp = dataeeg > treshold;

data1 = dataeeg;
data1(~mask_emg_tmp) = 0;
data1(mask_emg_tmp)  = 1;

smoothLengthSample = @(smoothLengthMS) round(srate * (smoothLengthMS/1000));
P = smoothLengthSample(500);
data1 = movmean(data1', P)';

% ======================================
% At least X% of electrodes must be affected
% ======================================
P2 = 95;
dataTmp = data1;
dataTmp(dataTmp == 0) = [];
treshold = prctile(dataTmp, P2);

mask_final = data1 > treshold;
mask_final = 100 * mean(mask_final,1);

% Mask1
% -> Many channels affected together with the EOG
% -> This happens due in large EMG/movement/blink artifacts
P3 = 15; % X% of elec
extremeMaskTmp2 = mask_final >= 2*P3;

% Mask
mask_final(mask_final < P3) = 0;
mask_final(extremeMaskTmp2) = 2*P3;

P = smoothLengthSample(1000);
mask_final = movmean(mask_final, P) > 0;

end