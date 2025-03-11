function EEG = preproc2_finalestimates(EEG,cfg)

fprintf('\n================================\n');
fprintf('Final data estimates\n');
fprintf('================================\n');

% Median voltage shift
fprintf('Estimating voltage range\n');
voltageShiftWithinEpoch = median(range(EEG.data,2),3);

% Channel correlation matrix
fprintf('Estimating channel correlation matrix\n');
EEG = check_channelcov(EEG);

% % Empirical frequency boundaries
% fprintf('Estimating empirical frequency boundaries\n');
% EEG = estimate_gedBounds(EEG);

% EMG leftovers
fprintf('Checking EMG leftovers\n');
slopesChannelsxEpochs = detect_emg(EEG,cfg.bch);
slopesChannelsxEpochs(slopesChannelsxEpochs < cfg.bch.muscleSlopeThreshold) = NaN;
slopesChannelsxEpochs = slopesChannelsxEpochs - cfg.bch.muscleSlopeThreshold;
BadEpochs = sum(slopesChannelsxEpochs,1,'omitnan');
muscleLeftover = mean(BadEpochs > 0);
fprintf('EMG leftovers: %1.2f\n',muscleLeftover);

% Log
EEG.ALSUTRECHT.epochRejections.MedianvoltageshiftwithinepochFinal = voltageShiftWithinEpoch;

% EEG.ALSUTRECHT.epochRejections.badTrialMuscleTreshhold2 = badTrialMuscleTreshhold;
EEG.ALSUTRECHT.epochRejections.muscle2 = muscleLeftover;
EEG.ALSUTRECHT.leftovers.muscle2 = EEG.ALSUTRECHT.epochRejections.muscle2;

end