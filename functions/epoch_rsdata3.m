function EEG = epoch_rsdata3(EEG,epochLength,epochOverlap)

fprintf('\nEpoching (L = %ds) resting-state data...\n',epochLength);

thisTask  = EEG.ALSUTRECHT.subject.task;
thisShift = (1-epochOverlap)*epochLength;

EEG = eeg_eegrej(EEG,EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
EEG = eeg_regepochs(EEG,'recurrence',thisShift,'eventtype',thisTask,'extractepochs','off');
EEG = pop_epoch(EEG,{thisTask},epochLength/2*[-1 1]);
EEG = eeg_checkset(EEG);

% check visually
% pop_eegplot(EEG,1,1,1);

% Log
EEG.ALSUTRECHT.blockinfo.goodTrls   = EEG.trials;
EEG.ALSUTRECHT.blockinfo.trlLength  = epochLength;
EEG.ALSUTRECHT.blockinfo.trlOverlap = epochOverlap;
EEG.ALSUTRECHT.blockinfo.dataLostbyEpoching = NaN;

end