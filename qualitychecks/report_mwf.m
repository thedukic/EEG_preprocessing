function [EEG, d] = report_mwf(EEG,EEGRAW)
% This function uses all artifact templates from the MWF cleaning steps
% to compute the SER and ARR cleaning efficacy metrics

% % Remove the same periods in the raw EEG as have been removed
% EEGRAW = eeg_eegrej(EEGRAW, EEGRAW.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);

% Combine all MWF masks into one
NoiseMaskFullLengthAll = EEG.ALSUTRECHT.MWF.R1.noiseMask + EEG.ALSUTRECHT.MWF.R2.noiseMask;
NoiseMaskFullLengthAll(NoiseMaskFullLengthAll>1) = 1;

% Average re-reference
EEG = do_reref(EEG,'aRegular');
EEGRAW = do_reref(EEGRAW,'aRegular');

% Calculate the artifact that has been removed
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
d = EEGRAW.data(chaneeg,:)-EEG.data(chaneeg,:);

% Calculate SER and ARR for each type of artifact in the MWF masks:
[EEG.ALSUTRECHT.MWF.cleaned.SER, EEG.ALSUTRECHT.MWF.cleaned.ARR] = mwf_performance(EEGRAW.data(chaneeg,:), d, NoiseMaskFullLengthAll);

% Visual check
% vis_artifacts(EEG,EEGRAW);

end