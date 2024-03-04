chaneeg = {EEGICA.chanlocs(strcmp({EEGICA.chanlocs.type},'EEG')).labels};
EEGICA2 = pop_select(EEGICA,'channel',chaneeg);

EEGICA2.chanlocs = EEGICA2.chanlocs(:)';
EEGICA2.allchan  = EEGICA2.chanlocs;

chanlabels = {EEGICA2.chanlocs.labels};

RELAX_cfg = [];
% 0 = data almost certainly has blinks, 1 = data might not have blinks, 2 = data definitely doesn't have blinks.
RELAX_cfg.ProbabilityDataHasNoBlinks = 0;
RELAX_cfg.BlinkElectrodes            = chanlabels([72:74, 76:84, 89:96])';
RELAX_cfg.HighPassFilter             = 0.25; % Sets the high pass filter. 1Hz is best for ICA decomposition if you're examining just oscillatory data, 0.25Hz seems to be the highest before ERPs are adversely affected by filtering
RELAX_cfg.LowPassFilter              = 80; % If you filter out data below 75Hz, you can't use the objective muscle detection method
RELAX_cfg.ms_per_sample              = 1000/EEGICA2.srate; 

% Epoch data, detect extremely bad data, delete channels if over the set threshold for proportion of data affected by extreme outlier for each electrode
[continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching(EEGICA2,RELAX_cfg);

% Mark extreme periods for exclusion from MWF cleaning, and deletion before wICA cleaning
[continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg);

% Use the continuous data to detect eye blinks and mark
% these in the EEG.event as well as in the mask. The output is
% continuous data but includes all the previous extreme period
% markings from the epoched data.
if RELAX_cfg.ProbabilityDataHasNoBlinks<2
    [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg); % use an IQR threshold method to detect and mark blinks
    if continuousEEG.RELAX.IQRmethodDetectedBlinks(1,1)==0 % If a participants doesn't show any blinks, make a note
        NoBlinksDetected{FileNumber,1}=FileName;
        warning('No blinks were detected - if blinks are expected then you should visually inspect the file');
    end
end

% Use epoched data and FFT to detect slope of log frequency log
% power, add periods exceeding muscle threshold to mask:
[continuousEEG, epochedEEG] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg);

EEGICA2 = continuousEEG;
% The following pads very brief lengths of mask periods
% in the template (without doing this, very short periods can
% lead to rank deficiency), and excludes extreme artifacts from the
% cleaning template (so the MWF cleaning step just ignores extreme
% artifacts in it's template - doesn't include them in either the
% clean or artifact mask, but does apply cleaning to them).

% If period has been marked as shorter than RELAX_cfg.MinimumArtifactDuration, then pad it out.
[EEGICA2] = RELAX_pad_brief_mask_periods(EEGICA2, RELAX_cfg, 'notblinks');

EEGICA2.RELAX.NoiseMaskFullLengthR1=EEGICA2.RELAXProcessing.Details.NoiseMaskFullLength;
EEGICA2.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal=mean(EEGICA2.RELAXProcessing.Details.NoiseMaskFullLength,'omitnan');
EEGICA2.RELAX.ProportionMarkedInMWFArtifactMaskTotalR1=EEGICA2.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal;

% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
EEGICA22 = RELAX_perform_MWF_cleaning(EEGICA2, RELAX_cfg);


vis_artifacts(EEGICA22,EEGICA2);