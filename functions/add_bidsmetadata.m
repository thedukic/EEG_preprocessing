function EEG = add_bidsmetadata(EEG, subject)
% ADD_BIDS_METADATA Updates the EEG structure with BIDS-compliant naming.
%
% Usage: EEG = add_bids_metadata(EEG, subject)

% 1. Extract core identifiers
subID     = subject.id;
groupName = subject.group;
taskName  = subject.task;
descTag   = 'epochs';

% 2. Standardise Session ID
% %02d automatically adds a leading zero for numbers 1-9
% and leaves 10+ alone. No length() checks required.
sesID = sprintf('%02s', subject.visit(2:end));

% 3. Generate BIDS Filename
% sub-<label>_ses-<label>_task-<label>_desc-<label>_eeg.set
fileName = sprintf('sub-%s_ses-%s_task-%s_desc-%s_eeg.set', ...
    subID, sesID, taskName, descTag);

% 4. Update EEGLAB metadata fields
EEG.setname = fileName;
EEG.subject = subID;
EEG.session = sesID;
EEG.group   = groupName;

% 5. Quality Status
EEG.etc.quality_status = EEG.ALSUTRECHT.automagicmetrics.status;

end