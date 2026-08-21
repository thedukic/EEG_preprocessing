function export_data(EEG, subject, i_file)
% EXPORT_DATA Saves pre-processed EEG data and lightweight QA metrics.
%
% This function handles the full dataset save and extracts the
% ALSUTRECHT metadata for faster cohort-level analysis.

% 1. Determine file paths and pre-processing stage
if i_file == 0
    data_name = fullfile(subject.data, subject.clnfile0);
    stage_label = 1;
else
    data_name = fullfile(subject.data, subject.clnfile{i_file});
    stage_label = 2;
end

% 2. Ensure the main pre-processing directory exists
if ~exist(subject.data, 'dir')
    mkdir(subject.data);
end

% Update the fields
EEG.ALSUTRECHT.subject = subject;
EEG.icaact = [];

% 3. Save the full EEG dataset
fprintf('\n%s: Saving preprocessed data (Part %d)...\n', subject.id, stage_label);
save(data_name, 'EEG');

% 4. Export lightweight QA metrics (only for post-cleaning stages)
if i_file > 0
    export_qa_metrics(EEG, subject, i_file);
end

end

function export_qa_metrics(EEG, subject, i_file)
% Helper function to keep the main export logic clean.
if ~exist(subject.qa, 'dir')
    mkdir(subject.qa);
end

% Verify metadata exists before saving
qa_name = fullfile(subject.qa, subject.qametrics{i_file});
qa_data = EEG.ALSUTRECHT;

save(qa_name, 'qa_data', '-v7.3');
fprintf('%s: QA metrics exported successfully.\n', subject.id);
end