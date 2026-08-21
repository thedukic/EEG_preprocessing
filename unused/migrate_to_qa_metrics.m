function migrate_to_qa_metrics(myPaths, subjects)
% migrate_to_qa_metrics(myPaths, subjects)
%
% A temporary utility to prepare existing preprocessed datasets for the
% refactored dashboard by extracting and saving lightweight QA files.

fprintf('\n==================================================================\n');
fprintf('Migrating existing datasets to lightweight QA format\n');
fprintf('==================================================================\n');

NSUB = length(subjects);

for i = 1:NSUB
    % Define the path to the fully preprocessed 'b' file
    subject = preproc_subject(subjects{i}, myPaths, 2);
    fileName = fullfile(subject.preproc, subject.clnfile{1});

    % Load the EEG struct
    temp = load(fileName, 'EEG');

    % Call the extraction function
    export_data(temp.EEG, subject, 1);
end

fprintf('\n==================================================================\n');
fprintf('Migration complete. You can now run report_final.\n');
fprintf('==================================================================\n');

end