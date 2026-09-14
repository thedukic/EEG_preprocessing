function subject = preproc_folders_subject(id, myPaths)
% -------------------------------------------------------------------------
% Define paths and files
% -------------------------------------------------------------------------
subject = [];
subject.id        = id;
subject.task      = myPaths.task;
subject.group     = myPaths.group;
subject.visit     = ['T' num2str(myPaths.visit)];
subject.codever   = myPaths.codever;
subject.mycodes   = myPaths.mycodes;

subject.rawdata   = fullfile(myPaths.rawdata, subject.id);
subject.preproc   = fullfile(myPaths.preproc, subject.id, subject.codever);
subject.data      = fullfile(subject.preproc, 'data');
subject.figures   = fullfile(subject.preproc, 'figures');
subject.qa        = fullfile(subject.preproc, 'qa');

% File name: part 1
subject.filename = [subject.id '_' subject.visit '_' subject.task];
subject.filename_clean_1  = [subject.filename '_cleandata_a.mat'];

% File name: part 2
% TODO: make it dynamical, do not assume epoching strategy
subject.filename_clean_2{1} = [subject.filename '_cleandata_b.mat'];
subject.filename_clean_2{2} = [subject.filename '_cleandata_c.mat'];
subject.qametrics{1}        = [subject.filename '_qametrics_b.mat'];
subject.qametrics{2}        = [subject.filename '_qametrics_c.mat'];

end