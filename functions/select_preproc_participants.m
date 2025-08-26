function [subjects, NSUB] = select_preproc_participants(myPaths)

if strcmpi(myPaths.group,'AFM')
    % Preprocess only C9 FCOs
    subjects = select_fcos(myPaths);
    % subjects = setdiff(subjects,list_participants('E:\3_PREPROCESSED_DATA\SART\AFM\T1',{}));

    % % Manually select specific participants
    % load("C:\DATA\MATLAB\myCodes\RS\files\subjects_ARPP21.mat","subjects");
    % subjects = [{'ALS39019'}, subjects];

    % % Manually select T2 C9 AFM (which are used in T1)
    % load('C:\DATA\MATLAB\myCodes\RS\files\subjectsAFM89.mat', 'subjectsAFM89');
    % pathBulk = 'Y:\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\1_RAW\EEG_DATA\AFM\T2';
    % subjects = list_participants(pathBulk,subjectsAFM89);
    % subjects = {'ALS34276','ALS34475','ALS38801'};

    % subjectsSSD = list_participants('E:\1_EEG_DATA\AFM\T2',[]);
    % m = ismember(subjects,subjectsSSD);
    % subjects(~m)
else
    % Preprocess all participants
    subjects = list_participants(myPaths.rawdata,{});
end

% Check (fails for DUB data)
% assert(all(contains(subjects,'ALS')));

% Report
NSUB = length(subjects);
fprintf('Processing %d %s participants...\n',NSUB,myPaths.group);

end