function [subjects, NSUB] = select_preproc_participants(myPaths)

if strcmpi(myPaths.group,'AFM')
    % Preprocess only C9 FCOs
    subjects = select_fcos(myPaths);
    % subjects = setdiff(subjects,list_participants('E:\3_PREPROCESSED_DATA\SART\AFM\T1',{}));

    % % Manually select specific participants
    % subjects = {'ALS27315'};

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