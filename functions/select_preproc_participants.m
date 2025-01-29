function [subjects, NSUB] = select_preproc_participants(myPaths)

if strcmpi(myPaths.group,'AFM')
    % Preprocess only C9 FCOs
    subjects = select_fcos(myPaths);

    % % Manually select specific participants
    % subjects = {'ALS27315'};

else
    % Preprocess all participants
    subjects = list_participants(myPaths.rawdata,{});
end

% Check
assert(all(contains(subjects,'ALS')));

% Report
NSUB = length(subjects);
fprintf('Processing %d %s participants....\n',NSUB,myPaths.group);

end