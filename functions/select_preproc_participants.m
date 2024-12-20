function [subjects, NSUB] = select_preproc_participants(myPaths)

if strcmpi(myPaths.group,'AFM')
    % Preprocess selected AFM participants, tho we could actually do them  all
    addpath('C:\DATA\MATLAB\myCodes\RS\common');
    addpath('C:\DATA\MATLAB\myCodes\Progeny');
    addpath('C:\DATA\MATLAB\myCodes\RS\external\FisherTest');
    myPaths.excl = {};
    subjects = select_participants([],'C9',myPaths);
    rmpath('C:\DATA\MATLAB\myCodes\RS\common');
    rmpath('C:\DATA\MATLAB\myCodes\Progeny');
    rmpath('C:\DATA\MATLAB\myCodes\RS\external\FisherTest');

    % Manually select specific participants
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