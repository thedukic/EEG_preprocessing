function subjects = select_fcos(myPaths)
% Preprocess selected C9 AFM participants

% Initialise
drivedata = 'E:';
myPaths.excpath = fullfile(drivedata,'2_OTHER_DATA\Excel\Utrecht\');
myPaths.gendata = [myPaths.excpath 'C9status.xlsx'];
myPaths.peddata = [myPaths.excpath 'Pedigrees.xlsx'];
myPaths.cogdata = [myPaths.excpath 'ECAS.txt'];
myPaths.nexdata = [myPaths.excpath 'NE.txt'];
myPaths.dmdata1 = [myPaths.excpath 'Table1.txt'];

% Add paths temporarily
addpath('C:\DATA\MATLAB\myCodes\RS\common');
addpath('C:\DATA\MATLAB\myCodes\Progeny');
addpath('C:\DATA\MATLAB\myCodes\RS\external\FisherTest');

% Find them
myPaths.excl = {};
subjects = select_participants([],'C9',myPaths);
close all

% Remove them
rmpath('C:\DATA\MATLAB\myCodes\RS\common');
rmpath('C:\DATA\MATLAB\myCodes\Progeny');
rmpath('C:\DATA\MATLAB\myCodes\RS\external\FisherTest');

end