function todolist = list_missing(pathraw,pathdone)

% pathraw =  'E:\1_EEG_DATA\AFM\T1';
% pathdone = 'E:\3_PREPROCESSED_DATA\RS\AFM\T1';

A = list_subjects(pathdone,[]);
B = list_subjects(pathraw,[]);
todolist = B(~ismember(B,A));


end