function subjects = list_subjects(folderpath,todolist)
%
% List subjects from the given folder
% v1
% SDukic, February 2024

FilesList0  = dir(folderpath);
FilesList0  = FilesList0([FilesList0(:).isdir]==1);
[~,sortID]  = sort({FilesList0.name});
FilesList00 = {FilesList0.name};
subjects = FilesList00(1,sortID);
subjects(1:2) = [];

% Maybe we dont want them all
if ~isempty(todolist)
    % TODO check if indices/logical vector, otherwise do strcmp
    subjects = subjects(todolist);
end

end