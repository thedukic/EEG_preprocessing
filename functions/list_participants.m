function subjects = list_participants(folderpath,todolist)
%
% List subjects from the given folder
% SDukic, March 2024
%

FilesList0  = dir(folderpath);
FilesList0  = FilesList0([FilesList0(:).isdir]==1);
[~,sortID]  = sort({FilesList0.name});
FilesList00 = {FilesList0.name};
subjects = FilesList00(1,sortID);
subjects(1:2) = [];

% Maybe we dont want to do them all!
if ~isempty(todolist)
    if islogical(todolist)
        assert(length(todolist)==length(subjects));
        subjects = subjects(todolist);
    elseif isnumeric(todolist)
        subjects = subjects(todolist);
    elseif iscell(todolist)
        if sum(ismember(subjects,todolist))==length(todolist)
            subjects = todolist;
        else
            error('Some of your to-do participants are not in the given folder.');
        end
    else
        error('Check your to-do participant list.');
    end
end

end