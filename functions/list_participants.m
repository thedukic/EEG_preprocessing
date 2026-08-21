function subjects = list_participants(path_folder,list_todo)
%
% List subjects from the given folder
% SDukic, June 2026
%

FilesList0  = dir(path_folder);
FilesList0  = FilesList0([FilesList0(:).isdir] == 1);
[~, sortID] = sort({FilesList0.name});
FilesList00 = {FilesList0.name};

if ~isempty(FilesList00)
    subjects = FilesList00(1, sortID);
    subjects(1:2) = [];

    % Maybe we dont want to do them all
    if ~isempty(list_todo)
        if islogical(list_todo)
            assert(length(list_todo) == length(subjects));
            subjects = subjects(list_todo);

        elseif isnumeric(list_todo)
            subjects = subjects(list_todo);

        elseif iscell(list_todo)
            todomask = ismember(subjects,list_todo);
            if sum(todomask) == length(list_todo)
                subjects = list_todo;
            else
                warning('Some of your to-do participants are not in the given folder.');
                subjects = subjects(todomask);
            end

        else
            error('Check your to-do participant list.');
        end
    end
else
    warning('Folder empty or not found.');
    subjects = {};
end

% mask = cellfun(@(x) contains(x,'ALS'), subjects);
% subjects = subjects(mask);

end