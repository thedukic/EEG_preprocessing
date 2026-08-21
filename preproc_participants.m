function myPathsOut = preproc_participants(i_group, i_visit, myPaths)

% Define
myPathsOut          = myPaths;
myPathsOut.task     = myPaths.task;
myPathsOut.group    = myPaths.group{i_group};
myPathsOut.visit    = myPaths.visit(i_visit);

% % Messy but it could be imporved if all data is in one folder
% if ismember(myPathsTmp.group , {'ALS','PLS','PMA','MND'})
%     myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata, 'ALS', ['T' num2str(myPathsTmp.visit)]);
% elseif ismember(myPathsTmp.group , {'AFM'})
%     myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata, 'AFM', ['T' num2str(myPathsTmp.visit)]);
% elseif ismember(myPathsTmp.group , {'CONTROL'})
%     myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata, 'CONTROL', ['T' num2str(myPathsTmp.visit)]);
% end
myPathsOut.rawdata  = fullfile(myPathsOut.rootrawdata, myPathsOut.group, ['T' num2str(myPathsOut.visit)]);
myPathsOut.preproc  = fullfile(myPathsOut.rootpreproc, myPathsOut.task, myPathsOut.group, ['T' num2str(myPathsOut.visit)]);

% Preprocess all participants
myPathsOut.subjects = list_participants(myPathsOut.rawdata, {});

% Select only the relevant participants
% eg. folder may have more participants but you want ALS only
myPathsOut.subjects = select_relevant(myPathsOut.subjects, myPathsOut);

% % Overrride
% myPathsOut.subjects = {'ALS36104'};
% load('C:\DATA\MATLAB\myCodes\preprocessing\files\list_c9_als.mat', 'list_als'); myPathsOut.subjects = list_als;

% % Check (but fails for DUB data)
% assert(all(contains(subjects, 'ALS')));

% Report
NSUB = length(myPathsOut.subjects);
fprintf('Processing %d %s participants...\n', NSUB, myPathsOut.group);

end

% =========================================================================
% HELPER FUNCTION
% =========================================================================
function subjects = select_relevant(subjects, myPaths)

% Column
subjects = subjects(:);

% Load
Table1 = readtable(myPaths.table1, 'VariableNamingRule', 'preserve');

% Categorial
Table1.ID        = categorical(Table1.ID);
Table1.GROUP     = categorical(Table1.GROUP);
Table1.DIAGNOSIS = categorical(Table1.DIAGNOSIS);
Table1.SUBTYPE   = categorical(Table1.SUBTYPE);
Table1.C9ORF72   = categorical(Table1.C9ORF72);
Table1.GENEPED   = categorical(Table1.("PED GENE1"));

% Check if a subgroup is needed
if isempty(myPaths.subgroup)
    myPaths.subgroup = myPaths.group;
end
fprintf('Checking group (subgroup): %s (%s)\n', myPaths.group, myPaths.subgroup);

% Select group
if strcmpi(myPaths.subgroup, 'ALS')
    Table1 = Table1(Table1.DIAGNOSIS == 'ALS' & Table1.SUBTYPE == 'ALS', :);

elseif strcmpi(myPaths.subgroup, 'CONTROL')
    Table1 = Table1(Table1.GROUP == 'CONTROL', :);

elseif strcmpi(myPaths.subgroup, 'AFM_C9ORF72')
    Table1 = Table1(Table1.GROUP == 'AFM' & Table1.GENEPED  == 'C9orf72' & Table1.C9ORF72 ~= 'NA', :);

elseif strcmpi(myPaths.subgroup, 'AFM_ARPP21')
    % load("C:\DATA\MATLAB\myCodes\RS\files\subjects_ARPP21.mat", "subjects");
    % subjects = [{'ALS39019'}, subjects]; % ALS39019 was an fco
    Table1 = Table1(Table1.GROUP == 'AFM' & Table1.GENEPED  == 'ARPP21', :);

elseif strcmpi(myPaths.subgroup, 'MND_C9ORF72')
    Table1 = Table1(Table1.GROUP == 'PATIENT' & (Table1.GENEPED  == 'C9orf72' | Table1.C9ORF72 == 'TRUE'), :);

elseif strcmpi(myPaths.subgroup, 'MND_SOD1')
    load("C:\DATA\MATLAB\myCodes\preprocessing\files\lists\list_sod1_als.mat", "list_sod1_als");
    Table1 = Table1(ismember(Table1.ID, list_sod1_als), :);

else
    error('No participants selected.');
end

% Select visit
Table1 = Table1(Table1.VISIT == myPaths.visit, :);

% =========================================================================
list_table1 = cellstr(Table1.ID);
list_notfound = find(~ismember(list_table1, subjects));

if ~isempty(list_notfound)
    for i = 1:length(list_notfound)
        % Extract the date for the current iteration
        visitDate = Table1.("DATE: VISIT")(list_notfound(i));

        % Only display if the date is today or in the past
        if visitDate <= datetime('today')
            fprintf('Not found: %s (%s)\n', list_table1{list_notfound(i)}, visitDate);
        end
    end
end

list_extrafound1 = ~ismember(subjects, list_table1);
if any(list_extrafound1)
    list_extrafound2 = subjects(list_extrafound1);
    fprintf('Extra found: %s\n', list_extrafound2{:});
end

Table1(list_notfound, :) = [];
subjects(list_extrafound1) = [];

end