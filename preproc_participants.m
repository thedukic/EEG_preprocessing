function myPathsOut = preproc_participants(i_group, i_visit, myPaths)

% #########################################################################
% Define (don't change)
% #########################################################################
myPathsOut       = myPaths;
myPathsOut.task  = myPaths.task;
myPathsOut.group = myPaths.group{i_group};
myPathsOut.visit = myPaths.visit(i_visit);

% % Messy but it could be imporved if all cohorts are in one folder
% if ismember(myPathsTmp.group, {'ALS','PLS','PMA','MND'})
%     myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata, 'ALS', ['T' num2str(myPathsTmp.visit)]);
% elseif ismember(myPathsTmp.group, {'AFM'})
%     myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata, 'AFM', ['T' num2str(myPathsTmp.visit)]);
% elseif ismember(myPathsTmp.group, {'CONTROL'})
%     myPathsTmp.rawdata  = fullfile(myPathsTmp.rootrawdata, 'CONTROL', ['T' num2str(myPathsTmp.visit)]);
% end
myPathsOut.rawdata  = fullfile(myPathsOut.rootrawdata, myPathsOut.group, ['T' num2str(myPathsOut.visit)]);
myPathsOut.preproc  = fullfile(myPathsOut.rootpreproc, myPathsOut.task, myPathsOut.group, ['T' num2str(myPathsOut.visit)]);

% #########################################################################
% Select a method for participant selection (change if needed)
% #########################################################################
% -------------------------------------------------------------------------
% Method A: Preprocess all participants in the given folder
% -------------------------------------------------------------------------
% myPathsOut.subjects = list_participants(myPathsOut.rawdata, {});

% -------------------------------------------------------------------------
% Method B: Select only the relevant participants (Utrecht datasets)
% eg. Folder may have more participants but you want C9+ patients only
% ------------------------------------------------------------------------
myPathsOut.subjects = list_participants(myPathsOut.rawdata, {});
myPathsOut.subjects = select_relevant(myPathsOut.subjects, myPathsOut);

% -------------------------------------------------------------------------
% Method C: Manually select
% -------------------------------------------------------------------------
% myPathsOut.subjects = {'ALS12345'};
% load(fullfile(myPaths.mycodes, 'files', 'list_c9_als.mat'), 'list_als'); myPathsOut.subjects = list_als;

% #########################################################################
% % Check (but fails for DUB data)
% assert(all(contains(subjects, 'ALS')));

% Report
num_subjects = length(myPathsOut.subjects);
fprintf('Processing %d %s participants...\n', num_subjects, myPathsOut.group);

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
    Table1_subsample = Table1(Table1.DIAGNOSIS == 'ALS' & Table1.SUBTYPE == 'ALS', :);
    % Table1_subsample = Table1(Table1.DIAGNOSIS == 'ALS', :);

elseif strcmpi(myPaths.subgroup, 'CONTROL')
    Table1_subsample = Table1(Table1.GROUP == 'CONTROL', :);

elseif strcmpi(myPaths.subgroup, 'AFM_C9ORF72')
    Table1_subsample = Table1(Table1.GROUP == 'AFM' & Table1.GENEPED  == 'C9orf72' & Table1.C9ORF72 ~= 'NA', :);

elseif strcmpi(myPaths.subgroup, 'AFM_ARPP21')
    % load("C:\DATA\MATLAB\myCodes\RS\files\subjects_ARPP21.mat", "subjects");
    % subjects = [{'ALS39019'}, subjects]; % ALS39019 was an fco
    Table1_subsample = Table1(Table1.GROUP == 'AFM' & Table1.GENEPED  == 'ARPP21', :);

elseif strcmpi(myPaths.subgroup, 'MND_C9ORF72')
    Table1_subsample = Table1(Table1.GROUP == 'PATIENT' & (Table1.GENEPED  == 'C9orf72' | Table1.C9ORF72 == 'TRUE'), :);

elseif strcmpi(myPaths.subgroup, 'MND_SOD1')
    load("C:\DATA\MATLAB\myCodes\preprocessing\files\lists\list_sod1_als.mat", "list_sod1_als");
    Table1_subsample = Table1(ismember(Table1.ID, list_sod1_als), :);

else
    error('No participants selected.');
end

% Select visit
Table1_subsample = Table1_subsample(Table1_subsample.VISIT == myPaths.visit, :);

% Report
[mask_missing, mask_extra] = report_table(Table1, Table1_subsample, subjects);

% Remove the ones that are missing
% Recorded data but not in the folder
Table1_subsample(mask_missing, :) = [];
subjects(mask_extra) = [];

end

function [mask_missing, mask_extra] = report_table(Table1, Table1_subsample, subjects)

% Standardise IDs to string arrays
all_ids_full      = string(Table1.ID);
all_ids_subsample = string(Table1_subsample.ID);
subjects_str      = string(subjects);

% =========================================================================
% 1. Missing subjects: In Subsample, passed visit date, but NOT in folder
% =========================================================================
is_past_or_today = Table1_subsample.("DATE: VISIT") <= datetime('today');
mask_missing     = ~ismember(all_ids_subsample, subjects_str) & is_past_or_today;
tbl_missing      = Table1_subsample(mask_missing, :);

Table_Missing = table();
if ~isempty(tbl_missing)
    Table_Missing.STATUS          = repmat("MISSING", height(tbl_missing), 1);
    Table_Missing.ID              = string(tbl_missing.ID);
    Table_Missing.VISIT_DATE      = string(tbl_missing.("DATE: VISIT"));
    Table_Missing.DIAGNOSIS_DATE  = string(tbl_missing.("DATE: DIAGNOSIS"));
    Table_Missing.GROUP           = string(tbl_missing.GROUP);
    Table_Missing.DIAGNOSIS       = string(tbl_missing.DIAGNOSIS);
    Table_Missing.SUBTYPE         = string(tbl_missing.SUBTYPE);
    Table_Missing.CONVERTER       = string(tbl_missing.CONVERTER);
else
    Table_Missing = table('Size', [0, 8], ...
        'VariableTypes', repmat({'string'}, 1, 8), ...
        'VariableNames', {'STATUS', 'ID', 'VISIT_DATE', 'DIAGNOSIS_DATE', 'GROUP', 'DIAGNOSIS', 'SUBTYPE', 'CONVERTER'});
end

% =========================================================================
% 2. Extra subjects: In folder, but NOT in Subsample
% =========================================================================
mask_extra = ~ismember(subjects_str, all_ids_subsample);
extra_ids  = subjects_str(mask_extra);
n_extra    = numel(extra_ids);

Table_Extra = table();
if n_extra > 0
    % Match extra folder IDs against the full master Table1
    [is_in_full, loc_in_full] = ismember(extra_ids, all_ids_full);

    % Initialise string arrays
    v_visit_dates = repmat("NOT_IN_MASTER_TABLE1", n_extra, 1);
    v_diag_dates  = repmat("NOT_IN_MASTER_TABLE1", n_extra, 1);
    v_group       = repmat("NOT_IN_MASTER_TABLE1", n_extra, 1);
    v_diag        = repmat("NOT_IN_MASTER_TABLE1", n_extra, 1);
    v_subt        = repmat("NOT_IN_MASTER_TABLE1", n_extra, 1);
    v_conv        = repmat("NOT_IN_MASTER_TABLE1", n_extra, 1);

    % Fill in metadata for IDs found in master Table1
    matched_idx = loc_in_full(is_in_full);
    v_visit_dates(is_in_full) = string(Table1.("DATE: VISIT")(matched_idx));
    v_diag_dates(is_in_full)  = string(Table1.("DATE: DIAGNOSIS")(matched_idx));
    v_group(is_in_full)       = string(Table1.GROUP(matched_idx));
    v_diag(is_in_full)        = string(Table1.DIAGNOSIS(matched_idx));
    v_subt(is_in_full)        = string(Table1.SUBTYPE(matched_idx));
    v_conv(is_in_full)        = string(Table1.CONVERTER(matched_idx));

    Table_Extra.STATUS          = repmat("EXTRA", n_extra, 1);
    Table_Extra.ID              = extra_ids(:);
    Table_Extra.VISIT_DATE      = v_visit_dates;
    Table_Extra.DIAGNOSIS_DATE  = v_diag_dates;
    Table_Extra.GROUP           = v_group;
    Table_Extra.DIAGNOSIS       = v_diag;
    Table_Extra.SUBTYPE         = v_subt;
    Table_Extra.CONVERTER       = v_conv;
else
    Table_Extra = table('Size', [0, 8], ...
        'VariableTypes', repmat({'string'}, 1, 8), ...
        'VariableNames', {'STATUS', 'ID', 'VISIT_DATE', 'DIAGNOSIS_DATE', 'GROUP', 'DIAGNOSIS', 'SUBTYPE', 'CONVERTER'});
end

% =========================================================================
% 3. Combine and Display
% =========================================================================
SummaryTable = [Table_Missing; Table_Extra];

disp('--- Missing and Extra Subjects Summary ---');
disp(SummaryTable);

end