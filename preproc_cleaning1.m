function preproc_cleaning1(myPaths,id)
%
% Script for EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, January 2025
% =========================================================================

% Load preprocessing settings
cfg = preproc_parameters;

% Define paths and files
subject         = [];
subject.id      = id;
subject.task    = myPaths.task;
subject.group   = myPaths.group;
subject.visit   = myPaths.visit;
subject.mycodes = myPaths.mycodes;
subject.rawdata = fullfile(myPaths.rawdata, subject.id);
subject.preproc = fullfile(myPaths.preproc, subject.id);
subject.icadata = fullfile(subject.preproc, upper(cfg.ica.type2));
subject.clnfile = [subject.id '_' myPaths.visit '_' myPaths.task '_cleandata_' myPaths.rnum 'a.mat'];

% Find datasets
[subject.datablocks, NBLK] = list_datasets(subject.rawdata,myPaths.task);

% Report
fprintf('\n');
disp('==================================================================');
fprintf('%s | %s | %s dataset | processing part 1\n',myPaths.group,subject.id,myPaths.task);
disp('==================================================================');
fprintf('\n');

% Load data
if ~isempty(subject.datablocks)
    if strcmpi(myPaths.task,'MT')
        EEG = pop_biosig(subject.datablocks,'channels',1:168);
        EEG = pop_select(EEG,'rmchannel',[131 132 143:160]);
    else
        EEG = pop_biosig(subject.datablocks,'channels',1:136);

        % Rare cases of RS datasets recorded with 168 channels
        if any([EEG(:).srate] == 2048)
            pop_editoptions('option_parallel', 0);
            maskfs2k = [EEG(:).srate] == 2048;
            EEG(maskfs2k) = pop_biosig(subject.datablocks(maskfs2k),'channels',1:168);
            EEG(maskfs2k) = pop_select(EEG(maskfs2k),'channel',[1:128 161:168]);
            EEG(maskfs2k) = pop_resample(EEG(maskfs2k),512);
            pop_editoptions('option_parallel', 1);
            assert(all([EEG(:).srate] == 512));
        end
    end
else
    warning([subject.id ' is missing ' myPaths.task ' data!']);
    return;
end

% If loaded, then make a folder
if exist(subject.preproc,'dir')~=7, mkdir(subject.preproc); end

% Open a report
t0 = datetime("now");
procTimeTags = {myPaths.proctime; strrep(strrep(char(t0),':','-'),' ','-')};
subject.fid  = fopen(fullfile(subject.preproc,['report_preprocess_v' myPaths.rnum '.txt']),'w+');

fprintf(subject.fid,'%s | %s | %s dataset\n\n',myPaths.group,subject.id,myPaths.task);
fprintf(subject.fid,'Code version %s\n',myPaths.rnum);
fprintf(subject.fid,'Started: %s\n',t0);

% Manually fix datasets
% C48 is very strage - low quality data?
if strcmp(subject.id,'C50') && strcmpi(myPaths.task,'SART')
    warning([subject.id 'has swapped C- and B- set. Fixing that now...']);
    for i = 1:NBLK
        EEG(i).data(33:96,:) = [EEG(i).data(65:96,:); EEG(i).data(33:64,:)];
    end
end
EEG = eeg_checkset(EEG,'loaddata');

% Channel locations
chanlocs = readlocs('biosemi128_eeglab.ced');
EEG = fix_chanlocs(EEG,chanlocs);

% Add subject info
for i = 1:NBLK
    EEG(i).ALSUTRECHT.subject = subject;
    EEG(i).ALSUTRECHT.cfg = cfg;
end

% Remove CMS drop-outs (blue light flashing)
EEG = detect_dropouts(EEG);

% Demean
EEG = pop_rmbase(EEG,[]);

% Resample to 256 Hz
EEG = pop_resample(EEG,256);

% Keep event info
EEG = extract_eventinfo(EEG,cfg.trg);

% Detect flat channels
EEG = remove_flatelectrodes(EEG,cfg.bch);

% Reference
EEG = do_reref(EEG,'aRobust');

% Filter highpass only
EEG = do_filtering(EEG,'highpass',cfg.flt);

% Remove line noise
EEG = reduce_linenoise(EEG);

% Remove line noise leftovers
EEG = reduce_linenoiseleftovers(EEG);

% Cut block ends
EEG = remove_datasetends(EEG);

% Make EXT bipolar
EEG = make_extbipolar(EEG);

% Mark where each RS block starts/ends
if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO') || strcmpi(myPaths.task,'EC')
    EEG = make_rsmasks(EEG);
end

% Merge datasets
EEG = pop_mergeset(EEG,1:NBLK);

% Check if EC has eye blinks; Currently too sensitive?
if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EC')
    EEG = check_eyesclosedeyeblinks(EEG);
end

% % Clean EMG
% if strcmpi(myPaths.task,'MT')
%     EEG = clean_emgdata(EEG);
% end

% Separate EMG before EEG cleaning
if strcmpi(myPaths.task,'MT')
    chanemg = {EEG(1).chanlocs(strcmp({EEG(1).chanlocs.type},'EMG')).labels};
    EMG = pop_select(EEG,'channel',chanemg);
    EEG = pop_select(EEG,'rmchannel',chanemg);
end

% Make a copy
% EEGRAW = EEG;

% Remove electrode wobbles and pops
% EEG = remove_electrodewobbles(EEG,lower(cfg.ica.type1));

% Remove noisy electrodes
EEG = remove_noisyelectrodes(EEG,cfg.bch);

% Log/report bad channel
EEG = report_badelectrodes(EEG);

% Remove extremely bad epochs
if ~strcmpi(myPaths.task,'MT')
    EEG = remove_extremeperiods(EEG);
else
    [EEG, EMG] = remove_extremeperiods(EEG,EMG);
end

% Reduce artifacts
EEG = reduce_artifacts(EEG,cfg.bch);

% Separate EXT channels before ICA
chanext = {EEG.chanlocs(strcmp({EEG.chanlocs.type},'EXT')).labels};
EXT = pop_select(EEG,'channel',chanext);
EEG = pop_select(EEG,'rmchannel',chanext);

% Interpolate bad electrodes
if ~isempty(EEG.ALSUTRECHT.badchaninfo.badElectrodes)
    EEG = pop_interp(EEG,chanlocs,'spherical');
end

% Reference
EEG = do_reref(EEG,'aRegular');

% ICA
EEG = do_ICA(EEG,cfg);

% Detect artifact ICs
EEG = detect_badcomponents2(EEG,EXT,cfg);

% Deal with artifact ICs
EEG = remove_badcomponents(EEG,cfg);

% Report artifact leftovers
EEG = report_leftovers(EEG,EXT,1,cfg);

% Report ICA
EEG = report_ICA(EEG);

% Interim data saving
% EEG = pop_saveset(EEG,'filename',[subject.icafile '.set'],'filepath',subject.preproc);

% Visually check the cleaning
% compare_visually(EEG,EEGRAW,cfg.trg);

% Merge channels
EEG = merge_eeglabsets(EEG,EXT);
if strcmpi(myPaths.task,'MT'), EEG = merge_eeglabsets(EEG,EMG); end
clearvars EXT EMG

% Save cleaned data
fprintf('\n%s: Saving the preprocessed data (part 1)...\n',subject.id);
save(fullfile(subject.preproc,subject.clnfile),'EEG','cfg','procTimeTags');

% Close the report
t1 = datetime("now");
dd = round(minutes(diff([t0 t1])));
fprintf(subject.fid,'\n\nFinished: %s\n',t1);
fprintf(subject.fid,'Running time: %d min.\n',dd);
fclose(subject.fid);

% Report
fprintf('Finished: %s\n',t1);
fprintf('Running time: %d min.\n\n',dd);

end