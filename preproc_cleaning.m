function issues_to_check = preproc_cleaning(myfolders,id)
%
% Script for EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, April 2024
% =========================================================================
% TODO
% 1. Parameter for trial length,
%    maybe good to save data with 2s and 5s or 10s
%

% Load preprocessing settings
cfg = preproc_parameters;

% Define paths and files
subject         = [];
subject.id      = id;
subject.task    = myfolders.task;
subject.group   = myfolders.group;
subject.visit   = myfolders.visit;
subject.rawdata = fullfile(myfolders.rawdata, subject.id);
subject.preproc = fullfile(myfolders.preproc, subject.id);
subject.icadata = fullfile(subject.preproc, upper(cfg.ica.type2));
subject.clnfile = [subject.id '_' myfolders.visit '_' myfolders.task '_cleandata_' cfg.rnum '.mat'];

% Find datasets
[subject.datablocks, NBLK] = list_datasets(subject.rawdata,myfolders.task);

% Report
fprintf('\n');
disp('==================================================================');
fprintf('%s | %s | %s dataset\n',myfolders.group,subject.id,myfolders.task);
disp('==================================================================');
fprintf('\n');

% Load data
if ~isempty(subject.datablocks)
    if strcmpi(myfolders.task,'MT')
        EEG = pop_biosig(subject.datablocks,'channels',1:168);
        EEG = pop_select(EEG,'rmchannel',[131 132 143:160]);
    else
        EEG = pop_biosig(subject.datablocks,'channels',1:136);

        % Rare cases of RS datasets recorded with 168 channels
        if any([EEG(:).srate] == 2048)
            maskfs2k = [EEG(:).srate] == 2048;
            EEG(maskfs2k) = pop_biosig(subject.datablocks(maskfs2k),'channels',1:168);
            EEG(maskfs2k) = pop_select(EEG(maskfs2k),'channel',[1:128 161:168]);
        end
    end
else
    warning([subject.id ' is missing ' myfolders.task ' data!']);
    EEG.ALSUTRECHT.subject.id   = subject.id;
    EEG.ALSUTRECHT.subject.task = subject.task;
    EEG = report_issues(EEG);
    issues_to_check = EEG.ALSUTRECHT.issues_to_check;
    return;
end

% If loaded, then make a folder
if exist(subject.preproc,'dir')~=7, mkdir(subject.preproc); end

% Open a report
t0 = datetime("now");
proctime = strrep(strrep(char(t0),':','-'),' ','-');
subject.fid = fopen(fullfile(subject.preproc,['report_preprocess_v' cfg.rnum '.txt']),'w+');

fprintf(subject.fid,'%s | %s | %s dataset\n\n',myfolders.group,subject.id,myfolders.task);
fprintf(subject.fid,'Code version %s\n',cfg.rnum);
fprintf(subject.fid,'Started: %s\n',t0);

% Manually fix datasets in some rare cases
if strcmp(subject.id,'C50') && strcmpi(myfolders.task,'SART')
    % C48 is very strage - low quality data?
    warning([subject.id 'has swapped C- and B- set. Fixing that now...']);
    for j = 1:NBLK
        EEG(j).data(33:96,:) = [EEG(j).data(65:96,:); EEG(j).data(33:64,:)];
    end
elseif strcmp(subject.id,'ALS26603') && strcmpi(myfolders.task,'MMN')
    warning([subject.id ' has MMN2 file currpted. Removing it now...']);
    EEG(2) = []; NBLK = length(EEG);

elseif strcmp(subject.id,'ALS26603') && strcmpi(myfolders.task,'SART')
    warning([subject.id ' has SART1 and SART3 files currpted. Removing some parts now...']);
    EEG(1) = pop_select(EEG(1),'rmtime',[0 5; 93.3 100.1]);
    EEG(3) = pop_select(EEG(3),'rmtime',[183 201]);
end
EEG = eeg_checkset(EEG,'loaddata');

% Channel locations
chanlocs = readlocs('biosemi128_eeglab.ced');
EEG = fix_chanlocs(EEG,chanlocs);

% Add subject info
for j = 1:NBLK
    EEG(j).ALSUTRECHT.subject = subject;
end

% Demean
EEG = pop_rmbase(EEG,[]);

% Resample to 256 Hz
EEG = pop_resample(EEG,256);

% Keep event info
EEG = extract_eventinfo(EEG);

% Detect flat channels on raw data
EEG = remove_flatelec(EEG,cfg.bch);

% Reference
EEG = do_reref(EEG,'aRobust');

% Remove the line noise here?
% EEG = reduce_linenoise(EEG);

% Filter
EEG = do_filtering(EEG,cfg);

% Remove line noise
EEG = reduce_linenoise(EEG);

% Check leftovers
EEG = detect_linenoiseleftovers(EEG);

% Cut block ends
EEG = remove_datasetends(EEG);

% Bipolar EXT
EEG = make_extbipolar(EEG);

% Mark where each RS block starts/ends
if strcmpi(myfolders.task,'RS') || strcmpi(myfolders.task,'EO') || strcmpi(myfolders.task,'EC')
    EEG = make_rsmasks(EEG);
end

% Merge datasets
EEG = pop_mergeset(EEG,1:NBLK);

% Clean EMG
if strcmpi(myfolders.task,'MT')
    % EEG = clean_emgdata(EEG);
end

% Separate EMG before EEG cleaning
if strcmpi(myfolders.task,'MT')
    chanemg = {EEG(1).chanlocs(strcmp({EEG(1).chanlocs.type},'EMG')).labels};
    EMG = pop_select(EEG,'channel',chanemg);
    EEG = pop_select(EEG,'rmchannel',chanemg);
end

% Make a copy
% EEGRAW = EEG;

% Remove electrode wobbles and pops
EEG = remove_electrodewobbles(EEG,lower(cfg.ica.type1));

% Detect and remove noisy electrodes
EEG = remove_noisyelec(EEG,cfg.bch);

% Log/report bad channel
EEG = report_badelectrodes(EEG);

% Detect extremely bad epochs
EEG = detect_extremelybadepochs2(EEG);

% Remove extremely bad epochs
if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)
    EEG = eeg_eegrej(EEG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);

    if strcmpi(myfolders.task,'MT')
        EMG = eeg_eegrej(EMG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    end
end

% % Check if EC has eye blinks; Too sensitive?
% if strcmpi(myfolders.task,'EC')
%     EEG = check_eyesclosedeyeblinks(EEG);
% end

% Do MWF and ASR
EEG = reduce_artifacts(EEG,cfg.bch);

% Separate EXT channels before ICA
chanext = {EEG.chanlocs(strcmp({EEG.chanlocs.type},'EXT')).labels};
EXT = pop_select(EEG,'channel',chanext);
EEG = pop_select(EEG,'rmchannel',chanext);

% Interpolate bad electrodes
if ~isempty(EEG.ALSUTRECHT.badchaninfo.badElectrodes)
    EEG = pop_interp(EEG,chanlocs);
end

% Common-average before ICA
EEG = do_reref(EEG,'aRegular');

% Epoch MT data
if strcmpi(myfolders.task,'MT')
    % Epoch MT data before ICA, as there is often a lot of noise in between trials
    EEG = pop_epoch(EEG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EXT = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EMG = pop_epoch(EMG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0] ,[]);
end

% ICA
EEG = do_ICA(EEG,cfg);

% Flag artifact ICs
EEG = detect_badICs(EEG,EXT,cfg);

% ICA log/report
EEG = report_ICA(EEG);

% Epoch MMN
if strcmpi(myfolders.task,'MMN')
    EEG = pop_epoch(EEG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0),cfg.trg.mmn{2},'epochinfo','yes');
    EXT = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0),cfg.trg.mmn{2},'epochinfo','yes');
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0] ,[]);
end
% Epoch SART: 1. wrt visual stimuli and 2. wrt response times
if strcmpi(myfolders.task,'SART')
    EEG0 = EEG; EXT0 = EXT;
    EEG  = pop_epoch(EEG0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0),cfg.trg.sart1{2},'epochinfo','yes');
    EXT  = pop_epoch(EXT0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0),cfg.trg.sart1{2},'epochinfo','yes');
    EEG  = pop_rmbase(EEG,[(EEG.xmin)*1000 0] ,[]);

    % This does not make sure that epochs are only from correct Go
    EEG2 = pop_epoch(EEG0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0),cfg.trg.sart2{2},'epochinfo','yes');
    EXT2 = pop_epoch(EXT0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0),cfg.trg.sart2{2},'epochinfo','yes');
    EEG2 = pop_rmbase(EEG2,[(EEG2.xmin)*1000 0] ,[]);
    clearvars EEG0 EXT0
elseif strcmpi(myfolders.task,'RS') || strcmpi(myfolders.task,'EO') || strcmpi(myfolders.task,'EC')
    % 2s w/ 0.75 overlap
    EEG = epoch_rsdata3(EEG,cfg.trg.rs{1},cfg.trg.rs{2});
    EXT = epoch_rsdata3(EXT,cfg.trg.rs{1},cfg.trg.rs{2});
end

% % ICLabel
% EEG = iclabel(EEG);
% if strcmpi(myfolders.task,'SART'), EEG2 = iclabel(EEG2); end
%
% % Flag artifact ICs
% EEG = pop_icflag(EEG,cfg.ica.iclabel);
% if strcmpi(myfolders.task,'SART'), EEG2 = pop_icflag(EEG2,cfg.ica.iclabel); end
%
% % ICA log/report
% EEG = report_ica(EEG);
% if strcmpi(myfolders.task,'SART'), EEG2 = report_ica(EEG2); end

% Interim data saving
% EEGM = pop_saveset(EEGM,'filename',[subject.icafile '.set'],'filepath',subject.preproc);

% Merge EXT
EEG = merge_eeglabsets(EEG,EXT);
if strcmpi(myfolders.task,'SART'), EEG2 = merge_eeglabsets(EEG2,EXT2); end

% Do wICA
EEG = do_wICA(EEG);
if strcmpi(myfolders.task,'SART'), EEG2 = do_wICA(EEG2); end

% Visually check the cleaning
% compare_visually(EEG,EEGRAW,cfg.trg);

% Report artifact leftovers
EEG = report_leftovers(EEG);

% Merge EMG
if strcmpi(myfolders.task,'MT')
    EEG = merge_eeglabsets(EEG,EMG);
end

% Record warnings about potential issues
EEG = report_issues(EEG);
issues_to_check = EEG.ALSUTRECHT.issues_to_check;
if strcmpi(myfolders.task,'SART'), EEG2.ALSUTRECHT.issues_to_check = issues_to_check; end

% Save cleaned data
fprintf('\n%s: Saving preprocessed data...\n',subject.id);
preprocReport = EEG.ALSUTRECHT;
if ~strcmpi(myfolders.task,'SART')
    save(fullfile(subject.preproc,subject.clnfile),'EEG','preprocReport','cfg','proctime');
else
    save(fullfile(subject.preproc,subject.clnfile),'EEG','EEG2','preprocReport','cfg','proctime');
end

% Close the report
t1 = datetime("now");
dd = round(minutes(diff([t0 t1])));
fprintf(subject.fid,'\n\nFinished: %s\n',t1);
fprintf(subject.fid,'Running time: %d min.\n',dd);
fclose(subject.fid);

fprintf('Finished: %s\n',t1);
fprintf('Running time: %d min.\n\n',dd);
end