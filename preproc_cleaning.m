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
subject.clnfile = [subject.id '_' myfolders.visit '_' myfolders.task '_cleandata_' cfg.rnum '.mat'];

% Find datasets
[subject.datablocks, NBLK] = list_datasets(subject.rawdata,myfolders.task);

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
subject.fid = fopen(fullfile(subject.preproc,['report_preprocess_v' cfg.rnum '.txt']),'w+');

% Log basic info
t0 = datetime("now");
proctime = strrep(strrep(char(t0),':','-'),' ','-');

fprintf('\n');
disp('==================================================================');
fprintf('%s | %s | %s dataset\n',myfolders.group,subject.id,myfolders.task);
disp('==================================================================');
fprintf('\n');

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

% Subject info
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
% Doing this first and then filtering might help line-noise removal
chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
for j = 1:NBLK
    EEG(j).data(chaneeg,:) = EEG(j).data(chaneeg,:) - trimmean(EEG(j).data(chaneeg,:),10);
    EEG(j).ref = 'average';
end

% Remove the line noise here?
% EEG = reduce_linenoise(EEG);

% Filter EEG only
chaneeg = find(strcmp({EEG(1).chanlocs.type},'EEG'));
if strcmpi(myfolders.task,'RS') || strcmpi(myfolders.task,'MT')
    % RS/MT
    EEG = filter_signal(EEG,cfg.flt.rsmt.lp,cfg.flt.rsmt.hp,chaneeg,'eeglab');
else
    % MMN/SART (ERP)
    EEG = filter_signal(EEG,cfg.flt.erp.lp,cfg.flt.erp.hp,chaneeg,'eeglab');
end

% Filter EXT only
chanext = find(strcmp({EEG(1).chanlocs.type},'EXT'));
EEG = filter_signal(EEG,cfg.flt.ext.lp,cfg.flt.ext.hp,chanext,'eeglab');

% Filter MT only
chanemg = find(strcmp({EEG(1).chanlocs.type},'EMG'));
if strcmpi(myfolders.task,'MT')
    EEG = filter_signal(EEG,cfg.flt.emg.lp,cfg.flt.emg.hp,chanemg,'eeglab');
end

% Cut block ends
EEG = remove_datasetends(EEG);

% Remove the line noise
EEG = reduce_linenoise(EEG);

% Bipolar EXT
EEG = make_extbipolar(EEG);

% Mark where each RS block starts/ends
if strcmpi(myfolders.task,'RS')
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
EEG = remove_electrodewobbles(EEG,lower(cfg.ica.type));

% Detect and remove noisy electrodes
EEG = remove_noisyelec(EEG,cfg.bch);

% Log/report bad channel
EEG = report_badelectrodes(EEG,'individual');

% Detecte extremely bad epochs
EEG = detect_extremelybadepochs2(EEG);

% % Check if EC has eye blinks; Too sensitive?
% if strcmpi(myfolders.task,'RS')
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
EEG.data = EEG.data - (sum(EEG.data,1)/(EEG.nbchan+1));
EEG.ref  = 'average';

% Remove extreme periods
if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3) && ~strcmpi(myfolders.task,'RS')
    EEG = eeg_eegrej(EEG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    EXT = eeg_eegrej(EXT, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);

    if strcmpi(myfolders.task,'MT')
        EMG = eeg_eegrej(EMG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    end
end

% Epoch MT data
if strcmpi(myfolders.task,'RS')
    % Ñot needed as now we can epoch RS data after ICA -> better!
    % EEG = epoch_rsdata(EEG,2*EEG.srate); % Epoch size 2 [s]
    % EXT = epoch_rsdata(EXT,2*EXT.srate); % Epoch size 2 [s]
elseif strcmpi(myfolders.task,'MT')
    % Epoch MT data before ICA, as there is often a lot of noise in between trials
    EEG = pop_epoch(EEG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EXT = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EMG = pop_epoch(EMG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0] ,[]);
end

% ICA
EEG = do_ICA(EEG,cfg);

% ICLabel
EEG = iclabel(EEG);

% Flag artifact ICs
EEG = pop_icflag(EEG,cfg.ica.iclabel);

% ICA log/report
EEG = report_ica(EEG);

% Epoch MMN
if strcmpi(myfolders.task,'MMN')
    EEG = pop_epoch(EEG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0),cfg.trg.mmn{2},'epochinfo','yes');
    EXT  = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0),cfg.trg.mmn{2},'epochinfo','yes');
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0] ,[]);
end
% Epoch SART: 1. wrt visual stimuli and 2. wrt response times
if strcmpi(myfolders.task,'SART')
    EEG0 = EEG; EXT0 = EXT;
    EEG = pop_epoch(EEG0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0),cfg.trg.sart1{2},'epochinfo','yes');
    EXT = pop_epoch(EXT0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0),cfg.trg.sart1{2},'epochinfo','yes');
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0] ,[]);

    % This does not make sure that epochs are only from correct Go
    EEG2 = pop_epoch(EEG0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0),cfg.trg.sart2{2},'epochinfo','yes');
    EXT2 = pop_epoch(EXT0,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0),cfg.trg.sart2{2},'epochinfo','yes');
    EEG2 = pop_rmbase(EEG2,[(EEG2.xmin)*1000 0] ,[]);
    clearvars EEG0 EXT0
elseif strcmpi(myfolders.task,'RS')
    % 2s w/ 0.75 overlap -> [1,5,9,...] -> [1:4:end] trials are unique!
    EEG = epoch_rsdata2(EEG,cfg.trg.rs{1},cfg.trg.rs{2});
    EXT = epoch_rsdata2(EXT,cfg.trg.rs{1},cfg.trg.rs{2});
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

% Report muscle and blink artifact leftovers
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
fprintf('\n%s: Saving preprocessed data...\n\n',subject.id);
if ~strcmpi(myfolders.task,'SART')
    save(fullfile(subject.preproc, subject.clnfile),'EEG','cfg','proctime');
else
    save(fullfile(subject.preproc, subject.clnfile),'EEG','EEG2','cfg','proctime');
end

% Close the report
t1 = datetime("now");
dd = round(minutes(diff([t0 t1])));
fprintf(subject.fid,'Finished: %s\n',t1);
fprintf(subject.fid,'Running time: %d min.\n',dd);
fclose(subject.fid);

fprintf('Finished: %s\n',t1);
fprintf('Running time: %d min.\n',dd);

end