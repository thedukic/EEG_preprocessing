function issues_to_check = preproc_cleaning(myfolders,id)
%
% Script for EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, March 2024
% =========================================================================
%

% Load preprocessing settings
cfg = preproc_parameters;

% Define paths and files
subject         = [];
subject.id      = id;
subject.rawdata = fullfile(myfolders.rawdata, subject.id);
subject.preproc = fullfile(myfolders.preproc, subject.id);
% subject.icapath = fullfile(subject.preproc, cfg.ica.type);
% subject.icafile = [subject.id '_' myfolders.visit '_' myfolders.task '_icadata_' cfg.rnum];
subject.clnfile = [subject.id '_' myfolders.visit '_' myfolders.task '_cleandata_' cfg.rnum '.mat'];

% Find datasets to load
[subject.datablocks, NBLK] = list_datasets(subject.rawdata,myfolders.task);

% Make folders
if exist(subject.preproc,'dir')~=7, mkdir(subject.preproc); end
% if exist(subject.icapath,'dir')~=7, mkdir(subject.icapath); end

% Open a report
subject.fid = fopen(fullfile(subject.preproc,['report_preprocess_v' cfg.rnum '.txt']),'w+');
% Log basic info
t0 = datetime("now");
proctime = strrep(strrep(char(t0),':','-'),' ','-');
fprintf('\n\n%s | %s | %s dataset\n',myfolders.group,subject.id,myfolders.task);
fprintf(subject.fid,'%s | %s | %s dataset\n\n',myfolders.group,subject.id,myfolders.task);
fprintf(subject.fid,'Code version %s\n',cfg.rnum);
fprintf(subject.fid,'Start: %s\n',t0);

% Load data
if ~isempty(subject.datablocks)
    if strcmpi(myfolders.task,'MT')
        EEG = pop_biosig(subject.datablocks,'channels',1:168);
        EEG = pop_select(EEG,'rmchannel',[131 132 143:160]);
    else
        EEG = pop_biosig(subject.datablocks,'channels',1:136);

        % Rare cases of RS datasets
        if any([EEG(:).srate] == 2048)
            EEG = pop_biosig(subject.datablocks,'channels',1:168);
            EEG = pop_select(EEG,'channel',[1:128 161:168]);
        end
    end
else
    warning([subject.id ' is missing ' myfolders.task ' data!']);
    EEG.ALSUTRECHT.subject.id = subject.id;
    EEG = report_issues(EEG,myfolders.task);
    issues_to_check = EEG.ALSUTRECHT.issues_to_check;

    return;
end

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
EEG = extract_eventinfo(EEG,myfolders.task);

% Detect flat channels on raw data
EEG = remove_flatelec(EEG,cfg.bch);

% Reference
% Doing this first and then filtering might help line-noise removal
% EEG = pop_reref(EEG,'LM');
chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
for j = 1:NBLK
    EEG(j).data(chaneeg,:) = EEG(j).data(chaneeg,:) - trimmean(EEG(j).data(chaneeg,:),10);
    EEG(j).ref = 'average';
end

% Maybe better to reorder:
% 1. Filter
% 2. Cute ends
% 3. Remove line noise

% Remove the line noise
EEG = reduce_linenoise(EEG);

% Filter
if strcmpi(myfolders.task,'RS') || strcmpi(myfolders.task,'MT')
    % RS/MT
    EEG = filter_signal(EEG,cfg.flt.rsmt.lp,cfg.flt.rsmt.hp,'eeglab');
else
    % MMN/SART (ERP)
    EEG = filter_signal(EEG,cfg.flt.erp.lp,cfg.flt.erp.hp,'eeglab');
end

% Bipolar EXT
EEG = make_extbipolar(EEG,myfolders.task);

% Cut block ends
EEG = remove_datasetends(EEG,myfolders.task);

% Mark where each RS block starts/ends
if strcmpi(myfolders.task,'RS')
    EEG = make_rsmasks(EEG);
end

% Merge datasets
EEG = pop_mergeset(EEG,1:NBLK);

% Separate EMG before EEG cleaning (EMG cleaning???)
if strcmpi(myfolders.task,'MT')
    chanemg = {EEG(1).chanlocs(strcmp({EEG(1).chanlocs.type},'EMG')).labels};
    EMG = pop_select(EEG,'channel',chanemg);
    EEG = pop_select(EEG,'rmchannel',chanemg);
end

% % Make a copy
% EEGRAW = EEG;

% Remove electrode wobbles and pops
EEG = remove_electrodewobbles(EEG,lower(cfg.ica.type));

% Detect and remove noisy electrodes
EEG = remove_noisyelec(EEG,cfg.bch);

% Log/report bad channel 
EEG = report_badelectrodes(EEG,'individual');

% Detecte extremely bad epochs
EEG = detect_extremelybadepochs2(EEG);

% Check if EC has eye blinks
% Estiamte blinks in EO and compair to detections in EC
if strcmpi(myfolders.task,'RS')
    EEG = check_eyesclosedeyeblinks(EEG); % Too sensitive?
end

% Do Multi-channel Wiener filtering
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
if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)
    EEG = eeg_eegrej(EEG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    EXT = eeg_eegrej(EXT, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);

    if strcmpi(myfolders.task,'MT')
        EMG = eeg_eegrej(EMG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    end
    if strcmpi(myfolders.task,'RS')
        EEG.ALSUTRECHT.blockinfo.rs_mask(:,EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = [];
        assert(size(EEG.data,2)==size(EEG.ALSUTRECHT.blockinfo.rs_mask,2));
    end
end

% Epoch MT data before ICA, as there is often a lot of noise in between trials
if strcmpi(myfolders.task,'MT')
    EEG = pop_epoch(EEG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EXT = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EMG = pop_epoch(EMG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0] ,[]);
end

% ICA
assert(getrank(EEG.data)>cfg.ica.icMax);
EEG = pop_runica(EEG,'icatype',lower(cfg.ica.type),'extended',1,'pca',cfg.ica.icMax,'maxsteps',1000); % ,'lrate',1e-5 //// 0.00065/log(cfgica.icMax)
EEG = eeg_checkset(EEG,'ica');

% ICLabel
EEG = iclabel(EEG);

% Flag artifact ICs
EEG = pop_icflag(EEG,cfg.ica.iclabel);

% ICA log/report
EEG = report_ica(EEG);

% % Epoch MMN/SART, but it would be good to epoch SART around RT too
% if strcmpi(myfolders.task,'MMN') || strcmpi(myfolders.task,'SART')
%     EEGM = pop_epoch(EEGM,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.(lower(myfolders.task)){1},'Uniformoutput',0),cfg.trg.(lower(myfolders.task)){2},'epochinfo','yes');
%     EXT  = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.(lower(myfolders.task)){1},'Uniformoutput',0),cfg.trg.(lower(myfolders.task)){2},'epochinfo','yes');
%     EEGM = pop_rmbase(EEGM,[(EEGM.xmin)*1000 0] ,[]);
% end

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
% compare_visually(EEG,EEGRAW,myfolders.task,cfg.trg);

% Report muscle and blink artifact leftovers
EEG = report_leftovers(EEG);

% Merge EMG
if strcmpi(myfolders.task,'MT')
    EEG = merge_eeglabsets(EEG,EMG);
end

% Record warnings about potential issues
EEG = report_issues(EEG,myfolders.task);
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
% fprintf('\n\n%s | %s | %s dataset\n',myfolders.group,subject.id,myfolders.task);
fprintf(subject.fid,'Finish: %s\n',t1);
fprintf(subject.fid,'Running time: %d min.\n',dd);
fclose(subject.fid);

end