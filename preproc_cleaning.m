function issues_to_check = preproc_cleaning(myfolders,id)
%
% Script for EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, February 2024
% =========================================================================
%

% Preprocessing settings
cfg = preproc_parameters;

% =========================================================================
subject         = [];
subject.id      = id;
subject.rawdata = fullfile(myfolders.rawdata, subject.id);
subject.preproc = fullfile(myfolders.preproc, subject.id);
% subject.icapath = fullfile(subject.preproc, cfg.ica.type);
subject.icafile = [subject.id '_' myfolders.visit '_' myfolders.task '_icadata_' cfg.rnum];
subject.clnfile = [subject.id '_' myfolders.visit '_' myfolders.task '_cleandata_' cfg.rnum '.mat'];

% Find datasets to load
[subject.datablocks, NBLK] = list_datasets(subject.rawdata,myfolders.task);

if exist(subject.preproc,'dir')~=7, mkdir(subject.preproc); end
% if exist(subject.icapath,'dir')~=7, mkdir(subject.icapath); end

% Prepare text file for the report
subject.fid = fopen(fullfile(subject.preproc,['report_preprocess_v' cfg.rnum '.txt']),'w+');

% Timestamp and basic info
t0 = datetime("now");
proctime = strrep(strrep(char(t0),':','-'),' ','-');
fprintf('\n\n%s | %s | %s dataset\n',myfolders.group,subject.id,myfolders.task);
fprintf(subject.fid,'%s | %s | %s dataset\n\n',myfolders.group,subject.id,myfolders.task);
fprintf(subject.fid,'Code version %s\n',cfg.rnum);
fprintf(subject.fid,'Start: %s\n',t0);

% Load data
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

% Separate EMG before EEG cleaning (EMG cleaning???)
if strcmpi(myfolders.task,'MT')
    chanemg = {EEG(1).chanlocs(strcmp({EEG(1).chanlocs.type},'EMG')).labels};
    EMG = pop_select(EEG,'channel',chanemg);
    EEG = pop_select(EEG,'rmchannel',chanemg);
end

% Merge datasets
EEGM = pop_mergeset(EEG,1:NBLK);
if strcmpi(myfolders.task,'MT')
    EMG = pop_mergeset(EMG,1:NBLK);
end

% Mark where each RS block starts/ends
if strcmpi(myfolders.task,'RS')
    EEGM = make_rsmasks(EEGM,cellfun(@(x,y) size(x,2),{EEG(:).data}));
end

% Remove electrode wobbles and pops
EEGM = remove_electrodewobbles(EEGM);

% Detect and remove noisy electrodes
EEGM = remove_noisyelec(EEGM,cfg.bch);

% Bad channel log/report
EEGM = report_badelectrodes(EEGM);

% Detecte extremely bad epochs
% EEGM = detect_extremelybadepochs(EEGM); % Too sensitive?
EEGM = detect_extremelybadepochs2(EEGM);

% Check if EC has eye blinks
if strcmpi(myfolders.task,'RS')
    EEGM = check_eyesclosedeyeblinks(EEGM);
end

% Do Multi-channel Wiener filtering
EEGM = reduce_artifacts(EEGM,cfg.bch);

% Separate EXT channels before ICA
chanext = {EEGM.chanlocs(strcmp({EEGM.chanlocs.type},'EXT')).labels};
EXT  = pop_select(EEGM,'channel',chanext);
EEGM = pop_select(EEGM,'rmchannel',chanext);

% Interpolate bad electrodes
if ~isempty(EEGM.ALSUTRECHT.badchaninfo.badElectrodes)
    EEGM = pop_interp(EEGM,chanlocs);
end

% Common-average before ICA
EEGM.data = EEGM.data - (sum(EEGM.data,1)/(EEGM.nbchan+1));
EEGM.ref  = 'average';

% Remove extreme periods, temporarily maybe?
if ~isempty(EEGM.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)
    EEGM = eeg_eegrej(EEGM, EEGM.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    EXT  = eeg_eegrej(EXT, EEGM.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    if strcmpi(myfolders.task,'MT')
        EMG = eeg_eegrej(EMG, EEGM.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    end
    if strcmpi(myfolders.task,'RS')
        EEGM.ALSUTRECHT.blockinfo.rs_mask(:,EEGM.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = [];
    end
end

% Epoch MT data before ICA, as there is often a lot of noise in between trials
if strcmpi(myfolders.task,'MT')
    EEGM = pop_epoch(EEGM,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EXT  = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EMG  = pop_epoch(EMG,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0),cfg.trg.mt{2},'epochinfo','yes');
    EEGM = pop_rmbase(EEGM,[(EEGM.xmin)*1000 0] ,[]);
end

% ICA
assert(getrank(EEGM.data)>cfg.ica.icMax);
EEGM = pop_runica(EEGM,'icatype',lower(cfg.ica.type),'extended',1,'pca',cfg.ica.icMax,'maxsteps',1000); % ,'lrate',1e-5 //// 0.00065/log(cfgica.icMax)
EEGM = eeg_checkset(EEGM,'ica');

% Epoch MMN/SART
if strcmpi(myfolders.task,'MMN') || strcmpi(myfolders.task,'SART')
    EEGM = pop_epoch(EEGM,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.(lower(myfolders.task)){1},'Uniformoutput',0),cfg.trg.(lower(myfolders.task)){2},'epochinfo','yes');
    EXT  = pop_epoch(EXT,arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.(lower(myfolders.task)){1},'Uniformoutput',0),cfg.trg.(lower(myfolders.task)){2},'epochinfo','yes');
    EEGM = pop_rmbase(EEGM,[(EEGM.xmin)*1000 0] ,[]);
end

% ICLabel
EEGM = iclabel(EEGM);

% Flag artifact ICs
EEGM = pop_icflag(EEGM,cfg.ica.iclabel);

% ICA log/report
EEGM = report_ica(EEGM);

% Save the dataset
% EEGM = pop_saveset(EEGM,'filename',[subject.icafile '.set'],'filepath',subject.preproc);

% Do wICA
EEGM = merge_eeglabsets(EEGM,EXT);
EEGM = do_wICA(EEGM);
% EEGM = do_CEEMDAN_ICA(EEGM);

% Merge EMG
if strcmpi(myfolders.task,'MT')
    EEGM = merge_eeglabsets(EEGM,EMG);
end

% =========================================================================
% Record warnings about potential issues
issues_to_check.aFileName = subject.id;

issues_to_check.FlatElectrodesDiscrepancy = EEGM.ALSUTRECHT.badchaninfo.flatElectrodesDiscrepancy;
if length(EEGM.ALSUTRECHT.badchaninfo.badElectrodes)/128>0.2
    issues_to_check.RejectedTooManyElectrodes = length(EEGM.ALSUTRECHT.badchaninfo.badElectrodes);
else
    issues_to_check.RejectedTooManyElectrodes = 0;
end
if EEGM.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier>0.2
    issues_to_check.HighProportionExcludedAsExtremeOutlier = EEGM.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier;
else
    issues_to_check.HighProportionExcludedAsExtremeOutlier = 0;
end
if EEGM.ALSUTRECHT.MWF.EMG.ProportionOfDataShowingMuscleActivityTotal > 0.5
    issues_to_check.HighProportionOfEMG = EEGM.ALSUTRECHT.MWF.EMG.ProportionOfDataShowingMuscleActivityTotal;
else
    issues_to_check.HighProportionOfEMG = 0;
end
MFWrounds = fields(EEGM.ALSUTRECHT.MWF);
for j = 1:length(MFWrounds)
    issues_to_check.(['MWF' MFWrounds{j} 'Status1']) =  EEGM.ALSUTRECHT.MWF.(MFWrounds{j}).status;

    if isnan(EEGM.ALSUTRECHT.MWF.(MFWrounds{j}).signalToErrorRatio) || isnan(EEGM.ALSUTRECHT.MWF.(MFWrounds{j}).artifactToResidueRatio)
        issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  false;
    else
        issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  true;
    end

    if EEGM.ALSUTRECHT.MWF.(MFWrounds{j}).proportionMarkedForMWF>0.6
        issues_to_check.(['MWF' MFWrounds{j} 'BadData']) =  EEGM.ALSUTRECHT.MWF.(MFWrounds{j}).proportionMarkedForMWF;
    else
        issues_to_check.(['MWF' MFWrounds{j} 'BadData']) =  0;
    end
end
% if ~isempty(tmpLabels)
%     issues_to_check.HighProportionOfBadDataMWF = strjoin(aField(tmpLabels),', ');
% else
%     issues_to_check.HighProportionOfBadDataMWF = 0;
% end

if EEGM.ALSUTRECHT.ica.proportionArtifactICsReducedbywICA>0.3
    issues_to_check.HighProportionOfArtifactICs = EEGM.ALSUTRECHT.ica.proportionArtifactICsReducedbywICA;
else
    issues_to_check.HighProportionOfArtifactICs = 0;
end
issues_to_check.DataTooShortForValidICA = EEGM.ALSUTRECHT.ica.DataLengthForValidICA;
if strcmpi(myfolders.task,'RS')
    if EEGM.ALSUTRECHT.blockinfo.ec_blinks>0
        issues_to_check.ECEyeBinksDetected = EEGM.ALSUTRECHT.blockinfo.ec_blinks;
    else
        issues_to_check.ECEyeBinksDetected = 0;
    end
end
EEGM.ALSUTRECHT.issues_to_check = issues_to_check;

% =========================================================================
% Save cleaned data
fprintf('\n%s: Saving preprocessed data...\n\n',subject.id);
save(fullfile(subject.preproc, subject.clnfile),'EEGM','cfg','proctime');

% Close the report
t1 = datetime("now");
dd = round(minutes(diff([t0 t1])));
fprintf(subject.fid,'\nFinish: %s\n',t1);
fprintf(subject.fid,'Running time: %d min.\n',dd);
fclose(subject.fid);

%%
% =========================================================================
%                              HELPER FUNCTIONS
% =========================================================================

% % =========================================================================
% % 1
% % =========================================================================
% function EEG = fix_noisemask(EEG,noisemask0)
% % Make a mask for each trial
% % eg
% % 111 000 000
% % 000 111 000
% % 000 000 111
% N = cellfun(@(x,y) size(x,2),{EEG(:).data});
% assert(size(EEG(1).data,1)<=136);
%
% NBLK = length(EEG);
% mask_trial = false(NBLK,sum(N));
% for j = 1:NBLK
%     if j == 1
%         mask_trial(j,1:N(1)) = true;
%     else
%         mask_trial(j,sum(N(1:j-1))+1:sum(N(1:j))) = true;
%     end
% end
%
% NCHN = EEG(1).nbchan+1;
% % noisemask = cell(1,NBLK);
% for j = 1:NBLK
%     % noisemask{j} = [noisemask0(mask_trial(j,:))];
%     EEG(j).data(NCHN,:) = [noisemask0(mask_trial(j,:))];
%     EEG(j).nbchan = NCHN;
%     EEG(j).chanlocs(NCHN).labels = 'NOISE';
%     EEG(j).chanlocs(NCHN).type   = 'MSK';
% end
% % Check
% EEG = eeg_checkset(EEG);
