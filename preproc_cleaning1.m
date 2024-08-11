function checkReport = preproc_cleaning1(myPaths,id)
%
% Script for EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
% =========================================================================
% SDukic edits
% v1, August 2024
% =========================================================================
% TODO
% 1. Targeted ECG cleaning
% 2.
%

% Load preprocessing settings
cfg = preproc_parameters;

% Define paths and files
subject         = [];
subject.id      = id;
subject.task    = myPaths.task;
subject.group   = myPaths.group;
subject.visit   = myPaths.visit;
subject.rawdata = fullfile(myPaths.rawdata, subject.id);
subject.preproc = fullfile(myPaths.preproc, subject.id);
subject.icadata = fullfile(subject.preproc, upper(cfg.ica.type2));
subject.clnfile = [subject.id '_' myPaths.visit '_' myPaths.task '_cleandata_' cfg.rnum '.mat'];

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
            maskfs2k = [EEG(:).srate] == 2048;
            EEG(maskfs2k) = pop_biosig(subject.datablocks(maskfs2k),'channels',1:168);
            EEG(maskfs2k) = pop_select(EEG(maskfs2k),'channel',[1:128 161:168]);
        end
    end
else
    warning([subject.id ' is missing ' myPaths.task ' data!']);
    EEG.ALSUTRECHT.subject.id   = subject.id;
    EEG.ALSUTRECHT.subject.task = subject.task;
    EEG = report_issues(EEG);
    checkReport = EEG.ALSUTRECHT.issues_to_check;
    return;
end

% If loaded, then make a folder
if exist(subject.preproc,'dir')~=7, mkdir(subject.preproc); end

% Open a report
t0 = datetime("now");
procTimeTags = {myPaths.proctime; strrep(strrep(char(t0),':','-'),' ','-')};
subject.fid  = fopen(fullfile(subject.preproc,['report_preprocess_v' cfg.rnum '.txt']),'w+');

fprintf(subject.fid,'%s | %s | %s dataset\n\n',myPaths.group,subject.id,myPaths.task);
fprintf(subject.fid,'Code version %s\n',cfg.rnum);
fprintf(subject.fid,'Started: %s\n',t0);

% Manually fix datasets in some rare cases
if strcmp(subject.id,'C50') && strcmpi(myPaths.task,'SART')
    % C48 is very strage - low quality data?
    warning([subject.id 'has swapped C- and B- set. Fixing that now...']);
    for j = 1:NBLK
        EEG(j).data(33:96,:) = [EEG(j).data(65:96,:); EEG(j).data(33:64,:)];
    end
elseif strcmp(subject.id,'ALS26603') && strcmpi(myPaths.task,'MMN')
    warning([subject.id ' has MMN2 file currpted. Removing it now...']);
    EEG(2) = []; NBLK = length(EEG);

elseif strcmp(subject.id,'ALS26603') && strcmpi(myPaths.task,'SART')
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
    EEG(j).ALSUTRECHT.cfg = cfg;
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

% Filter highpass only
EEG = do_filtering(EEG,1,cfg.flt);

% Remove line noise
EEG = reduce_linenoise(EEG);

% Check leftovers
EEG = detect_linenoiseleftovers(EEG);

% Cut block ends
EEG = remove_datasetends(EEG);

% Bipolar EXT
EEG = make_extbipolar(EEG);

% Mark where each RS block starts/ends
if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO') || strcmpi(myPaths.task,'EC')
    EEG = make_rsmasks(EEG);
end

% Merge datasets
EEG = pop_mergeset(EEG,1:NBLK);

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
EEG = remove_electrodewobbles(EEG,lower(cfg.ica.type1));

% Detect and remove noisy electrodes
EEG = remove_noisyelec(EEG,cfg.bch);

% Log/report bad channel
EEG = report_badelectrodes(EEG);

% Detect extremely bad epochs
EEG = detect_extremelybadepochs2(EEG);

% Check if EC has eye blinks; Too sensitive?
% if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EC')
%     EEG = check_eyesclosedeyeblinks(EEG);
% end

% Make a copy
EEGRAW = EEG;

% Reduce artifacts
EEG = reduce_artifacts(EEG,cfg.bch);

% Report MWF metrics
EEG = report_mwf(EEG,EEGRAW);
clearvars EEGRAW

% Remove extremely bad epochs
if ~isempty(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3)
    EEG = eeg_eegrej(EEG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);

    if strcmpi(myPaths.task,'MT')
        EMG = eeg_eegrej(EMG, EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs3);
    end
end

% Separate EXT channels before ICA
chanext = {EEG.chanlocs(strcmp({EEG.chanlocs.type},'EXT')).labels};
EXT = pop_select(EEG,'channel',chanext);
EEG = pop_select(EEG,'rmchannel',chanext);

% Interpolate bad electrodes
if ~isempty(EEG.ALSUTRECHT.badchaninfo.badElectrodes)
    EEG = pop_interp(EEG,chanlocs,'spherical');
end

% Filter lowpass only
EEG = do_filtering(EEG,2,cfg.flt);

% Common-average before ICA
EEG = do_reref(EEG,'aRegular');

% ICA
EEG = do_ICA(EEG,cfg);

% Flag artifact ICs
EEG = detect_badICs(EEG,EXT,cfg);

% Do wICA
EEG = do_wICA(EEG,EXT,cfg);

% Interim data saving
% EEG = pop_saveset(EEG,'filename',[subject.icafile '.set'],'filepath',subject.preproc);

% Visually check the cleaning
% compare_visually(EEG,EEGRAW,cfg.trg);

% Report ICA
EEG = report_ICA(EEG);

% Report artifact leftovers
EEG = report_leftovers(EEG,EXT,cfg);

% Epoch
if strcmpi(myPaths.task,'MMN')
    condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG,condLabel,cfg.trg.mmn{2},'epochinfo','yes');
    EXT = pop_epoch(EXT,condLabel,cfg.trg.mmn{2},'epochinfo','yes');

elseif strcmpi(myPaths.task,'SART')
    EEG0 = EEG; EXT0 = EXT;

    % SART wrt visual stimuli
    condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0);
    EEG  = pop_epoch(EEG0,condLabel,cfg.trg.sart1{2},'epochinfo','yes');
    EXT  = pop_epoch(EXT0,condLabel,cfg.trg.sart1{2},'epochinfo','yes');

    % SART wrt response times
    condLabel2 = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0);
    EEG2 = pop_epoch(EEG0,condLabel2,cfg.trg.sart2{2},'epochinfo','yes');
    EXT2 = pop_epoch(EXT0,condLabel2,cfg.trg.sart2{2},'epochinfo','yes');

    EEG.ALSUTRECHT.SART.type  = 'StimulusLocked';
    EXT.ALSUTRECHT.SART.type  = 'StimulusLocked';
    EEG2.ALSUTRECHT.SART.type = 'ResponseLocked';
    EXT2.ALSUTRECHT.SART.type = 'ResponseLocked';

    clearvars EEG0 EXT0

elseif strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO') || strcmpi(myPaths.task,'EC')
    % 2s w/ 0.75 overlap
    % EEG = epoch_rsdata3(EEG,cfg.trg.rs{1},cfg.trg.rs{2}); % OK if proc EO/EC only
    % EXT = epoch_rsdata3(EXT,cfg.trg.rs{1},cfg.trg.rs{2});
    EEG = epoch_rsdata2(EEG,cfg.trg.rs{1},cfg.trg.rs{2});   % OK if proc EO+EC together
    EXT = epoch_rsdata2(EXT,cfg.trg.rs{1},cfg.trg.rs{2});

elseif strcmpi(myPaths.task,'MT')
    % Maybe before ICA, as there is often a lot of noise before/after the contrations
    condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0);
    EEG = pop_epoch(EEG,condLabel,cfg.trg.mt{2},'epochinfo','yes');
    EXT = pop_epoch(EXT,condLabel,cfg.trg.mt{2},'epochinfo','yes');
    EMG = pop_epoch(EMG,condLabel,cfg.trg.mt{2},'epochinfo','yes');
end

% Merge
EEG = merge_eeglabsets(EEG,EXT);
if strcmpi(myPaths.task,'SART'), EEG2 = merge_eeglabsets(EEG2,EXT2); end
if strcmpi(myPaths.task,'MT'),   EEG = merge_eeglabsets(EEG,EMG); end
clearvars EXT EXT2 EMG

% Record warnings about potential issues
EEG = report_issues(EEG);
if strcmpi(myPaths.task,'SART'), EEG2.ALSUTRECHT.issues_to_check = EEG.ALSUTRECHT.issues_to_check; end

% Save cleaned data
fprintf('\n%s: Saving the preprocessed data (part 1)...\n',subject.id);
preprocReport = EEG.ALSUTRECHT;
if ~strcmpi(myPaths.task,'SART')
    EEG.icaact = [];
    save(fullfile(subject.preproc,subject.clnfile),'EEG','preprocReport','cfg','procTimeTags');
else
    EEG.icaact  = [];
    EEG2.icaact = [];
    save(fullfile(subject.preproc,subject.clnfile),'EEG','EEG2','preprocReport','cfg','procTimeTags');
end

% Return
checkReport = EEG.ALSUTRECHT.issues_to_check;
checkReport = rmfield(checkReport,'Medianvoltageshiftwithinepoch');

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