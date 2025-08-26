function preproc_cleaning1(myPaths,id)
% =========================================================================
%
% Script for EEG data preprocessing
% ALS Centre, University Medical Centre Utrecht
%
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

% Find files
subject.datablocks = list_datasets(subject.rawdata,myPaths.task);

% Report
fprintf('\n==================================================================\n');
fprintf('%s | %s | %s dataset | processing part 1\n',myPaths.group,subject.id,myPaths.task);
fprintf('==================================================================\n');

% Load data
if ~isempty(subject.datablocks)
    % Report if only EC or EO are found for RS proc
    EEG = load_biosemidata(subject,myPaths);
else
    warning([subject.id ' is missing ' myPaths.task ' data. Skipping...']); return;
end

% Make a folder
if exist(subject.preproc,'dir') ~= 7
    mkdir(subject.preproc);
else
    warning('The folder already exists. The new and the old files might be mixed:'); warning('%s',subject.preproc);
end

% Open a report
t0 = datetime("now");
procTimeTags = {myPaths.proctime; strrep(strrep(char(t0),':','-'),' ','-')};
subject.fid  = fopen(fullfile(subject.preproc,['report_preprocess_v' myPaths.rnum '.txt']),'w+');

fprintf(subject.fid,'%s | %s | %s dataset\n\n',myPaths.group,subject.id,myPaths.task);
fprintf(subject.fid,'Code version %s\n',myPaths.rnum);
fprintf(subject.fid,'Started: %s\n',t0);

% Manually fix some datasets
EEG = do_manualfix(EEG,subject,myPaths);

% Fix events
EEG = fix_events(EEG);

% Add subject/channel info
EEG = add_info(EEG,subject,cfg);

% Remove CMS drop-outs
EEG = detect_dropouts(EEG);

% Resample
EEG = do_resampling(EEG,256);

% Detect flat channels
EEG = remove_flatelectrodes(EEG,cfg.bch);

% Remove extremely bad epochs
[EEG, flagExclude] = remove_extremeperiods2(EEG);

% Check if it is worth continuing 
if flagExclude, warning([subject.id ' has very noisy data. Skipping...']); return; end

% Keep event info
EEG = extract_eventinfo(EEG,cfg.trg);

% Filter highpass only
EEG = do_filtering(EEG,'highpass',cfg.flt);

% Reference
EEG = do_reref(EEG,'aRobust');

% Remove line noise
EEG = reduce_linenoise(EEG);

% Remove line noise leftovers
EEG = reduce_linenoiseleftovers(EEG);

% Remove other spectral peaks
EEG = reduce_spectrapeaks(EEG);

% % Cut block ends - not needed
% EEG = remove_datasetends(EEG);

% Make EXT bipolar
EEG = make_extbipolar(EEG);

% Mark where each RS block starts/ends
if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO') || strcmpi(myPaths.task,'EC')
    EEG = make_rsmasks(EEG);
end

% Merge blocks
EEG = merge_eeglabblocks(EEG);

% Reduce sparse artifacts
EEG = reduce_sparseartifacts(EEG);

% % TODO: Check if EC-RS has eye blinks - Currently too sensitive
% if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EC')
%     EEG = check_eyesclosedeyeblinks(EEG);
% end

% % TODO: Clean EMG
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

% Remove noisy electrodes
EEG = remove_noisyelectrodes(EEG,cfg.bch);

% Log/report bad channel
EEG = report_badelectrodes(EEG);

% % Remove extremely bad epochs - already done
% if ~strcmpi(myPaths.task,'MT')
%     EEG = remove_extremeperiods(EEG);
% else
%     [EEG, EMG] = remove_extremeperiods(EEG,EMG);
% end

% Reduce artifacts
if cfg.mwf.do
    EEG = reduce_artifacts(EEG,cfg.bch);
else
    fprintf('MWF cleaning: Skipped! Turned off by the user.');
end

% Separate EXT channels before ICA
chanext = {EEG.chanlocs(strcmp({EEG.chanlocs.type},'EXT')).labels};
EXT = pop_select(EEG,'channel',chanext);
EEG = pop_select(EEG,'rmchannel',chanext);

% Interpolate bad electrodes
EEG = do_channelinterp(EEG,'spherical');

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

% Merge channels/sets
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