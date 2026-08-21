
params = checkBlinkerDefaults(struct(), getBlinkerDefaults(EEG));

folder_blink = fullfile(subject.preproc, 'blinks');
if exist(folder_blink,'dir') ~= 7, mkdir(folder_blink); end

params.subjectID  = '';
params.experiment = '';
params.uniqueName = '';
params.task       = '';
params.fileName   = '';
% params.startDate  = '';
% params.startTime  = '';
params.blinkerSaveFile       = fullfile(folder_blink, 'AllUnrefBlinkSummary.mat');
params.blinkerDumpDir        = folder_blink;
params.dumpBlinkerStructures = true;
params.dumpBlinkImages       = true;
params.dumpBlinkPositions    = true;
% params.keepSignals           = false;
% params.showMaxDistribution   = false;
params.verbose               = true;
params.signalLabels          = {'VEOGS', 'VEOGI'};
params.excludeLabels         = {'HEOGL', 'HEOGR', 'LM', 'RM', 'ECGL', 'ECGR'};

[EEG2, com, blinks, blinkFits, blinkProperties, blinkStatistics, params] = pop_blinker(EEG, params);



x = EEG.icawinv(:, 1);
mytopoplot(x, [], '');





% Step 6. Detect outlier channels (05/05/2026).
originalEEG  = EEG;
chIdx_withoutCz = 1:EEG.nbchan;

EEG_B4ChRej = EEG;  % Keep original channel indexing.

numEpochs = 1000;
permillWindowLength = floor(EEG.pnts/numEpochs); % epoch length must be >=3.
if permillWindowLength < 3
    error('Data too short: consider using less than 1000 epochs at line 32.');
end
trimmedData = EEG.data(chIdx_withoutCz, 1:permillWindowLength*1000);
trimmedData3D = reshape(trimmedData, [size(trimmedData,1) permillWindowLength numEpochs]);

chStdMatrix = squeeze(std(trimmedData3D, 0, 2)); % ch x epoch.
chStdMedian = median(chStdMatrix, 2); % ch.
flatChannelIdx = find(chStdMedian < 0.1);
nonflatChIdx   = setdiff(chIdx_withoutCz, flatChannelIdx);
chStdMedianLog = log(chStdMedian);
chRejLogThreshold = median(chStdMedianLog(nonflatChIdx)) + mad(chStdMedianLog(nonflatChIdx),1)*1.4826*5; % mad(X,1) for using median.
outlierChannelIdx = find(chStdMedianLog > chRejLogThreshold);
badChannelIdx = sort([flatChannelIdx;outlierChannelIdx]);

EEG = pop_select(EEG, 'nochannel', badChannelIdx);
EEG.etc.badChannelIdx = badChannelIdx;









eeglab_path = 'C:\DATA\MATLAB\EEG\2_OTHER_DATA\mri\template\eeg\BEM';
templateChannelFilePath = fullfile(eeglab_path, 'electrodes.mat');
hdmFilePath             = fullfile(eeglab_path, 'headmodel.mat');

EEG = pop_dipfit_settings(EEG, 'hdmfile', hdmFilePath, 'coordformat', 'MNI',...
    'mrifile', '[EEGLABroot]/eeglab/plugins/dipfit2.3/standard_BEM/standard_mri.mat',...
    'chanfile', templateChannelFilePath, 'coord_transform', coordinateTransformParameters,...
    'chansel', 1:EEG.nbchan);

EEG = pop_multifit(EEG, 1:EEG.nbchan,'threshold', 100, 'dipplot','off','plotopt',{'normlen' 'on'});