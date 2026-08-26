function EEG = detect_badcomponents(EEG,EXT,cfg)

fprintf('\nDetecting bad ICs...\n\n');

% Double-check
EEG = eeg_checkset(EEG,'ica');

% Make sure IC activations are present
% EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
% EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
EEG.icaact = [];
ICAdata = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
ICAdata = ICAdata(:,:);

NICA = size(ICAdata,1);

%% ========================================================================
% 1. ICLabel
EEG = iclabel(EEG);

% Flag artifact ICs
EEG = pop_icflag(EEG,cfg.ica.iclabel);

% Log IClabel info
EEG.ALSUTRECHT.ica.ICLabel.bics = find(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.ICLabel.clss = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel.pvec, EEG.ALSUTRECHT.ica.ICLabel.cvec] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);

%% ========================================================================
% 2. Correlation with external channels (ECG / VEOG / HEOG)
chanecg  = find(strcmp({EXT.chanlocs.labels},'ECG')); % if it was recorded!
chanveog = find(strcmp({EXT.chanlocs.labels},'VEOG'));
chanheog = find(strcmp({EXT.chanlocs.labels},'HEOG'));
% maskExt  = [chanecg, chanveog, chanheog];

% If the signal is not recorded, try to estimate it using
% an approximate reconstruciton (ICA-like) weights
if isempty(chanecg)
    warning('ECG signals was not recorded, but we can try to estimate it from the EEG data...'); % disp(extLabels);
    load(fullfile(EEG.ALSUTRECHT.subject.mycodes,'files','ECGsignalweights'),'ECGweights');

    % chaneeg = strcmp({EEG.chanlocs.type},'EEG');
    % ECGdata = ECGweights*EEG.data(chaneeg,:);
    EXTdataECG = ECGweights'*EEG.data(:,:);

    EXTTMP = EXT;
    NCHNECG = EXT.nbchan+1;
    EXTTMP.data(NCHNECG,:)          = EXTdataECG;
    EXTTMP.chanlocs(NCHNECG).labels = 'ECG';
    EXTTMP.chanlocs(NCHNECG).type   = 'EXT';

    ECGsignalLabel = 'Estimated';
else
    EXTdataECG = EXT.data(chanecg,:);
    EXTTMP = EXT;
    ECGsignalLabel = 'Recorded';
end

EXTdataEOG = EXT.data([chanveog chanheog],:);
extLabels  = {'ECG','VEOG','HEOG'};
% extLabels = {EXT.chanlocs(maskExt).labels};
% disp(extLabels);

% Temporarily filter for better detection
% It is fine that it will be double-filtered
% EOG: 1-10 Hz
% ECG: 8-16 Hz / 10-20 Hz
[blEOG, alEOG] = butter(4,10/(EEG.srate/2),'low');
[blECG, alECG] = butter(4,35/(EEG.srate/2),'low');
[bhEOG, ahEOG] = butter(4,1/(EEG.srate/2),'high');
[bhECG, ahECG] = butter(4,1/(EEG.srate/2),'high');

assert(isstable(blEOG,alEOG));
assert(isstable(blECG,alECG));
assert(isstable(bhEOG,ahEOG));
assert(isstable(bhECG,ahECG));

% Filter ICA EOG
ICAdataEOG = do_filteringcore(blEOG,alEOG,ICAdata,EEG.event,EEG.srate);
ICAdataEOG = do_filteringcore(bhEOG,ahEOG,ICAdataEOG,EEG.event,EEG.srate);

% Filter ICA ECG
ICAdataECG = do_filteringcore(blECG,alECG,ICAdata,EEG.event,EEG.srate);
ICAdataECG = do_filteringcore(bhECG,ahECG,ICAdataECG,EEG.event,EEG.srate);

% Transpose
ICAdataEOG = ICAdataEOG';
ICAdataECG = ICAdataECG';

% Filter EXT EOG
EXTdataEOG = do_filteringcore(blEOG,alEOG,EXTdataEOG,EEG.event,EEG.srate);
EXTdataEOG = do_filteringcore(bhEOG,ahEOG,EXTdataEOG,EEG.event,EEG.srate);

% Filter EXT ECG
EXTdataECG = do_filteringcore(blECG,alECG,EXTdataECG,EEG.event,EEG.srate);
EXTdataECG = do_filteringcore(bhECG,ahECG,EXTdataECG,EEG.event,EEG.srate);

% Combine: ECG / VEOG / HEOG
EXTdata = [EXTdataECG; EXTdataEOG]';

% Correlation
% 1. ECG
[ECGmask, ECGbadEpoch] = detect_ecg(EXTTMP,200,ECGsignalLabel);

Nhbeats = size(ECGbadEpoch,1);
corrECG = NaN(NICA,Nhbeats);

for i = 1:Nhbeats
    thisChunk = ECGbadEpoch(i,1):ECGbadEpoch(i,2);
    corrECG(:,i) = corr(ICAdataECG(thisChunk,:).^2,abs(EXTdata(thisChunk,1)));
    % corrECG(:,i) = corr(ICAdataECG(thisChunk,:),EXTdata(thisChunk,1));
end
% figure; imagesc(corrECG);
corrECG = abs(zscore(mean(corrECG,2)));
% figure; plot(corrECG);

% % Cross-trial phase statistics
% % https://mne.tools/stable/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA.find_bads_ecg
% % https://sci-hub.st/https://ieeexplore.ieee.org/document/4536072
% [blECG, alECG] = butter(4,16/(EEG.srate/2),'low');
% [bhECG, ahECG] = butter(4,8/(EEG.srate/2),'high');
% ICAdataECG = do_filteringcore(blECG,alECG,ICAdata,EEG.event,EEG.srate);
% ICAdataECG = do_filteringcore(bhECG,ahECG,ICAdataECG,EEG.event,EEG.srate);
%
% V = my_ctps(ICAdataECG, ECGbadEpoch);
% % figure; plot(zscore(V));

% 2. VEOG
[noiseMask, eyeBlinksEpochs] = detect_veog(EXT,200,false);

Nblinks = size(eyeBlinksEpochs,1);
corrVEOG = NaN(NICA,Nblinks);

for i = 1:Nblinks
    thisChunk = eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2);
    corrVEOG(:,i) = corr(ICAdataEOG(thisChunk,:),EXTdata(thisChunk,2));
end
% figure; imagesc(corrVEOG);
corrVEOG = abs(zscore(mean(corrVEOG,2)));

% 3. HEOG
% [saccadesMask, saccadesEpochs] = detect_heog(EXT,200);
corrHEOG = abs(zscore(corr(ICAdataEOG,EXTdata(:,3))));

% Combine: ECG / VEOG / HEOG
corrMat = [corrECG, corrVEOG, corrHEOG];

% Tresholds: ECG / VEOG / HEOG (R or Zscore)
% corrTreshold = [0.6 0.6 0.6];
corrTreshold = [3 3 3];
[badIC, badICtype] = find(corrMat>corrTreshold);

% figure; imagesc(corrMat); clim([0 1]); colorbar;
% xticks([1 2 3]); xticklabels(extLabels); colormap(brewermap(128,'PuRd'));
% EEG = pop_saveset(EEG,'filename',['TMP.set'],'filepath',EEG.ALSUTRECHT.subject.preproc);

% Prevent false positive HEOG
% Here, HEOG ICs are easily confused by broad L-R dipolar (brain) ICs
% Remove HEOG detections if they are likely brain ICs
maskHEOG = find(badICtype==3);
falseHEOGIC = EEG.ALSUTRECHT.ica.ICLabel.cvec(badIC(maskHEOG)) == 1 & EEG.ALSUTRECHT.ica.ICLabel.pvec(badIC(maskHEOG)) > 0.6;
falseHEOGIC = maskHEOG(falseHEOGIC);

badIC(falseHEOGIC) = [];
badICtype(falseHEOGIC) = [];

% Report
NEXT = length(extLabels);
for i = 1:NEXT
    if any(badICtype==i)
        fprintf('%s ICs (N = %d) were idetified using EXT channel correlation (R>%1.1f).\n',extLabels{i},sum(badICtype==i),corrTreshold(i));
    else
        fprintf('No %s ICs were idetified using EXT channel correlation (R>%1.1f).\n',extLabels{i},corrTreshold(i));
    end
end

% Log
EEG.ALSUTRECHT.ica.extra1.corr = single(corrMat);
EEG.ALSUTRECHT.ica.extra1.bics = badIC;
EEG.ALSUTRECHT.ica.extra1.cvec = badICtype;
EEG.ALSUTRECHT.ica.extra1.clss = extLabels;

% Combine ICLabel and ECG correlation
% ICsMostLikelyHeart = EEG.ALSUTRECHT.ica.ICLabel.cvec == 4 & EEG.ALSUTRECHT.ica.ICLabel.pvec > 0.5;
ICsMostLikelyHeart = false(NICA,1);
ICsMostLikelyHeart(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 4) = true;
ICsMostLikelyHeart(EEG.ALSUTRECHT.ica.extra1.bics(EEG.ALSUTRECHT.ica.extra1.cvec==1)) = true;

%% ========================================================================
% 3A. Detect EMG ICs using freq slopes
% options.muscleFreqIn    = [7, 70];
options.muscleFreqIn = cfg.bch.muscleSlopeFreq;
options.Freq_to_compute = [1, 100];

% Calculate pwelch to enable detection of log-freq log-power slopes,
% indicative of muscle activity
[pxx,fp] = pwelch(ICAdata',size(ICAdata,2),[],size(ICAdata,2),EEG.srate);
FFTout = pxx';
fp = fp';

% Calculate FFT bins
freq = options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
fftBins = zeros(size(FFTout,1),size(freq,2)); % preallocate
for a = 1:size(freq,2)
    [~, index1]=min(abs(fp-((freq(1,a)-0.25))));
    [~, index2]=min(abs(fp-((freq(1,a)+0.25))));
    fftBins(:,a) = mean(FFTout(:,index1:index2),2); % creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
end

% figure;
% nexttile; plot(fp,pxx(:,1));
% nexttile; plot(freq,fftBins(1,:));

% better muscle comp_number identification:
options.muscleFreqEx = [50-2 50+2];

muscleRatio = NaN(1,NICA);
for compNum = 1:NICA
    % Define frequencies to include in the analysis
    if ~isempty(options.muscleFreqIn)
        [~,fin1] = min(abs(options.muscleFreqIn(1) - freq));
        [~,fin2] = min(abs(options.muscleFreqIn(2) - freq));
        freqHz = freq(1,fin1:fin2);
        freqPow = fftBins(compNum,fin1:fin2);
    else
        freqHz = freq;
        freqPow = fftBins(compNum,:);
    end
    % Define frequencies to exclude from fit
    if ~isempty(options.muscleFreqEx)
        [~,fex1] = min(abs(options.muscleFreqEx(1) - freqHz));
        [~,fex2] = min(abs(options.muscleFreqEx(2) - freqHz));
        freqHz(fex1:fex2) = [];
        freqPow(fex1:fex2) = [];
    end
    % Fit linear regression to log-log data
    p = polyfit(log(freqHz),log(freqPow),1);
    % Store the slope
    muscleRatio(compNum) = p(1);
end

% Detect EMG components
% muscleSlopeThreshold = cfg.bch.muscleSlopeThreshold;
muscleSlopeThreshold = -0.5;
ICsMostLikelyMuscle  = muscleRatio >= muscleSlopeThreshold;

% Add IClabel EMG ICs as well
ICsMostLikelyMuscleICLabel = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics) == 2);
ICsMostLikelyMuscle(ICsMostLikelyMuscleICLabel) = true;

if any(ICsMostLikelyMuscle)
    fprintf('\nMuscle ICs (N = %d) were idetified using ICLabel and >=%1.1f slope.\n',sum(ICsMostLikelyMuscle),muscleSlopeThreshold);
else
    fprintf('\nNo muscle ICs were idetified using ICLabel and >=%1.1f slope.\n',muscleSlopeThreshold);
end

% 3B. Use icablinkmetrics for eyeblinks
% EOGdata = mean(EEG.data(ismember({EEG.chanlocs.labels},cfg.ica.blinkchans),:),1);
EOGdata = EXTdata(:,2)';

% Put bandpassed ICA data
EEG.icaact = ICAdataEOG';

try
    icablinkmetricsout = icablinkmetrics(EEG,'ArtifactChannel',EOGdata,'Alpha',0.001,'VisualizeData','False');
    if any(icablinkmetricsout.identifiedcomponents>0)
        fprintf('icablinkmetrics has identified %d eye IC(s).\n',length(icablinkmetricsout.identifiedcomponents));
        % ICsMostLikelyNotBrain(icablinkmetricsout.identifiedcomponents) = true;
        % ICsMostLikelyEye(icablinkmetricsout.identifiedcomponents)      = true;
    else
        fprintf('icablinkmetrics has not identified any eye ICs.\n');
        icablinkmetricsout.identifiedcomponents = [];
        % icablinkmetricsout.metrics.corr_Pvalue  = [];
        % icablinkmetricsout.metrics.conv_Pvalue  = [];
        % icablinkmetricsout.metrics.perc_Pvalue  = [];
    end
catch
    warning('icablinkmetrics has failed...');
    icablinkmetricsout.identifiedcomponents = [];
    icablinkmetricsout.metrics.corr_Pvalue  = [];
    icablinkmetricsout.metrics.conv_Pvalue  = [];
    icablinkmetricsout.metrics.perc_Pvalue  = [];
end

% Put back
EEG.icaact = [];

% Log
EEG.ALSUTRECHT.ica.extra2.musleSlope = single(muscleRatio(:));
EEG.ALSUTRECHT.ica.extra2.blinkPvals = single([icablinkmetricsout.metrics.corr_Pvalue; icablinkmetricsout.metrics.conv_Pvalue; icablinkmetricsout.metrics.perc_Pvalue]');
EEG.ALSUTRECHT.ica.extra2.ICsMostLikelyMuscle = ICsMostLikelyMuscle(:);
EEG.ALSUTRECHT.ica.extra2.icablinkcomps       = icablinkmetricsout.identifiedcomponents(:);
EEG.ALSUTRECHT.ica.extra2.bics = [find(ICsMostLikelyMuscle(:)); icablinkmetricsout.identifiedcomponents(:)];
EEG.ALSUTRECHT.ica.extra2.cvec = [ones(sum(ICsMostLikelyMuscle),1); 2*ones(length(icablinkmetricsout.identifiedcomponents),1)];
EEG.ALSUTRECHT.ica.extra2.clss = {'Muscle','Eye'};

%% 4. Final eyeblink IC detection method
% Based on topoplot correlations

lbls1 = EEG.etc.ic_classification.ICLabel.classes(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics));
lbls2 = EEG.ALSUTRECHT.ica.extra1.clss(EEG.ALSUTRECHT.ica.extra1.cvec);
lbls3 = EEG.ALSUTRECHT.ica.extra2.clss(EEG.ALSUTRECHT.ica.extra2.cvec);
allBadIClbls = [lbls1(:); lbls2(:); lbls3(:)];
allBadICindx = [EEG.ALSUTRECHT.ica.ICLabel.bics(:); EEG.ALSUTRECHT.ica.extra1.bics(:); EEG.ALSUTRECHT.ica.extra2.bics(:)];
allBadICMeth = [ones(length(lbls1(:)),1); 2*ones(length(lbls2(:)),1); 3*ones(length(lbls3(:)),1)];
assert(length(allBadIClbls)==length(allBadICindx));

blinkMask = ismember(allBadIClbls,{'Eye','VEOG'});
allBadICMethTmp = allBadICMeth(blinkMask);
allBadICindxTmp = allBadICindx(blinkMask);
allBadICindxTmpUnq = unique(allBadICindxTmp);

blinksICs = allBadICindxTmpUnq;
% blinksICs = [];
% for i = 1:length(allBadICindxTmpUnq)
%     maskTmp = allBadICindxTmp == allBadICindxTmpUnq(i);
%     if sum(maskTmp)>1
%         blinksICs = [blinksICs, allBadICindxTmpUnq(i)];
%     end
% end

% blinksICs = [];
% for i = 1:length(allBadICindx)
%     if strcmpi(allBadIClbls{i},'VEOG') && any(EEG.ALSUTRECHT.ica.extra2.icablinkcomps == allBadICindx(i))
%         blinksICs = [blinksICs, allBadICindx(i)];
%     end
% end
% if ~isempty(blinksICs), blinksICs = min(blinksICs); end

blinksICsNew  = [];
distBlinksICs = [];

% figure;
if ~isempty(blinksICs)
    corrMat = single(abs(corr(EEG.icawinv,EEG.icawinv(:,blinksICs))));
    [blinksICs, tmp] = find(corrMat>0.8);
    blinksICs = unique(blinksICs);

    % figure; imagesc(corrMat); clim([0 1]); colorbar;
    % colormap(brewermap(128,'PuRd'));

    maskChanBlink = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);

    for i = 1:length(blinksICs)
        W0 = zscore(EEG.icawinv(:,blinksICs(i)));
        W = abs(W0);
        D = max(W(maskChanBlink))-max(W(~maskChanBlink));
        % W = W0; D = abs(mean(W(maskChanBlink))-mean(W(~maskChanBlink)));
        if D>1
            blinksICsNew  = [blinksICsNew, blinksICs(i)];
            distBlinksICs = [distBlinksICs, D];
            % mytopoplot(W0,maskChanBlink,['IC' num2str(blinksICs(i)) ', D = ' num2str(distBlinksICs(end))],nexttile); colorbar;
        end
    end
    corrIC = mean(corrMat(blinksICsNew,:),2);

    % figure; mytopoplot(W,maskChanBlink,'',nexttile); colorbar;
else
    corrMat = NaN;
    corrIC  = [];
end

if any(blinksICsNew)
    fprintf('\nEye ICs (N = %d) were idetified using topoplot correlation.\n',length(blinksICsNew));
else
    fprintf('\nNo eye ICs were idetified using topoplot correlation.\n');
end

% % Fix as it is likely false positive
% EEG.ALSUTRECHT.ica.extra2.bics = [find(ICsMostLikelyMuscle(:)); blinksICsNew(:)];
% EEG.ALSUTRECHT.ica.extra2.cvec = [ones(sum(ICsMostLikelyMuscle),1); 2*ones(length(blinksICsNew),1)];

% Log
EEG.ALSUTRECHT.ica.extra3.corr0 = corrMat;
EEG.ALSUTRECHT.ica.extra3.corr  = corrIC;
EEG.ALSUTRECHT.ica.extra3.dist  = distBlinksICs(:);
EEG.ALSUTRECHT.ica.extra3.bics  = blinksICsNew(:);
EEG.ALSUTRECHT.ica.extra3.cvec  = ones(length(blinksICsNew),1);
EEG.ALSUTRECHT.ica.extra3.clss  = {'Blink'};

%% ========================================================================
% Final log
lbls1 = EEG.etc.ic_classification.ICLabel.classes(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics));
lbls2 = EEG.ALSUTRECHT.ica.extra1.clss(EEG.ALSUTRECHT.ica.extra1.cvec);
lbls3 = EEG.ALSUTRECHT.ica.extra2.clss(EEG.ALSUTRECHT.ica.extra2.cvec);
lbls4 = EEG.ALSUTRECHT.ica.extra3.clss(EEG.ALSUTRECHT.ica.extra3.cvec);

prob1 = EEG.ALSUTRECHT.ica.ICLabel.pvec(EEG.ALSUTRECHT.ica.ICLabel.bics);
prob2 = EEG.ALSUTRECHT.ica.extra1.corr(sub2ind(size(EEG.ALSUTRECHT.ica.extra1.corr), EEG.ALSUTRECHT.ica.extra1.bics, EEG.ALSUTRECHT.ica.extra1.cvec));
prob3Muscle = EEG.ALSUTRECHT.ica.extra2.musleSlope(EEG.ALSUTRECHT.ica.extra2.bics(EEG.ALSUTRECHT.ica.extra2.cvec==1));
prob3Blinks = mean(EEG.ALSUTRECHT.ica.extra2.blinkPvals(EEG.ALSUTRECHT.ica.extra2.bics(EEG.ALSUTRECHT.ica.extra2.cvec==2)),2);
prob3       = [prob3Muscle; prob3Blinks];
prob4       = EEG.ALSUTRECHT.ica.extra3.corr;

meth1  = repmat({'P'},length(prob1),1);
meth2  = repmat({'R'},length(prob2),1);
meth3a = repmat({'Slope'},length(prob3Muscle),1);
meth3b = repmat({'P'},length(prob3Blinks),1);
meth3  = [meth3a; meth3b];
meth4  = repmat({'R'},length(prob4),1);

EEG.ALSUTRECHT.ica.combi.method = [ones(length(EEG.ALSUTRECHT.ica.ICLabel.bics(:)),1); ...
    2*ones(length(EEG.ALSUTRECHT.ica.extra1.bics(:)),1); ...
    3*ones(length(EEG.ALSUTRECHT.ica.extra2.bics(:)),1); ...
    4*ones(length(EEG.ALSUTRECHT.ica.extra3.bics(:)),1)];
EEG.ALSUTRECHT.ica.combi.bics   = [EEG.ALSUTRECHT.ica.ICLabel.bics(:); EEG.ALSUTRECHT.ica.extra1.bics(:); EEG.ALSUTRECHT.ica.extra2.bics(:); EEG.ALSUTRECHT.ica.extra3.bics(:)];
EEG.ALSUTRECHT.ica.combi.prbs   = [prob1(:); prob2(:); prob3(:); prob4(:)];
EEG.ALSUTRECHT.ica.combi.lbls   = [lbls1(:); lbls2(:); lbls3(:); lbls4(:)];
EEG.ALSUTRECHT.ica.combi.meth   = [meth1(:); meth2(:); meth3(:); meth4(:)];

% Labels:
% 1 'Brain'
% 2 'Muscle'
% 3 'Eye'
% 4 'Heart'
% 5 'Line Noise'
% 6 'Channel Noise'
% 7 'Other'

EEG.ALSUTRECHT.ica.combi.report = EEG.ALSUTRECHT.ica.ICLabel.cvec;
EEG.ALSUTRECHT.ica.combi.report(ICsMostLikelyHeart)                                                = 4; % ECG       -> (4) Heart
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra2.ICsMostLikelyMuscle)                     = 2; % Slope     -> (2) Muscle
% EEG.ALSUTRECHT.ica.combi.report(ICsMostLikelyBlink)                                              = 3; % Blinks-> (3) Eye
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra1.bics(EEG.ALSUTRECHT.ica.extra1.cvec>=2)) = 3; % VEOG/HEOG -> (3) Eye
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra2.bics(EEG.ALSUTRECHT.ica.extra2.cvec==2)) = 3; % Blinks    -> (3) Eye
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra3.bics(EEG.ALSUTRECHT.ica.extra3.cvec==1)) = 3; % Blinks    -> (3) Eye

assert(length(EEG.ALSUTRECHT.ica.combi.bics)==length(EEG.ALSUTRECHT.ica.combi.prbs));
assert(length(EEG.ALSUTRECHT.ica.combi.bics)==length(EEG.ALSUTRECHT.ica.combi.lbls));

%% ========================================================================
% Only blinks ICs
ICsMostLikelyBlink = false(NICA,1);
ICsMostLikelyBlink(EEG.ALSUTRECHT.ica.extra3.bics) = true;

% Only EMG ICs
ICsMostLikelyMuscle = EEG.ALSUTRECHT.ica.extra2.ICsMostLikelyMuscle;

% Remove these as they are likely complex artifact ICs
ICsMostLikelyComplex = ICsMostLikelyMuscle&ICsMostLikelyBlink;

% TMP = EEG.etc.ic_classification.ICLabel.classifications(ICsMostLikelyBlink,2);
% for i = 1:20
%     pdN = fitdist(ICAdata(i,:)','normal');
%     pdT = fitdist(ICAdata(i,:)','tLocationScale'); % tLocationScale / Stable
%
%     % figure; hold on;
%     % hdata = histogram(ICAdata(i,:),'binwidth',0.1,'Normalization','pdf');
%     [n,edges,bin] = histcounts(ICAdata(i,:),'binwidth',0.5,'Normalization','pdf');
%     hdata = [];
%     hdata.Values = n;
%     hdata.BinLimits = edges([1 end]);
%
%     % myLim = round([min(ICAdata(1,:))-5, max(ICAdata(1,:))+5]);
%     xgrid = linspace(hdata.BinLimits(1),hdata.BinLimits(2),length(hdata.Values));
%     pdfEstN = pdf(pdN,xgrid);
%     pdfEstT = pdf(pdT,xgrid);
%
%     % plot(xgrid,pdfEstN);
%     % plot(xgrid,pdfEstT);
%     % legend({'Data','Normal','Student'});
%
%     RsqN = get_Rsq(hdata.Values,pdfEstN);
%     RsqT = get_Rsq(hdata.Values,pdfEstT);
%     if RsqT>RsqN, disp(i);  end
% end

% Other ICs for wICA
ICsforwICA = false(NICA,1);
ICsforwICA(unique(EEG.ALSUTRECHT.ica.combi.bics)) = true;
ICsforwICA(ICsMostLikelyMuscle|ICsMostLikelyComplex) = false;

% Log
EEG.ALSUTRECHT.ica.ICsMostLikelyBlink   = ICsMostLikelyBlink;
EEG.ALSUTRECHT.ica.ICsMostLikelyMuscle  = ICsMostLikelyMuscle;
EEG.ALSUTRECHT.ica.ICsMostLikelyComplex = ICsMostLikelyComplex;
EEG.ALSUTRECHT.ica.ICsforwICA           = ICsforwICA;

end