function EEG = detect_badICs(EEG,EXT,cfg)

%% ========================================================================
% ICLabel
EEG = iclabel(EEG);

% Flag artifact ICs
EEG = pop_icflag(EEG,cfg.ica.iclabel);

% Log IClabel info
EEG.ALSUTRECHT.ica.ICLabel.bics = find(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.ICLabel.clss = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel.pvec, EEG.ALSUTRECHT.ica.ICLabel.cvec] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);

% % Manually add heart ICs
% K = find(contains(EEG.ALSUTRECHT.ica.ICLabel_clss,'Heart'));
% ICHeart = EEG.ALSUTRECHT.ica.ICLabel_cvec==K;
% if any(ICHeart)
%     if EEG.ALSUTRECHT.ica.ICLabel_pvec(ICHeart)>0.5
%         EEG.reject.gcompreject(ICHeart) = true;
%     end
% end

%% ========================================================================
% Make sure IC activations are present
if isempty(EEG.icaact)
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end

chanveog = strcmp({EXT.chanlocs.labels},'VEOG');
chanheog = strcmp({EXT.chanlocs.labels},'HEOG');
chanecg  = strcmp({EXT.chanlocs.labels},'ECG'); % if it was recorded!
maskExt  = logical(chanveog+chanheog+chanecg);
dataExt  = EXT.data(maskExt,:);
% disp({EXT.chanlocs(maskExt).labels});
extLabels = {EXT.chanlocs(maskExt).labels};

% Temporarily filter for better detection
% It is fine that it will be now double-filtered
[bl, al] = butter(2,25/(EEG.srate/2),'low');
[bh, ah] = butter(2,1/(EEG.srate/2),'high');
assert(isstable(bl,al));
assert(isstable(bh,ah));

ICAdata = filtfilt(bl,al,EEG.icaact');
ICAdata = filtfilt(bh,ah,ICAdata);
EXTdata = filtfilt(bl,al,dataExt');
EXTdata = filtfilt(bh,ah,EXTdata);

% ICAdata = robust_zscore(ICAdata);
% EXTdata = robust_zscore(EXTdata);
% ICAdata = zscore(ICAdata);
% EXTdata = zscore(EXTdata);

corrMat = abs(corr(ICAdata,EXTdata));

% figure; imagesc(corrMat); clim([0 1]); colorbar;
% xticks([1 2 3]); xticklabels(extLabels); colormap(brewermap(128,'PuRd'));

% ECG / VEOG / HEOG
if length(extLabels)==3
    [badIC, badICtype] = find(corrMat>[0.2 0.2 0.9]);
else
    % VEOG / HEOG
    assert(sum(chanveog|chanheog)==2); assert(sum(chanecg)==0);
    [badIC, badICtype] = find(corrMat>[0.2 0.9]);
end

% Log
EEG.ALSUTRECHT.ica.extra1.corr = single(corrMat);
EEG.ALSUTRECHT.ica.extra1.bics = badIC;
EEG.ALSUTRECHT.ica.extra1.cvec = badICtype;
EEG.ALSUTRECHT.ica.extra1.clss = extLabels;

%% ========================================================================
% 1. Detect EMG ICs using freq slopes
options.muscleFreqIn    = [7, 70];
options.Freq_to_compute = [1, 100];

NICA = size(EEG.icaact,1);

% Resize EEG.icaact if required
if size(EEG.icaact,3) > 0
    eegData = reshape(EEG.icaact,NICA,[]);
else
    eegData = EEG.icaact;
end

% Calculate pwelch to enable detection of log-freq log-power slopes,
% indicative of muscle activity
[pxx,fp] = pwelch(eegData',size(eegData,2),[],size(eegData,2),EEG.srate);
FFTout = pxx';
fp = fp';

% Calculate FFT bins
freq = options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
fftBins = zeros(size(FFTout,1),size(freq,2)); % preallocate
for a=1:size(freq,2)
    [~, index1]=min(abs(fp-((freq(1,a)-0.25))));
    [~, index2]=min(abs(fp-((freq(1,a)+0.25))));
    fftBins(:,a) = mean(FFTout(:,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
end

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

ICsMostLikelyMuscle = muscleRatio>=cfg.bch.muscleSlopeThreshold;
% ICsMostLikelyMuscle = muscleRatio>=-0.5;
% ICsMostLikelyMuscle = (muscle_ICs==1);

% 2. Use icablinkmetrics for eyeblinks
if exist('icablinkmetrics','file') == 2
    try
        % icablinkmetricsout = icablinkmetrics(EEG0,'ArtifactChannel',EEG0.data(strcmp({EEG0.chanlocs.labels},'VEOG'),:),'Alpha',0.001,'VisualizeData','False');
        icablinkmetricsout = icablinkmetrics(EEG,'ArtifactChannel',mean(EEG.data(ismember({EEG.chanlocs.labels},cfg.ica.blinkchans),:),1),'Alpha',0.001,'VisualizeData','False');
        if any(icablinkmetricsout.identifiedcomponents>0)
            fprintf('icablinkmetrics has identified %d eye component(s).\n',length(icablinkmetricsout.identifiedcomponents));
            % ICsMostLikelyNotBrain(icablinkmetricsout.identifiedcomponents) = true;
            % ICsMostLikelyEye(icablinkmetricsout.identifiedcomponents)      = true;
        else
            fprintf('icablinkmetrics has not identified any eye components.\n');
            icablinkmetricsout.identifiedcomponents = [];
        end
    catch
        warning('icablinkmetrics has failed...');
    end
end

% Log
EEG.ALSUTRECHT.ica.extra2.musleSlope = single(muscleRatio(:));
EEG.ALSUTRECHT.ica.extra2.blinkPvals = single([icablinkmetricsout.metrics.corr_Pvalue; icablinkmetricsout.metrics.conv_Pvalue; icablinkmetricsout.metrics.perc_Pvalue]');
EEG.ALSUTRECHT.ica.extra2.ICsMostLikelyMuscle = ICsMostLikelyMuscle(:);
EEG.ALSUTRECHT.ica.extra2.icablinkcomps       = icablinkmetricsout.identifiedcomponents(:);
EEG.ALSUTRECHT.ica.extra2.bics = [find(ICsMostLikelyMuscle(:)); icablinkmetricsout.identifiedcomponents(:)];
EEG.ALSUTRECHT.ica.extra2.cvec = [ones(sum(ICsMostLikelyMuscle),1); 2*ones(length(icablinkmetricsout.identifiedcomponents),1)];
EEG.ALSUTRECHT.ica.extra2.clss = {'Muscle','Eye'};

%% Final eyeblink IC detection method
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

if ~isempty(blinksICs)
    corrMat = single(abs(corr(EEG.icawinv,EEG.icawinv(:,blinksICs))));
    [blinksICs, tmp] = find(corrMat>0.75);
    blinksICs = unique(blinksICs);

    % figure; imagesc(corrMat); clim([0 1]); colorbar;
    % colormap(brewermap(128,'PuRd'));

    maskChanBlink = ismember({EEG.chanlocs(:).labels},cfg.ica.blinkchans);

    for i = 1:length(blinksICs)
        W = abs(zscore(EEG.icawinv(:,blinksICs(i))));
        D = max(W(maskChanBlink))-max(W(~maskChanBlink));
        if D>0.5
            blinksICsNew  = [blinksICsNew, blinksICs(i)];
            distBlinksICs = [distBlinksICs, D];
        end
    end

    % figure; mytopoplot(W,maskChanBlink,'',nexttile); colorbar;
else
    corrMat = NaN;
end

% Fix as it is likely false positive
EEG.ALSUTRECHT.ica.extra2.bics = [find(ICsMostLikelyMuscle(:)); blinksICsNew(:)];
EEG.ALSUTRECHT.ica.extra2.cvec = [ones(sum(ICsMostLikelyMuscle),1); 2*ones(length(blinksICsNew),1)];

% Log
EEG.ALSUTRECHT.ica.extra3.corr = corrMat;
EEG.ALSUTRECHT.ica.extra3.dist = distBlinksICs(:);
EEG.ALSUTRECHT.ica.extra3.bics = blinksICsNew(:);
EEG.ALSUTRECHT.ica.extra3.cvec = ones(length(blinksICsNew),1);
EEG.ALSUTRECHT.ica.extra3.clss = {'Blink'};

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
prob4 = EEG.ALSUTRECHT.ica.extra3.corr(EEG.ALSUTRECHT.ica.extra3.bics);

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
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra1.bics(EEG.ALSUTRECHT.ica.extra1.cvec==1)) = 4; % ECG       -> (4) Heart
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra1.bics(EEG.ALSUTRECHT.ica.extra1.cvec>=2)) = 3; % VEOG/HEOG -> (3) Eye
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra2.bics(EEG.ALSUTRECHT.ica.extra2.cvec==1)) = 2; % Slope     -> (2) Muscle
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra2.bics(EEG.ALSUTRECHT.ica.extra2.cvec==2)) = 3; % Blinks    -> (3) Eye
EEG.ALSUTRECHT.ica.combi.report(EEG.ALSUTRECHT.ica.extra3.bics(EEG.ALSUTRECHT.ica.extra3.cvec==1)) = 3; % Blinks    -> (3) Eye

assert(length(EEG.ALSUTRECHT.ica.combi.bics)==length(EEG.ALSUTRECHT.ica.combi.prbs));
assert(length(EEG.ALSUTRECHT.ica.combi.bics)==length(EEG.ALSUTRECHT.ica.combi.lbls));

end