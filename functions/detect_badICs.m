function EEG = detect_badICs(EEG,EXT,cfg)

fprintf('\nDetecting bad ICs...\n\n');

% Double-check
EEG = eeg_checkset(EEG,'ica');

% Make sure IC activations are present
EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);

NICA = size(EEG.icaact,1);

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
maskExt  = [chanecg, chanveog, chanheog];

dataExt  = EXT.data(maskExt,:);
extLabels = {'ECG','VEOG','HEOG'};
NEXT = length(extLabels);
% extLabels = {EXT.chanlocs(maskExt).labels};
% disp(extLabels);

% Temporarily filter for better detection
% It is fine that it will be now double-filtered
[bl, al] = butter(4,30/(EEG.srate/2),'low');
[bh, ah] = butter(4,2/(EEG.srate/2),'high');
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

% ECG was not recorded
fprintf('\n');
if isempty(chanecg)
    warning('ECG was not recorded!'); % disp(extLabels);
    EXTdata = [NaN*EXTdata(:,1), EXTdata];
else
    % wt = modwt(EXTdata(:,1),5);
    % wtrec = zeros(size(wt));
    % wtrec(4:5,:) = wt(4:5,:);
    % ECG = imodwt(wtrec,'sym4');
    %
    % wt = modwt(ICAdata,5);
    % wtrec = zeros(size(wt));
    % wtrec(4:5,:,:) = wt(4:5,:,:);
    % ICAECG = imodwt(wtrec,'sym4');
    %
    % ECG = abs(ECG).^2;
    % ICAECG = abs(ICAECG).^2;
    %
    % corrECG = abs(corr(ECG',ICAECG));
end

% Correlation
corrECG = abs(corr(ICAdata.^2,abs(EXTdata(:,1))));
corrMat = abs(corr(ICAdata,EXTdata));
corrMat(:,1) = corrECG;

% Tresholds: ECG / VEOG / HEOG
corrTreshold = [0.3 0.6 0.6];
[badIC, badICtype] = find(corrMat>corrTreshold);

% figure; imagesc(corrMat); clim([0 1]); colorbar;
% xticks([1 2 3]); xticklabels(extLabels); colormap(brewermap(128,'PuRd'));
% EEG = pop_saveset(EEG,'filename',['TMP.set'],'filepath',EEG.ALSUTRECHT.subject.preproc);

% Prvent false positive HEOG
% Here, HEOG ICs are easily confused by broad L-R dipolar (brain) ICs
% Remove HEOG detections if they are likely brain ICs
maskHEOGIC = badIC(badICtype==3);
falseHEOGIC = EEG.ALSUTRECHT.ica.ICLabel.cvec(maskHEOGIC) == 1 & EEG.ALSUTRECHT.ica.ICLabel.pvec(maskHEOGIC) > 0.6;

badIC(falseHEOGIC) = [];
badICtype(falseHEOGIC) = [];

% Report
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
options.muscleFreqIn    = [7, 70];
options.Freq_to_compute = [1, 100];

% Calculate pwelch to enable detection of log-freq log-power slopes,
% indicative of muscle activity
eegData = EEG.icaact(:,:);

[pxx,fp] = pwelch(eegData',size(eegData,2),[],size(eegData,2),EEG.srate);
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
try
    % EOGdata = mean(EEG.data(ismember({EEG.chanlocs.labels},cfg.ica.blinkchans),:),1);
    EOGdata = dataExt(2,:);
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
%     pdN = fitdist(EEG.icaact(i,:)','normal');
%     pdT = fitdist(EEG.icaact(i,:)','tLocationScale'); % tLocationScale / Stable
%
%     % figure; hold on;
%     % hdata = histogram(EEG.icaact(i,:),'binwidth',0.1,'Normalization','pdf');
%     [n,edges,bin] = histcounts(EEG.icaact(i,:),'binwidth',0.5,'Normalization','pdf');
%     hdata = [];
%     hdata.Values = n;
%     hdata.BinLimits = edges([1 end]);
%
%     % myLim = round([min(EEG.icaact(1,:))-5, max(EEG.icaact(1,:))+5]);
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