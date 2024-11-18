function EEG = detect_badcomponents2(EEG,EXT,cfg)

fprintf('\nDetecting bad ICs...\n\n');

% Double-check
EEG = eeg_checkset(EEG,'ica');

% Make sure IC activations are present
EEG.icaact = [];
ICAdata    = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
NICA       = size(ICAdata,1);

%% ========================================================================
% 1. ICLabel
fprintf('\n================================\n');
fprintf('ICLabel\n');
fprintf('================================\n');

EEG = iclabel(EEG);

% Flag artifact ICs
EEG = pop_icflag(EEG,cfg.ica.iclabel);

% Log IClabel info
EEG.ALSUTRECHT.ica.ICLabel.bics = find(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.ICLabel.clss = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel.pvec, EEG.ALSUTRECHT.ica.ICLabel.cvec] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);

%% ========================================================================
% 2. Correlation with external channels (ECG / VEOG / HEOG)
fprintf('\n================================\n');
fprintf('ECG/VEOG/HEOG: Correlations with the external electrodes\n');
fprintf('================================\n');

chanecg  = find(strcmp({EXT.chanlocs.labels},'ECG')); % if it was recorded!
chanveog = find(strcmp({EXT.chanlocs.labels},'VEOG'));
chanheog = find(strcmp({EXT.chanlocs.labels},'HEOG'));
% maskExt  = [chanecg, chanveog, chanheog];

% If the signal is not recorded, try to estimate it using
% an approximate reconstruciton (ICA-like) weights
if isempty(chanecg)
    warning('ECG signals was not recorded. Some heart IC check cannot be done.');
    EXTdataECG = []; EXTTMP = [];

    % warning('ECG signals was not recorded, but we can try to estimate it from the EEG data...');
    % load(fullfile(EEG.ALSUTRECHT.subject.mycodes,'files','Heartweights'),'Heartweights');
    % EXTdataECG = Heartweights'*EEG.data(:,:);
    %
    % EXTTMP = EXT;
    % EXTTMP.nbchan                         = EXTTMP.nbchan+1;
    % EXTTMP.data(EXTTMP.nbchan,:)          = EXTdataECG;
    % EXTTMP.chanlocs(EXTTMP.nbchan).labels = 'ECG';
    % EXTTMP.chanlocs(EXTTMP.nbchan).type   = 'EXT';
    % ECGsignalLabel = 'Estimated';
else
    EXTdataECG = EXT.data(chanecg,:);

    EXTTMP = EXT;
    ECGsignalLabel = 'Recorded';
end

EXTdataEOG = EXT.data([chanveog chanheog],:);
extLabels  = {'ECG','VEOG','HEOG'};

% Temporarily filter for better detection
% It is fine that it will be double-filtered
% EOG: 1-10 Hz
% ECG: 8-16 Hz / 10-20 Hz
[blEOG, alEOG] = butter(4,10/(EEG.srate/2),'low');
[blECG, alECG] = butter(4,30/(EEG.srate/2),'low');
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
if ~isempty(EXTdataECG)
    EXTdataECG = do_filteringcore(blECG,alECG,EXTdataECG,EEG.event,EEG.srate);
    EXTdataECG = do_filteringcore(bhECG,ahECG,EXTdataECG,EEG.event,EEG.srate);
else
    EXTdataECG = NaN*EXTdataEOG(1,:);
end

% Combine: ECG / VEOG / HEOG
EXTdata = [EXTdataECG; EXTdataEOG]';

% Correlation
% 1. ECG
corrECG = corr(ICAdataECG.^2,abs(EXTdata(:,1)));

% % Targeted
% [ECGmask, ECGepoch] = detect_ecg(EXTTMP,250,ECGsignalLabel);
% if ~isnan(ECGmask)
%     Nhbeats = size(ECGepoch,1);
%     corrECG = NaN(NICA,Nhbeats);
%
%     for i = 1:Nhbeats
%         thisChunk = ECGepoch(i,1):ECGepoch(i,2);
%         % corrECG(:,i) = corr(ICAdataECG(thisChunk,:).^2,abs(EXTdata(thisChunk,1)));
%         corrECG(:,i) = corr(ICAdataECG(thisChunk,:),EXTdata(thisChunk,1));
%     end
%     % figure; imagesc(corrECG);
%     corrECG = mean(corrECG,2);
% else
%     corrECG = NaN(NICA,1);
% end

% 2. VEOG
corrVEOG = corr(ICAdataEOG,EXTdata(:,2));

% % Targeted
% [VEOGmask, VEOGepochs] = detect_veog(EXT,400,false);
%
% Nblinks = size(VEOGepochs,1);
% corrVEOG = NaN(NICA,Nblinks);
%
% for i = 1:Nblinks
%     thisChunk = VEOGepochs(i,1):VEOGepochs(i,2);
%     corrVEOG(:,i) = corr(ICAdataEOG(thisChunk,:),EXTdata(thisChunk,2));
% end
% % figure; imagesc(corrVEOG);
% corrVEOG = mean(corrVEOG,2);

% 3. HEOG
corrHEOG = corr(ICAdataEOG,EXTdata(:,3));

% Combine: ECG / VEOG / HEOG
corrMat = [corrECG, corrVEOG, corrHEOG];

% Tresholds: ECG / VEOG / HEOG
% R
corrTreshold = [0.4 0.6 0.4];
% Zscore
% corrTreshold = [3 3 3];

for i = 1:3
    switch i
        case 1
            % ECG, only 1 IC
            corrTmp = abs(corrMat(:,i));
            % corrTmp = abs(zscore(corrMat(:,i)));
            badICECG = find(corrTmp>corrTreshold(i));

            if length(badICECG)>1
                [~, badICTmpTmp] = max(corrTmp(badICECG));
                badICECG = badICECG(badICTmpTmp);
            end

        case 2
            % VEOG, can be multiple
            corrTmp = abs(corrMat(:,i));
            % corrTmp = abs(zscore(corrMat(:,i)));
            badICVEOG = find(corrTmp>corrTreshold(i));

            % % VEOG, iterative removal
            % badICVEOG = [];  % Store indices of bad components
            %
            % corrTmp = corrMat(:,i);  % Initialize the correlation data for this iteration
            % NICA = length(corrTmp);  % Get the length of the data
            %
            % flagGo = true;
            % while flagGo
            %     % Exclude previously identified bad components
            %     indx = 1:NICA;
            %     indx(badICVEOG) = [];  % Remove bad ICs from index list
            %
            %     % Calculate the Z-score for remaining components
            %     corrTmpZ = corrTmp(indx);  % Remove bad ICs from the correlation data
            %     corrTmpZ = (corrTmpZ - mean(corrTmpZ)) ./ std(corrTmpZ);  % Z-score of remaining data
            %
            %     % Find any components that exceed the threshold
            %     badICTmpTmp = indx(abs(corrTmpZ) > corrTreshold(i));  % Get indices of bad components in original data
            %
            %     % If we found new bad components, add them to the list
            %     if ~isempty(badICTmpTmp)
            %         badICVEOG = [badICVEOG badICTmpTmp];  % Add the newly found bad ICs to the total list
            %     else
            %         flagGo = false;  % No more components exceed the threshold, exit loop
            %     end
            % end

        case 3
            % HEOG, only 1 IC
            corrTmp = abs(corrMat(:,i));
            % corrTmp = abs(zscore(corrMat(:,i)));
            badICHEOG = find(corrTmp>corrTreshold(i));

            if length(badICHEOG)>1
                [maxValue, badICTmpTmp] = max(corrTmp(badICHEOG));
                badICHEOG = badICHEOG(badICTmpTmp);
            end

            if ~isempty(badICHEOG)
                % Prevent false positive HEOG
                % 1. Here, HEOG ICs are easily confused by broad L-R dipolar (brain) ICs
                falseHEOGIC1 = EEG.ALSUTRECHT.ica.ICLabel.cvec(badICHEOG) == 1 & EEG.ALSUTRECHT.ica.ICLabel.pvec(badICHEOG) > 0.6;

                % 2. VEOG ICs are sometimes very correlated with the HEOG signal
                load(fullfile(EEG.ALSUTRECHT.subject.mycodes,'files','Blinkweights'),'Blinkweights'); % template
                falseHEOGIC2 = abs(corr(EEG.icawinv(:,badICHEOG),Blinkweights,"type","Spearman")) > 0.85;

                % Either/or should be true if false positive
                badICHEOG(falseHEOGIC1 | falseHEOGIC2) = [];
            end
    end
end

% Combine: ECG / VEOG / HEOG
badIC = [badICECG(:); badICVEOG(:); badICHEOG(:)];
badICtype = [ones(length(badICECG),1); 2*ones(length(badICVEOG),1); 3*ones(length(badICHEOG),1)];

% Report
NEXT = length(extLabels);
for i = 1:NEXT
    if any(badICtype==i)
        fprintf('%s ICs (N = %d) were identified using EXT channel correlation (Z>%1.1f).\n',extLabels{i},sum(badICtype==i),corrTreshold(i));
        str = strjoin(string(badIC(badICtype==i)),', ');
        fprintf('%s ICs: %s\n',extLabels{i},str);
    else
        fprintf('No %s ICs were identified using EXT channel correlation (Z>%1.1f).\n',extLabels{i},corrTreshold(i));
    end
end

% Log
EEG.ALSUTRECHT.ica.corr.corr = single(corrMat);
EEG.ALSUTRECHT.ica.corr.bics = badIC;
EEG.ALSUTRECHT.ica.corr.cvec = badICtype;
EEG.ALSUTRECHT.ica.corr.clss = extLabels;

%% ========================================================================
% 3. Cross-trial phase statistics for ECG detection
% https://mne.tools/stable/generated/mne.preprocessing.ICA.html#mne.preprocessing.ICA.find_bads_ecg
% https://sci-hub.st/https://ieeexplore.ieee.org/document/4536072
fprintf('\n================================\n');
fprintf('ECG: Cross-trial phase statistics\n');
fprintf('================================\n');

if ~isempty(EXTTMP)
    % Detect ECG/QRS
    [ECGmask, ECGepoch, ~, ~, plusEstimate] = detect_ecg(EXTTMP,[-250 450],ECGsignalLabel);

    if ~isnan(ECGmask)
        [blECG, alECG] = butter(4,16/(EEG.srate/2),'low');
        [bhECG, ahECG] = butter(4,8/(EEG.srate/2),'high');
        ICAdataECG = do_filteringcore(blECG,alECG,ICAdata,EEG.event,EEG.srate);
        ICAdataECG = do_filteringcore(bhECG,ahECG,ICAdataECG,EEG.event,EEG.srate);

        [V, pK] = my_ctps(ICAdataECG, ECGepoch);
        % figure; plot(V);
        V = zscore(V);

        ctpsTreshold1 = 5;  % 0.3-0.4?
        ctpsTreshold2 = 20; % 20 original paper
        badIC1 = find(V>ctpsTreshold1);
        badIC2 = find(pK>=ctpsTreshold2);

        % Sometimes there are weird false positives
        load(fullfile(EEG.ALSUTRECHT.subject.mycodes,'files','Heartweights'),'Heartweights'); % template
        if any(badIC1)
            falseECGIC = abs(corr(EEG.icawinv(:,badIC1),Heartweights,"type","Spearman")) < 0.85;
            badIC1(falseECGIC) = [];
        end
        if any(badIC2)
            falseECGIC = abs(corr(EEG.icawinv(:,badIC2),Heartweights,"type","Spearman")) < 0.85;
            badIC2(falseECGIC) = [];
        end

        % Combine
        badIC = unique([badIC1(:); badIC2(:)]);

        if ~isempty(badIC1)
            fprintf('ECG ICs (N = %d) were identified using cross-trial phase statistics (Z>%1.1f).\n',length(badIC1),ctpsTreshold1);
            str = strjoin(string(badIC1),', ');
            fprintf('ECG ICs: %s\n',str);
        else
            fprintf('No ECG ICs were identified using cross-trial phase statistics (Z>%1.1f).\n',ctpsTreshold1);
        end
        if ~isempty(badIC2)
            fprintf('ECG ICs (N = %d) were identified using cross-trial phase statistics (pK>=%1.1f).\n',length(badIC2),ctpsTreshold2);
            str = strjoin(string(badIC2),', ');
            fprintf('ECG ICs: %s\n',str);
        else
            fprintf('No ECG ICs were identified using cross-trial phase statistics (pK>=%1.1f).\n',ctpsTreshold2);
        end

    else
        fprintf('Skipping as the estimated/recorded ECG signal is noisy.\n');
        V = NaN(NICA,1);
        badIC = [];
    end
else
    fprintf('Skipping as the ECG signal is not recorded.\n');
    V = NaN(NICA,1);
    badIC = [];
    plusEstimate = NaN;
end

% Log
EEG.ALSUTRECHT.ica.ctps.corr = V;
EEG.ALSUTRECHT.ica.ctps.bics = badIC;
EEG.ALSUTRECHT.ica.ctps.cvec = ones(length(EEG.ALSUTRECHT.ica.ctps.bics),1);
EEG.ALSUTRECHT.ica.ctps.clss = 'ECG';

EEG.ALSUTRECHT.subject.plusEstimate = plusEstimate;

%% ========================================================================
% 4. Detect EMG ICs using freq slopes
fprintf('\n================================\n');
fprintf('EMG: Power slopes\n');
fprintf('================================\n');

% options.muscleFreqIn       = [7, 70];
% options.muscleSlopeThreshold = -0.5;
options.muscleFreqIn         = cfg.bch.muscleSlopeFreq;
options.muscleSlopeThreshold = cfg.bch.muscleSlopeThreshold;
options.Freq_to_compute      = [1 100];
options.muscleFreqEx         = 50 + 2*[-1 1]; % Line freq +-bandwith

% Calculate pwelch to enable detection of log-freq log-power slopes, indicative of muscle activity
[pow,frefAll] = pwelch(ICAdata',size(ICAdata,2),[],size(ICAdata,2),EEG.srate);
pow = pow';
frefAll = frefAll';

% Calculate FFT bins
freq = options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
fftBins = zeros(size(pow,1),size(freq,2)); % preallocate
for a = 1:size(freq,2)
    [~, index1] = min(abs(frefAll-((freq(1,a)-0.25))));
    [~, index2] = min(abs(frefAll-((freq(1,a)+0.25))));
    fftBins(:,a) = mean(pow(:,index1:index2),2); % creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
end

% figure;
% nexttile; plot(fp,pxx(:,1));
% nexttile; plot(freq,fftBins(1,:));

% Better muscle comp_number identification:
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
badIC = find(muscleRatio > options.muscleSlopeThreshold);

if ~isempty(badIC)
    fprintf('Muscle ICs (N = %d) were identified using ICLabel and power slope (>%1.1f).\n',length(badIC),options.muscleSlopeThreshold);
    str = strjoin(string(badIC),', ');
    fprintf('EMG ICs: %s\n',str);
else
    fprintf('No muscle ICs were identified using ICLabel and power slope (>%1.1f).\n',options.muscleSlopeThreshold);
end

% Log
EEG.ALSUTRECHT.ica.emg.slope = muscleRatio(:);
EEG.ALSUTRECHT.ica.emg.bics  = badIC(:);
EEG.ALSUTRECHT.ica.emg.cvec  = ones(sum(EEG.ALSUTRECHT.ica.emg.bics),1);
EEG.ALSUTRECHT.ica.emg.clss  = 'EMG';

%% ========================================================================
% 5. Use icablinkmetrics for eyeblinks
fprintf('\n================================\n');
fprintf('VEOG: icablinkmetrics plugin\n');
fprintf('================================\n');

% Blink / VEOG channel
% EOGdata = mean(EEG.data(ismember({EEG.chanlocs.labels},cfg.ica.blinkchans),:),1);
EOGdata = EXTdata(:,2)';

% Put bandpassed ICA data
EEG.icaact = ICAdataEOG';

try
    icablinkmetricsout = icablinkmetrics(EEG,'ArtifactChannel',EOGdata,'Alpha',0.001,'VisualizeData','False');
    if any(icablinkmetricsout.identifiedcomponents>0)
        fprintf('Blink ICs (N = %d) were identified using icablinkmetrics.\n',length(icablinkmetricsout.identifiedcomponents));
        str = strjoin(string(icablinkmetricsout.identifiedcomponents),', ');
        fprintf('VEOG ICs: %s\n',str);
    else
        fprintf('No blink ICs were identified using icablinkmetrics.\n');
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
EEG.ALSUTRECHT.ica.icablinkmetrics.pval = single([icablinkmetricsout.metrics.corr_Pvalue; icablinkmetricsout.metrics.conv_Pvalue; icablinkmetricsout.metrics.perc_Pvalue]');
EEG.ALSUTRECHT.ica.icablinkmetrics.bics = icablinkmetricsout.identifiedcomponents(:);
EEG.ALSUTRECHT.ica.icablinkmetrics.cvec = ones(length(EEG.ALSUTRECHT.ica.icablinkmetrics.bics),1);
EEG.ALSUTRECHT.ica.icablinkmetrics.clss = 'Blink';

%% ========================================================================
% 6. Use spatial characteristics for channel pops and wobbles
fprintf('\n================================\n');
fprintf('Channel: Spatial characteristics\n');
fprintf('================================\n');

% Estimate spatial smootheness
[spatialSmoothness, badIC, spatialTreshold] = estimate_spatialsmoothnes(EEG);

if ~isempty(badIC)
    fprintf('Channel ICs (N = %d) were identified using spatial characteristics (Var1<%1.1f & Var2>%d).\n',length(badIC),spatialTreshold);
    str = strjoin(string(badIC),', ');
    fprintf('Channel ICs: %s\n',str);
else
    fprintf('No channel ICs were identified using spatial characteristics (Var1<%1.1f & Var2>%d).\n',spatialTreshold);
end

% Log
EEG.ALSUTRECHT.ica.spatialSmoothness.pval = spatialSmoothness;
EEG.ALSUTRECHT.ica.spatialSmoothness.bics = badIC(:);
EEG.ALSUTRECHT.ica.spatialSmoothness.cvec = ones(length(EEG.ALSUTRECHT.ica.spatialSmoothness.bics),1);
EEG.ALSUTRECHT.ica.spatialSmoothness.clss = 'Channel';

%% ========================================================================
% 7. Use blink IC template for blinks
fprintf('\n================================\n');
fprintf('VEOG: Blink IC template\n');
fprintf('================================\n');

% Load the template
load(fullfile(EEG.ALSUTRECHT.subject.mycodes,'files','Blinkweights'),'Blinkweights');

% Smoothen the ICs to remove 'noise' and improve the correlation
icawinvSmooth = estimate_invlaplacian(EEG.icawinv,EEG.chanlocs,1);

% Lower to capture imperfect ICs, but leads to false positives
corrMat = abs(corr(icawinvSmooth,Blinkweights,"type","Spearman"));
blinkTreshold = 0.95;
badIC = find(corrMat>blinkTreshold);

% Prevent false positive VEOG
% Here, VEOG ICs are confused by broad A-P dipolar (brain) ICs
falseVEOGIC = EEG.ALSUTRECHT.ica.ICLabel.cvec(badIC) == 1 & EEG.ALSUTRECHT.ica.ICLabel.pvec(badIC) > 0.5;
badIC(falseVEOGIC) = [];

if ~isempty(badIC)
    fprintf('VEOG ICs (N = %d) were identified using template correlation (R>%1.2f).\n',length(badIC),blinkTreshold);
    str = strjoin(string(badIC),', ');
    fprintf('VEOG ICs: %s\n',str);
else
    fprintf('No blink ICs were identified using template correlation (R>%1.2f).\n',blinkTreshold);
end

% Log
EEG.ALSUTRECHT.ica.blinkTemplateCorr.pval = corrMat;
EEG.ALSUTRECHT.ica.blinkTemplateCorr.bics = badIC(:);
EEG.ALSUTRECHT.ica.blinkTemplateCorr.cvec = ones(length(EEG.ALSUTRECHT.ica.blinkTemplateCorr.bics),1);
EEG.ALSUTRECHT.ica.blinkTemplateCorr.clss = 'Blink';

%% ========================================================================
% 8. Use Saccade IC template for Saccade
fprintf('\n================================\n');
fprintf('HEOG: Saccade IC template\n');
fprintf('================================\n');

% Load the template
load(fullfile(EEG.ALSUTRECHT.subject.mycodes,'files','Saccadeweights'),'Saccadeweights');

% Lower to capture imperfect ICs, but leads to false positives
corrMat = abs(corr(icawinvSmooth,Saccadeweights,"type","Pearson"));
saccadeTreshold = 0.95;
badIC = find(corrMat>saccadeTreshold);

% figure; 
% mytopoplot(Saccadeweights,[],[],nexttile);
% mytopoplot(icawinvSmooth(:,18),[],[],nexttile);

if ~isempty(badIC)
    fprintf('HEOG ICs (N = %d) were identified using template correlation (R>%1.2f).\n',length(badIC),saccadeTreshold);
    str = strjoin(string(badIC),', ');
    fprintf('HEOG ICs: %s\n',str);
else
    fprintf('No Saccade ICs were identified using template correlation (R>%1.2f).\n',saccadeTreshold);
end

% Log
EEG.ALSUTRECHT.ica.saccadeTemplateCorr.pval = corrMat;
EEG.ALSUTRECHT.ica.saccadeTemplateCorr.bics = badIC(:);
EEG.ALSUTRECHT.ica.saccadeTemplateCorr.cvec = ones(length(EEG.ALSUTRECHT.ica.saccadeTemplateCorr.bics),1);
EEG.ALSUTRECHT.ica.saccadeTemplateCorr.clss = 'Saccade';

%% ========================================================================
% Final log

% *Blink ICs
blink1 = EEG.ALSUTRECHT.ica.corr.bics(EEG.ALSUTRECHT.ica.corr.cvec==2);
blink2 = EEG.ALSUTRECHT.ica.icablinkmetrics.bics;
blink3 = EEG.ALSUTRECHT.ica.blinkTemplateCorr.bics;
blinkAll = unique([blink1(:); blink2(:); blink3(:)]);
ICsMostLikelyBlink = false(NICA,1);
ICsMostLikelyBlink(blinkAll) = true;

% *Saccades ICs
saccade1 = EEG.ALSUTRECHT.ica.corr.bics(EEG.ALSUTRECHT.ica.corr.cvec==3);
saccade2 = EEG.ALSUTRECHT.ica.saccadeTemplateCorr.bics;
saccadeAll = unique([saccade1(:); saccade2(:)]);
ICsMostLikelySaccade = false(NICA,1);
ICsMostLikelySaccade(saccadeAll) = true;

% *Eye ICs
eye1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics)==3);
ICsMostLikelyEyeICLabel = false(NICA,1);
ICsMostLikelyEyeICLabel(eye1) = true;

ICsMostLikelyEye = ICsMostLikelyBlink | ICsMostLikelySaccade | ICsMostLikelyEyeICLabel;

% *Muscle ICs
muscle1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics)==2);
muscle2 = EEG.ALSUTRECHT.ica.emg.bics;
muscleAll = unique([muscle1(:); muscle2(:)]);
ICsMostLikelyMuscle = false(NICA,1);
ICsMostLikelyMuscle(muscleAll) = true;

% *Complex ICs
ICsMostLikelyComplex = ICsMostLikelyMuscle & ICsMostLikelyEye;

% Update
ICsMostLikelyMuscle(ICsMostLikelyComplex)     = false;
ICsMostLikelyBlink(ICsMostLikelyComplex)      = false;
ICsMostLikelySaccade(ICsMostLikelyComplex)    = false;
ICsMostLikelyEyeICLabel(ICsMostLikelyComplex) = false;
ICsMostLikelyEye = ICsMostLikelyBlink | ICsMostLikelySaccade | ICsMostLikelyEyeICLabel;

% *Channel ICs
% channel1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics)==6);
% channel2 = EEG.ALSUTRECHT.ica.spatialSmoothness.bics;
% channelAll = unique([channel1(:); channel2(:)]);
ICsMostLikelyChannel = false(NICA,1);
% ICsMostLikelyChannel(channelAll) = true;

% % These are likely muscle?
% tmp = zscore(EEG.icawinv(:,ICsMostLikelyMuscle));
% % figure; imagesc(tmp);
% ICsMostLikelyChannelWrong = sum(tmp>3,1) > 0 & sum(tmp<-3,1) >0;

% % Sometimes muscle ICs are marked as channel ICs
% ICsMostLikelyChannelWrong = ICsMostLikelyChannel & ICsMostLikelyMuscle;
% if any(ICsMostLikelyChannelWrong)
%     ICsMostLikelyChannel(ICsMostLikelyChannelWrong) = false;
% end

% % Maybe do not bother with "higher" bad channel ICs, low variance anyway
% ICsMostLikelyChannel(35:end) = false;

% *Heart ICs
heart1 = EEG.ALSUTRECHT.ica.ICLabel.bics(EEG.ALSUTRECHT.ica.ICLabel.cvec(EEG.ALSUTRECHT.ica.ICLabel.bics)==4);
heart2 = EEG.ALSUTRECHT.ica.corr.bics(EEG.ALSUTRECHT.ica.corr.cvec==1);
heart3 = EEG.ALSUTRECHT.ica.ctps.bics;
heartAll = unique([heart1(:); heart2(:); heart3(:)]);
ICsMostLikelyHeart = false(NICA,1);
ICsMostLikelyHeart(heartAll) = true;

% Sometimes muscle or channel ICs are marked as heart ICs
ICsMostLikelyHeartWrong = ICsMostLikelyHeart & (ICsMostLikelyMuscle |  ICsMostLikelyComplex | ICsMostLikelyChannel);
if any(ICsMostLikelyHeartWrong)
    ICsMostLikelyHeart(ICsMostLikelyHeartWrong) = false;
end

% % *wICA ICs
% ICsforwICA = false(NICA,1);
% ICsforwICA(ICsMostLikelyEye | ICsMostLikelyHeart | ICsMostLikelyChannel) = true;
% ICsforwICA(ICsMostLikelyMuscle | ICsMostLikelyComplex) = false;

% There should be no overlap
assert(max(sum([ICsMostLikelyEye, ICsMostLikelyMuscle, ICsMostLikelyComplex],2))==1);

% Log
EEG.ALSUTRECHT.ica.ICsMostLikelyBlink   = ICsMostLikelyBlink;
EEG.ALSUTRECHT.ica.ICsMostLikelySaccade = ICsMostLikelySaccade;
EEG.ALSUTRECHT.ica.ICsMostLikelyEye     = ICsMostLikelyEye;
EEG.ALSUTRECHT.ica.ICsMostLikelyMuscle  = ICsMostLikelyMuscle;
EEG.ALSUTRECHT.ica.ICsMostLikelyComplex = ICsMostLikelyComplex;
EEG.ALSUTRECHT.ica.ICsMostLikelyHeart   = ICsMostLikelyHeart;
EEG.ALSUTRECHT.ica.ICsMostLikelyChannel = ICsMostLikelyChannel;
% EEG.ALSUTRECHT.ica.ICsforwICA           = ICsforwICA;

% Labels:
% 1 'Brain'
% 2 'Muscle'
% 3 'Eye'
% 4 'Heart'
% 5 'Line Noise'
% 6 'Channel Noise'
% 7 'Other'

EEG.ALSUTRECHT.ica.combi.report                       = EEG.ALSUTRECHT.ica.ICLabel.cvec;
EEG.ALSUTRECHT.ica.combi.report(ICsMostLikelyMuscle)  = 2;
EEG.ALSUTRECHT.ica.combi.report(ICsMostLikelyEye)     = 3;
EEG.ALSUTRECHT.ica.combi.report(ICsMostLikelyHeart)   = 4;
EEG.ALSUTRECHT.ica.combi.report(ICsMostLikelyChannel) = 6;

% assert(length(EEG.ALSUTRECHT.ica.combi.bics)==length(EEG.ALSUTRECHT.ica.combi.prbs));
% assert(length(EEG.ALSUTRECHT.ica.combi.bics)==length(EEG.ALSUTRECHT.ica.combi.lbls));

end