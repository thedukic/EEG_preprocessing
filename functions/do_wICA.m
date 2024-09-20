function EEG = do_wICA(EEG,EXT,cfg)
%
% Function is still in the development
% Based on: RELAX_wICA_on_ICLabel_artifacts
%
% Input:
%   K = Threshold multiplier for wavelet thresholding.
%       Higher thresh -> Less strict
%   L = Level set for stationary wavelet transform.
%       Higher levels give better freq resolution, but less temp resolution
%   W = Wavelet family to use.
%       Type "wavenames" to see a list of possible wavelets
%
% More info:
% https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2018.00097/full
% Given that the magnitude of artifacts can be far greater than
% that of neurophysiological signals, the component time series
% whose amplitudes are large enough to survive the wavelet-thresholding
% are taken as the artifact timeseries.
%
% Treshold of 0   -> whole IC rejected
% Higher treshold ->
% 1. Less of the IC is considered as noise
% 2. Focuses on low amplitudes of the IC,
%    so it is good for eye and channel artifacts
%
% SDukic, July 2024
% =========================================================================

adjustEye = false;

% Make sure IC activations are present
if isempty(EEG.icaact)
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end

NICA = size(EEG.icaact,1);
dataICs = reshape(double(EEG.icaact),NICA,[]);

%% NEW CODE

% Labels:
% 1 'Brain'
% 2 'Muscle'
% 3 'Eye'
% 4 'Heart'
% 5 'Line Noise'
% 6 'Channel Noise'
% 7 'Other'
% ICsMostLikelyNotBrain = (I==2 | I==3)'; % Muscle/Eye
% ICsMostLikelyNotBrain = (I>=2 & I<=6)';
% ICsMostLikelyNotBrain = (I>1)';
% ICsMostLikelyEye = (I==3)';

ICsMostLikelyNotBrain = false(NICA,1);
ICsMostLikelyBlink    = false(NICA,1);

% All bad ICs
ICsMostLikelyNotBrain(unique(EEG.ALSUTRECHT.ica.combi.bics)) = true;
% Only blinks ICs
ICsMostLikelyBlink(EEG.ALSUTRECHT.ica.extra3.bics) = true;
% Only EMG ICs
ICsMostLikelyMuscle = EEG.ALSUTRECHT.ica.extra2.ICsMostLikelyMuscle;

% Remvoe these as they are likely very complex ICs
% assert(~any(ICsMostLikelyMuscle&ICsMostLikelyBlink));
ICsMostLikelyComplex = ICsMostLikelyMuscle&ICsMostLikelyBlink;

% Report
fprintf('\nBlink ICs, N = %d.\n',sum(ICsMostLikelyBlink));
fprintf('Muscle ICs, N = %d.\n',sum(ICsMostLikelyMuscle));
fprintf('Complex ICs, N = %d (will be removed completely).\n',sum(ICsMostLikelyComplex));

% %%
% options.muscleFreqIn    = [7, 70];
% options.Freq_to_compute = [1, 100];
%
% % Calculate pwelch to enable detection of log-freq log-power slopes, indicative of muscle activity
% % Resize EEG.icaact if required
% if size(EEG.icaact,3) > 0
%     eegData = reshape(EEG.icaact,size(EEG.icaact,1),[]);
% else
%     eegData = EEG.icaact;
% end
%
% [pxx,fp] = pwelch(eegData',size(eegData,2),[],size(eegData,2),EEG.srate);
% FFTout = pxx';
% fp = fp';
%
% % Calculate FFT bins
% freq = options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
% fftBins = zeros(size(FFTout,1),size(freq,2)); % preallocate
% for a=1:size(freq,2)
%     [~, index1]=min(abs(fp-((freq(1,a)-0.25))));
%     [~, index2]=min(abs(fp-((freq(1,a)+0.25))));
%     fftBins(:,a) = mean(FFTout(:,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
% end
%
% %% better muscle comp_number identification:
% comps = size(EEG.icaact,1);
% options.muscleFreqEx = [50-2 50+2];
%
% muscleRatio = NaN(1,NICA);
% for compNum = 1:comps
%     % Define frequencies to include in the analysis
%     if ~isempty(options.muscleFreqIn)
%         [~,fin1] = min(abs(options.muscleFreqIn(1) - freq));
%         [~,fin2] = min(abs(options.muscleFreqIn(2) - freq));
%         freqHz = freq(1,fin1:fin2);
%         freqPow = fftBins(compNum,fin1:fin2);
%     else
%         freqHz = freq;
%         freqPow = fftBins(compNum,:);
%     end
%     % Define frequencies to exclude from fit
%     if ~isempty(options.muscleFreqEx)
%         [~,fex1] = min(abs(options.muscleFreqEx(1) - freqHz));
%         [~,fex2] = min(abs(options.muscleFreqEx(2) - freqHz));
%         freqHz(fex1:fex2) = [];
%         freqPow(fex1:fex2) = [];
%     end
%     % Fit linear regression to log-log data
%     p = polyfit(log(freqHz),log(freqPow),1);
%     % Store the slope
%     muscleRatio(compNum) = p(1);
% end
% ICsMostLikelyMuscle = muscleRatio>=cfg.bch.muscleSlopeThreshold;
% % ICsMostLikelyMuscle = (muscle_ICs==1);
%
% %% Use icablinkmetrics to double-check for blink components that ICLabel might have missed:
% if exist('icablinkmetrics','file') == 2
%     try
%         % icablinkmetricsout = icablinkmetrics(EEG0,'ArtifactChannel',EEG0.data(strcmp({EEG0.chanlocs.labels},'VEOG'),:),'Alpha',0.001,'VisualizeData','False');
%         icablinkmetricsout = icablinkmetrics(EEG,'ArtifactChannel',mean(EEG.data(ismember({EEG.chanlocs.labels},cfg.ica.blinkchans),:),1),'Alpha',0.001,'VisualizeData','False');
%         if any(icablinkmetricsout.identifiedcomponents>0)
%             fprintf('icablinkmetrics has identified %d eye component(s).\n',length(icablinkmetricsout.identifiedcomponents));
%             ICsMostLikelyNotBrain(icablinkmetricsout.identifiedcomponents) = true;
%             ICsMostLikelyEye(icablinkmetricsout.identifiedcomponents)      = true;
%         else
%             fprintf('icablinkmetrics has not identified any eye components.\n');
%         end
%     catch
%         warning('icablinkmetrics has failed...');
%     end
% end

%% Padding
NTPT = size(dataICs,2);
assert(NTPT==EEG.pnts);

check_padding_required = mod(NTPT,2^5);
if check_padding_required ~= 0
    padding = zeros(1, (2^5)-check_padding_required);
else
    padding = [];
end

%% Perform wavelet thresholding
fprintf('\nWavelet tresholding:\n');
artifact_comp = zeros(NICA,NTPT+length(padding));

cnt = 0;
for i = 1:NICA
    % Do this step on nonmuscle and noncomplex ICs only!
    if ICsMostLikelyNotBrain(i) && ~ICsMostLikelyMuscle(i) && ~ICsMostLikelyComplex(i)
        cnt = cnt+1;
        labels = EEG.ALSUTRECHT.ica.combi.lbls(EEG.ALSUTRECHT.ica.combi.bics==i);
        if length(labels)>1, labels = {strjoin(labels,'/')}; end
        fprintf('%d) IC%d: %s - Wavelet tresholding... ',cnt,i,labels{1});

        % Pad the component with zeros if required
        if ~isempty(padding)
            padded_comp = [dataICs(i,:), padding];
        else
            padded_comp = dataICs(i,:);
        end

        [wavelet_threshold,threshold_type,~] = ddencmp('den','wv',padded_comp); % automatically obtain wavelet enhancement threshold
        if ICsMostLikelyBlink(i)
            % Increase threshold for blink components based on optimal results in our informal testing
            wavelet_threshold = wavelet_threshold*1.5;
        else
            wavelet_threshold = wavelet_threshold*1;
        end
        wavelet_transform = swt(padded_comp,5,'coif5'); % apply stationary wavelet transform to each component to reduce neural contribution to component
        thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold); % remove negligible values by applying thresholding
        artifact_comp(i,:) = iswt(thresholded_wavelet_transform,'coif5'); % use inverse wavelet transform to obtained the wavelet transformed component

        fprintf('Used treshold %1.2f\n',wavelet_threshold);
        clearvars thresholded_wavelet_transform padded_comp wavelet_threshold threshold_type wavelet_transform
    end
end

%% Padding
% Pad non-artifact components with 0s in the same way that the artifact components were padded:
% if sum(ICsMostLikelyNotBrain)==0
%     artifact_comp(1,:) = zeros(1,size(EEG.data,2));
%     artifact_comp      = [artifact_comp(:,:), padding];
% end

% for i = 1:NICA
%     if ~ICsMostLikelyNotBrain(i)
%         artifact_comp(i,:) = zeros(1,size(artifact_comp,2));
%     end
% end

% Remove padding
if ~isempty(padding)
    artifact_comp = artifact_comp(:,1:end-numel(padding));
end

%% Obtain muscle artifact for subtraction by highpass filtering data instead of wICA:
fprintf('Adjusting the musle components (N = %d) by keeping only >15 Hz...\n',sum(ICsMostLikelyMuscle));

[z1, p1] = butter(2, 15./(EEG.srate/2), 'high');
if any(ICsMostLikelyMuscle)
    artifact_comp(ICsMostLikelyMuscle,:) = filtfilt(z1,p1,dataICs(ICsMostLikelyMuscle,:)')';
end

%% Restrict wICA cleaning of blink components to just blink periods
if adjustEye
    fprintf('Adjusting the eye components (N = %d) by targeting only the eye artifacts...\n',sum(ICsMostLikelyBlink));

    moving_mean_length     = round(200/(1000/EEG.srate));
    blink_length_threshold = round(100/(1000/EEG.srate));

    if any(ICsMostLikelyBlink)
        assert(size(EEG.data,2)==size(EXT.data,2));
        maskEyeBlinks = detect_eog(EXT,400,false); % 400 ms
    end

    [z1, p1] = butter(2, [0.5 25]./(EEG.srate/2), 'bandpass');
    maskBlinksICA = zeros(NICA,NTPT);

    for i = 1:NICA
        if ICsMostLikelyBlink(i)
            IC_filtered = filtfilt(z1,p1,dataICs(i,:)')';
            [blink_periods,~,~] = isoutlier(IC_filtered,'median',ThresholdFactor=2);
            ix_blinkstart = find(diff(blink_periods)==1)+1;  % indices where BlinkIndexMetric goes from 0 to 1
            ix_blinkend   = find(diff(blink_periods)==-1);   % indices where BlinkIndexMetric goes from 1 to 0

            % Use an IQR threshold method to detect and mark blinks
            % [EEG, ~] = RELAX_blinks_IQR_method(EEG, EEG, RELAX_cfg);
            % blink_periods(EEG.RELAX.eyeblinkmask==1)=1; +/ 400ms
            assert(length(blink_periods)==length(maskEyeBlinks));
            blink_periods(maskEyeBlinks) = 1;

            if ~isempty(ix_blinkstart)
                if ix_blinkend(1,1)<ix_blinkstart(1,1); ix_blinkend(:,1)=[]; end % if the first downshift occurs before the upshift, remove the first value in end
                if ix_blinkend(1,size(ix_blinkend,2))<ix_blinkstart(1,size(ix_blinkstart,2)); ix_blinkstart(:,size(ix_blinkstart,2))=[];end % if the last upshift occurs after the last downshift, remove the last value in start
                BlinkThresholdExceededLength=ix_blinkend-ix_blinkstart; % length of consecutive samples where blink threshold was exceeded
                BlinkRunIndex = find(BlinkThresholdExceededLength<round(blink_length_threshold/(1000/EEG.srate))); % find locations where blink threshold was not exceeded by more than X ms
                % find latency of the max voltage within each period where the blink
                % threshold was exceeded:
                if size(BlinkRunIndex,2)>0
                    for x=1:size(BlinkRunIndex,2)
                        o=ix_blinkstart(BlinkRunIndex(x));
                        c=ix_blinkend(BlinkRunIndex(x));
                        if c-o<round(blink_length_threshold/(1000/EEG.srate))
                            blink_periods(1,o:c)=0;
                        end
                    end
                end
            end
            padded_blink_periods=double(blink_periods);
            for c=flip(1:size(padded_blink_periods,2)-(moving_mean_length+1))
                if padded_blink_periods(1,c)==1
                    padded_blink_periods(1,c:c+moving_mean_length)=1;
                end
            end
            for c=(moving_mean_length+1):size(padded_blink_periods,2)
                if padded_blink_periods(1,c)==1
                    padded_blink_periods(1,c-moving_mean_length:c)=1;
                end
            end
            maskBlinksICA(i,:) = movmean(padded_blink_periods,[moving_mean_length moving_mean_length]);
            artifact_comp(i,:) = artifact_comp(i,:).*maskBlinksICA(i,:);
        end
    end
else
    fprintf('Non-targeted apporach is used for eye components.\n');
end

%% Remove these completely
artifact_comp(ICsMostLikelyComplex,:) = dataICs(ICsMostLikelyComplex,:);

%% Remove artifact and reconstruct data:
artifacts = EEG.icawinv*artifact_comp;
artifacts = reshape(artifacts,size(artifacts,1),NTPT,EEG.trials);

% chaneeg = strcmp({EEG.chanlocs.type},'EEG');
% EEGNEW = EEG;
% EEGNEW.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;
% vis_artifacts(EEGNEW,EEG);

chaneeg = strcmp({EEG.chanlocs.type},'EEG');
EEG.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;

%% Log
% % wICA
% eyeICs    = find(ICsMostLikelyEye(:));
% muscleICs = find(ICsMostLikelyMuscle(:));
% EEG.ALSUTRECHT.ica.wica.bics = [eyeICs; muscleICs];
%
% % Final
% EEG.ALSUTRECHT.ica.final.bics = [EEG.ALSUTRECHT.ica.combi.bics(:); EEG.ALSUTRECHT.ica.wica.bics(:)];
% EEG.ALSUTRECHT.ica.final.lbls = [EEG.ALSUTRECHT.ica.combi.lbls; repmat({'Eye2'},length(eyeICs),1); repmat({'Muscle2'},length(muscleICs),1)];
% assert(length(EEG.ALSUTRECHT.ica.final.bics)==length(EEG.ALSUTRECHT.ica.final.lbls));

NBIC = length(unique(EEG.ALSUTRECHT.ica.combi.bics));
labels = EEG.ALSUTRECHT.ica.combi.lbls;

EEG.ALSUTRECHT.ica.proportionArtifactICsReducedbywICA    = NBIC./NICA;
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICABlink   = sum(ICsMostLikelyBlink);
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAmuscle  = sum(ICsMostLikelyMuscle);
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAheart   = sum(ismember(labels,{'ECG','Heart'}));
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAchannel = sum(ismember(labels,'Channel Noise'));

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA cleaning\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of bad ICs:    %d\n',NBIC);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of blink IC:   %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICABlink);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of muscle IC:  %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAmuscle);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of heart IC:   %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAheart);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of channel IC: %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAchannel);

end