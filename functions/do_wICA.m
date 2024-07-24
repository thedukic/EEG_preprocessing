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

% Make sure IC activations are present
if isempty(EEG.icaact)
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end
NICA = size(EEG.icaact,1);
dataICs = reshape(EEG.icaact, size(EEG.icaact,1), []);

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
check_padding_required = mod(size(dataICs,2),2^5);
if check_padding_required ~=0
    padding = zeros(1,(2^5)-check_padding_required);
else
    padding = [];
end

%% Perform wavelet thresholding on eye movements (and also other components if selected), identified by ICLabel:
disp('Using targeted approach to clean artifacts:');
cnt = 0;
for i = 1:NICA
    if ICsMostLikelyNotBrain(i)
        cnt = cnt+1;
        labels = EEG.ALSUTRECHT.ica.combi.lbls(EEG.ALSUTRECHT.ica.combi.bics==i);
        if length(labels)>1, labels = {strjoin(labels,'/')}; end
        fprintf('%d) IC%d: %s - Wavelet tresholding... ',cnt,i,labels{1});

        if ~isempty(padding)
            padded_comp = [dataICs(i,:),padding]; % pad the component with zeros if required
        else
            padded_comp = dataICs(i,:);
        end

        [wavelet_threshold,threshold_type,~] = ddencmp('den','wv',padded_comp); % automatically obtain wavelet enhancement threshold
        if ICsMostLikelyBlink(i)
            wavelet_threshold = wavelet_threshold*2; % increase threshold for blink components based on optimal results in our informal testing
        else
            wavelet_threshold = wavelet_threshold*1;
        end
        wavelet_transform = swt(padded_comp,5,'coif5'); % apply stationary wavelet transform to each component to reduce neural contribution to component
        thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold); % remove negligible values by applying thresholding
        artifact_comp(i,:)  = iswt(thresholded_wavelet_transform,'coif5'); % use inverse wavelet transform to obtained the wavelet transformed component

        fprintf('Used treshold %1.2f\n',wavelet_threshold);
        clearvars thresholded_wavelet_transform padded_comp wavelet_threshold threshold_type wavelet_transform
    end
end

%% Padding
% Pad non-artifact components with 0s in the same way that the artifact components were padded:
if sum(ICsMostLikelyNotBrain)==0
    artifact_comp(1,:) = zeros(1,size(EEG.data,2));
    artifact_comp      = [artifact_comp(:,:),padding];
end

for i = 1:NICA
    if ~ICsMostLikelyNotBrain(i)
        artifact_comp(i,:) = zeros(1,size(artifact_comp,2));
    end
end

% Remove padding
if ~isempty(padding)
    artifact_comp = artifact_comp(:,1:end-numel(padding));
end

%% Restrict wICA cleaning of blink components to just blink periods
fprintf('Adjusting the eye components (N = %d)...\n',sum(ICsMostLikelyBlink));

moving_mean_length     = round(200/(1000/EEG.srate));
blink_length_threshold = round(100/(1000/EEG.srate));
clearvars M

if any(ICsMostLikelyBlink)
    assert(size(EEG.data,2)==size(EXT.data,2));
    eyeBlinksMask = detect_eog(EXT,400); % 400 ms
end

for i = 1:NICA
    if ICsMostLikelyBlink(i)
        [z1, p1] = butter(2, [0.5 25]./(EEG.srate/2), 'bandpass');
        dataIn = dataICs(i,:)';
        dataFilt1 = filtfilt(z1,p1,double(dataIn));
        IC_filtered = dataFilt1';
        [blink_periods,~,~] = isoutlier(IC_filtered,'median',ThresholdFactor=2);
        ix_blinkstart = find(diff(blink_periods)==1)+1;  % indices where BlinkIndexMetric goes from 0 to 1
        ix_blinkend   = find(diff(blink_periods)==-1);   % indices where BlinkIndexMetric goes from 1 to 0

        % [EEG, ~] = RELAX_blinks_IQR_method(EEG, EEG, RELAX_cfg); % use an IQR threshold method to detect and mark blinks
        % blink_periods(EEG.RELAX.eyeblinkmask==1)=1; +/ 400ms
        assert(length(blink_periods)==length(eyeBlinksMask));
        blink_periods(eyeBlinksMask) = 1;

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
        M(i,:) = movmean(padded_blink_periods,[moving_mean_length moving_mean_length]);
        artifact_comp(i,:) = (artifact_comp(i,:).*M(i,:));
    end
end

%% Obtain muscle artifact for subtraction by highpass filtering data instead of wICA:
fprintf('Adjsuting the musle components (N = %d)...\n',sum(ICsMostLikelyMuscle));

for i = 1:NICA
    if ICsMostLikelyMuscle(i)
        [z1, p1] = butter(2, 15./(EEG.srate/2), 'high');
        dataIn = dataICs(i,:)';
        dataFilt1 = filtfilt(z1,p1,double(dataIn));
        artifact_comp(i,:) = dataFilt1';
    end
end

%% Remove artifact and reconstruct data:
artifacts = EEG.icawinv*artifact_comp;
artifacts = reshape(artifacts,size(artifacts,1),EEG.pnts,EEG.trials);

chaneeg = strcmp({EEG.chanlocs.type},'EEG');
EEG.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;

%% OLD CODE
% if ~isempty(ICsArtifact)
%     % Wavelet thresholding only on (some) artifact components identified by ICLabel
%     % Doesnt seems to work on EMG artifacts
%     NBIC = length(ICsArtifact);
%     fprintf('\nPerforming wavelet thresholding on %d bad ICs detected by ICLabel...\n',NBIC);
%
%     % % Select EOG
%     % chanveog = strcmp({EEG.chanlocs.labels},'VEOG');
%     % chanheog = strcmp({EEG.chanlocs.labels},'HEOG');
%     % dataveog = EEG.data(chanveog,:);
%     % dataheog = EEG.data(chanheog,:);
%     % % figure; multisignalplot([dataveog; dataheog],EEG.srate,'r');
%     %
%     % [bl, al] = butter(4,20/(EEG.srate/2),'low');
%     % assert(isstable(bl,al));
%     % dataveog = filtfilt(bl,al,dataveog);
%     % dataheog = filtfilt(bl,al,dataheog);
%
%     % Wavelet padding, 2^level
%     modulus = mod(size(IC,2),2^L);
%     if modulus~=0
%         extra = zeros(1,(2^L)-modulus);
%         IC = [IC, repmat(extra,NICA,1)];
%         % dataveog = [dataveog, extra];
%         % dataheog = [dataheog, extra];
%     else
%         extra = [];
%     end
%
%
%     % wIC = zeros(size(IC));
%     % for i = 1:NBIC
%     %     label = EEG.ALSUTRECHT.ica.ICLabel_clss{EEG.ALSUTRECHT.ica.ICLabel_cvec(ICsArtifact(i))};
%     %
%     %     if strcmpi(label,'Muscle') || strcmpi(label,'Heart')
%     %         % Muscle >20 Hz
%     %         % Heart  0-150 Hz
%     %         % -> Not suitable for wavelet tresholding
%     %         % Maybe useful EMD?
%     %         fprintf('%d. %s - Removing completely...\n',i,label);
%     %         wIC(ICsArtifact(i),:) = IC(ICsArtifact(i),:);
%     %
%     %         % [imf,residual,info] = emd(X,'Interpolation','pchip','MaxNumIMF',5)
%     %
%     %     elseif strcmpi(label,'Channel Noise') || strcmpi(label,'Eye')
%     %         fprintf('%d. %s - Wavelet tresholding... ',i,label);
%     %
%     %         [thresh,sorh,~] = ddencmp('den','wv',IC(ICsArtifact(i),:)); % get automatic threshold value
%     %         thresh = thresh*K;                                          % multiply threshold by scalar
%     %         swc    = swt(IC(ICsArtifact(i),:),L,W);                     % use stationary wavelet transform (SWT) to wavelet transform the ICs
%     %         Y      = wthresh(swc,sorh,thresh);                          % threshold the wavelet to remove small values
%     %         wIC(ICsArtifact(i),:) = iswt(Y,W);                          % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
%     %
%     %         fprintf('Used treshold %1.2f\n',thresh);
%     %
%     %         % elseif strcmpi(label,'Eye')
%     %         %     fprintf('%d. %s - Wavelet tresholding... ',i,label);
%     %         %
%     %         %     WLT ='bior3.1';
%     %         %     [thresh,sorh,~] = ddencmp('den','wv',IC(ICsArtifact(i),:)); % get automatic threshold value
%     %         %     swc = swt(IC(ICsArtifact(i),:),LVL,WLT);                    % use stationary wavelet transform (SWT) to wavelet transform the ICs
%     %         %     Y = wthresh(swc,sorh,thresh);                               % threshold the wavelet to remove small values
%     %         %     wIC(ICsArtifact(i),:) = iswt(Y,WLT);                        % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
%     %         %
%     %         %     fprintf('Used treshold %1.2f\n',thresh);
%     %         %
%     %         %     % VEOG and HEOG movements
%     %         %     fprintf('%d. %s - Wavelet tresholding...\n',i,label);
%     %         %     [rho0v, p0v] = corr(IC(ICsArtifact(i),:)',dataveog);
%     %         %     [rho0h, p0h] = corr(IC(ICsArtifact(i),:)',dataheog);
%     %         %
%     %         %     % What if this fails?
%     %         %     if abs(rho0v)>=abs(rho0h)
%     %         %         fprintf('This eye IC is similar to VEOG: r = %1.2f, p = %1.2f\n',rho0v,p0v);
%     %         %         dataeog = dataveog;
%     %         %     else
%     %         %         fprintf('This eye IC is similar to HEOG: r = %1.2f, p = %1.2f\n',rho0h,p0h);
%     %         %         dataeog = dataheog;
%     %         %     end
%     %         %
%     %         %     [thresh,sorh,~] = ddencmp('den','wv',IC(ICsArtifact(i),:));
%     %         %     swc = swt(IC(ICsArtifact(i),:),LVL,WLT);
%     %         %
%     %         %     r = NaN(NK,1);
%     %         %     p = NaN(NK,1);
%     %         %     for j = 1:NK
%     %         %         thresh2 = thresh*K(j);
%     %         %         Y = wthresh(swc,sorh,thresh2);
%     %         %         tmp = iswt(Y,WLT);
%     %         %         % tmp = tmp(1:end-numel(extra));
%     %         %         [r(j), p(j)] = corr(tmp',dataeog);
%     %         %     end
%     %         %     % figure; plot(K,abs(r));
%     %         %     [a,b] = max(abs(r));
%     %         %     thresh = thresh*K(b);
%     %         %     fprintf('Max{r} : r = %1.2f, p = %1.2f\n',a,p(b));
%     %         %     fprintf('Final treshold %1.2f; used factor %1.2f\n',thresh,K(b));
%     %         %
%     %         %     Y = wthresh(swc,sorh,thresh);
%     %         %     wIC(ICsArtifact(i),:) = iswt(Y,WLT);
%     %         %     % clearvars y thresh sorh swc
%     %     end
%     %
%     %     % figure; multisignalplot(swc,EEG.srate,'r');
%     %     % figure; multisignalplot(Y,EEG.srate,'r');
%     %     % % figure; multisignalplot([sig; wIC(ICsArtifact(i),:)],EEG.srate,'r');
%     %     % tmp = [];
%     %     % tmp.data  = wIC(ICsArtifact(i),:);
%     %     % tmp.srate = EEG.srate;
%     %     % [pow, freq] = checkpowerspectrum(tmp,1,[]);
%     % end
%
%     % Remove extra padding
%     if ~isempty(extra)
%         wIC = wIC(:,1:end-numel(extra));
%     end
%
%     % % Plot the ICs vs. wICs
%     % figure;
%     % subplot(3,1,1);
%     % multisignalplot(IC,EEG.srate,'r');
%     % title('ICs');
%     % subplot(3,1,2);
%     % multisignalplot(wIC,EEG.srate,'r');
%     % title('wICs')
%     % subplot(3,1,3);
%     % multisignalplot(IC-wIC,EEG.srate,'r');
%     % title('Difference (IC - wIC)');
%
%     % Reconstruct artifacts and subtract
%     artifacts = EEG.icawinv*wIC;
%     artifacts = reshape(artifacts,size(artifacts,1),EEG.pnts,EEG.trials);
%
%     chaneeg   = strcmp({EEG.chanlocs.type},'EEG');
%     EEG.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;
%
%     % % Visualise
%     % % vis_artifacts(EEGNEW,EEG);
%     %
%     % NOISE = EEG;
%     % NOISE.data(chaneeg,:,:) = reshape(artifacts,size(artifacts,1),EEG.pnts,EEG.trials);
%     % vis_artifacts(EEG,NOISE);
%     % % [pow, freq] = checkpowerspectrum(NOISE,1:5,[]);
% else
%     fprintf('Skipping wavelet thresholding since there ar no obvious bad ICs detected...\n');
%     NBIC = 0;
% end

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