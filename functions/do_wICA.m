function EEG = do_wICA(EEG,EXT,flagEyeIC)
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
% SDukic, October 2024
% =========================================================================

% Double-check
EEG = eeg_checkset(EEG,'ica');

% Make sure IC activations are present
EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);

%% Extract bad ICs
ICsMostLikelyBlink   = EEG.ALSUTRECHT.ica.ICsMostLikelyBlink;
ICsMostLikelyMuscle  = EEG.ALSUTRECHT.ica.ICsMostLikelyMuscle;
ICsMostLikelyComplex = EEG.ALSUTRECHT.ica.ICsMostLikelyComplex;
ICsforwICA           = EEG.ALSUTRECHT.ica.ICsforwICA;

% Report
fprintf('\nwICA ICs,    N = %d.\n',sum(ICsforwICA));
fprintf('Blink ICs,   N = %d (targeted-wavelet tresholding).\n',sum(ICsMostLikelyBlink));
fprintf('Muscle ICs,  N = %d (simple >15 Hz filtering).\n',sum(ICsMostLikelyMuscle));
fprintf('Complex ICs, N = %d (removed completely).\n',sum(ICsMostLikelyComplex));
% fprintf('Other ICs (HEOG, ECG, channel noise), N = %d (wavelet tresholding).\n',sum(~)); % change based on flagEyeIC

%% Padding
dataICs = double(EEG.icaact);
assert(ismatrix(dataICs));

[NICA, NTPT] = size(dataICs);
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
    % Do this step on nonmuscle and noncomplex ICs only.
    % 1. EMG ICs are filtered out >15Hz
    % 2. Complex ICs are completely removed

    % -> Eye / Heart / Channel
    if ICsforwICA(i)
        cnt = cnt+1;
        thisLabel = EEG.ALSUTRECHT.ica.ICLabel.clss{EEG.ALSUTRECHT.ica.combi.report(i)};
        fprintf('%d) IC%d: %s... ',cnt,i,thisLabel);

        % Pad the component with zeros if required
        if ~isempty(padding)
            padded_comp = [dataICs(i,:), padding];
        else
            padded_comp = dataICs(i,:);
        end

        [wavelet_threshold,threshold_type,~] = ddencmp('den','wv',padded_comp); % automatically obtain wavelet enhancement threshold

        % % Increase threshold for blink components based on optimal results in our informal testing
        % if ICsMostLikelyBlink(i)
        %     wavelet_threshold = 2*wavelet_threshold;
        % else
        %     wavelet_threshold = wavelet_threshold;
        % end
        wavelet_threshold = 0.8*wavelet_threshold;

        wavelet_transform = swt(padded_comp,5,'coif5'); % apply stationary wavelet transform to each component to reduce neural contribution to component
        thresholded_wavelet_transform = wthresh(wavelet_transform,threshold_type,wavelet_threshold); % remove negligible values by applying thresholding
        artifact_comp(i,:) = iswt(thresholded_wavelet_transform,'coif5'); % use inverse wavelet transform to obtained the wavelet transformed component

        fprintf('Used treshold %1.2f\n',wavelet_threshold);
        clearvars thresholded_wavelet_transform padded_comp wavelet_threshold threshold_type wavelet_transform
    end
end

%% Remove padding
if ~isempty(padding)
    artifact_comp = artifact_comp(:,1:end-numel(padding));
end

%% Obtain muscle artifact for subtraction by highpass filtering data instead of wICA:
if any(ICsMostLikelyMuscle)
    fprintf('\nAdjusting the musle components (N = %d) by keeping only >15 Hz...\n',sum(ICsMostLikelyMuscle));
    [bh, ah] = butter(2, 15/(EEG.srate/2), 'high');
    % artifact_comp(ICsMostLikelyMuscle,:) = filtfilt(bh,ah,dataICs(ICsMostLikelyMuscle,:)')';
    artifact_comp(ICsMostLikelyMuscle,:) = do_filteringcore(bh,ah,dataICs(ICsMostLikelyMuscle,:),EEG.event,EEG.srate);

else
    fprintf('No musle components were found...\n');
end

%% Restrict wICA cleaning of blink components to just blink periods
if any(ICsMostLikelyBlink)
    switch flagEyeIC
        case 'remove'
            fprintf('\nRemoving the eye components (N = %d) completely...\n',sum(ICsMostLikelyBlink));
            artifact_comp(ICsMostLikelyBlink,:) = dataICs(ICsMostLikelyBlink,:);

        case 'target'
            fprintf('\nAdjusting the eye components (N = %d) by targeting only the eye artifacts...\n',sum(ICsMostLikelyBlink));

            moving_mean_length     = round(200/(1000/EEG.srate));
            blink_length_threshold = round(100/(1000/EEG.srate));

            assert(size(EEG.data,2)==size(EXT.data,2));
            maskEyeBlinks = detect_veog(EXT,400,false); % 400 ms

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

                    % figure; hold on;
                    % plot(EXT.data(2,1:20*256));
                    % plot(200*maskBlinksICA(i,1:20*256));
                    % plot(EXT.data(3,1:20*256));
                end
            end

        case 'nontarget'
            fprintf('\nNon-targeted apporach is used for eye components.\n');
    end
end

%% Remove complex ones completely
if any(ICsMostLikelyComplex)
    fprintf('\nRemoving the complex components (N = %d) completely...\n',sum(ICsMostLikelyComplex));
    artifact_comp(ICsMostLikelyComplex,:) = dataICs(ICsMostLikelyComplex,:);
else
    fprintf('\nNo complex components were found.\n');
end

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
NBIC = sum(ICsMostLikelyBlink | ICsMostLikelyMuscle | ICsMostLikelyComplex | ICsforwICA);

EEG.ALSUTRECHT.ica.numberArtifactICs                   = NBIC;
EEG.ALSUTRECHT.ica.proportionArtifactICs               = NBIC./NICA;
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICA      = sum(ICsforwICA);
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICABlink = sum(ICsMostLikelyBlink);
EEG.ALSUTRECHT.ica.numberArtifactICsMuscle             = sum(ICsMostLikelyMuscle);
EEG.ALSUTRECHT.ica.numberArtifactICsComplex            = sum(ICsMostLikelyComplex);


fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA cleaning\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of bad ICs:    %d\n',EEG.ALSUTRECHT.ica.numberArtifactICs);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of total wIC   %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICA);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of blink wIC:  %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICABlink);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of muscle IC:  %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsMuscle);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of complex IC: %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsComplex);

% Not needed
EEG.icaact = [];

end