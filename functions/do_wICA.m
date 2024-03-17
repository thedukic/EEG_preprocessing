function EEG = do_wICA(EEG)
%
% See: RELAX_wICA_on_ICLabel_artifacts
%   K = threshold multiplier...multiplies the computed threshold from
%       "ddencmp" by this number. Higher thresh multipliers = less
%       "background" (or low amp. signal) is kept in the wICs.
%   L = level set for stationary wavelet transform. Higher levels give
%       better frequency resolution, but less temporal resolution.
%       Default = 5
%   wavename = wavelet family to use. type "wavenames" to see a list of
%       possible wavelets. (default = "coif5");
%
% https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2018.00097/full
% Given that the magnitude of artifacts can be far greater than
% that of neurophysiological signals, the component time series
% whose amplitudes are large enough to survive the wavelet-thresholding
% are taken as the artifact timeseries.
%
% Treshold of 0 => whole IC rejected
% Higher treshold =>
% 1. less of the IC is considered as noise
% 2. focuses on low amplitudes of the IC,
%    so it is good for eye and channel artifacts
%
% SDukic, March 2024
%

% Make sure IC activations are present
if isempty(EEG.icaact)
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end
NICA = size(EEG.icaact,1);
IC   = reshape(EEG.icaact,NICA,[]);

% labels = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};
ICsArtifact = find(EEG.reject.gcompreject);
% ICsArtifact = ICsArtifact(~strcmpi(EEG.ALSUTRECHT.ica.ICLabel_clss(EEG.ALSUTRECHT.ica.ICLabel_cvec(ICsArtifact)),'Muscle'));
% ICsArtifact = ICsArtifact(strcmpi(EEG.ALSUTRECHT.ica.ICLabel_clss(EEG.ALSUTRECHT.ica.ICLabel_cvec(ICsArtifact)),'Eye'));

if ~isempty(ICsArtifact)
    % Wavelet thresholding only on (some) artifact components identified by ICLabel
    % Doesnt seems to work on EMG artifacts
    NBIC = length(ICsArtifact);
    fprintf('\nPerforming wavelet thresholding on %d bad ICs detected by ICLabel...\n',NBIC);

    % % Select EOG
    % chanveog = strcmp({EEG.chanlocs.labels},'VEOG');
    % chanheog = strcmp({EEG.chanlocs.labels},'HEOG');
    % dataveog = EEG.data(chanveog,:);
    % dataheog = EEG.data(chanheog,:);
    % % figure; multisignalplot([dataveog; dataheog],EEG.srate,'r');
    %
    % [bl, al] = butter(4,20/(EEG.srate/2),'low');
    % assert(isstable(bl,al));
    % dataveog = filtfilt(bl,al,dataveog);
    % dataheog = filtfilt(bl,al,dataheog);

    % Wavelet padding, 2^level
    LVL = 5;
    modulus = mod(size(IC,2),2^LVL);
    if modulus~=0
        extra = zeros(1,(2^LVL)-modulus);
        IC = [IC, repmat(extra,NICA,1)];
        % dataveog = [dataveog, extra];
        % dataheog = [dataheog, extra];
    else
        extra = [];
    end

    % dataveog = dataveog';
    % dataheog = dataheog';
    % K = 0.5:0.01:2;
    % NK = length(K);

    wIC = zeros(size(IC));
    for i = 1:NBIC
        label = EEG.ALSUTRECHT.ica.ICLabel_clss{EEG.ALSUTRECHT.ica.ICLabel_cvec(ICsArtifact(i))};

        if strcmpi(label,'Muscle') || strcmpi(label,'Heart')
            fprintf('%d. %s - Removing completely...\n',i,label);
            wIC(ICsArtifact(i),:) = IC(ICsArtifact(i),:);

        elseif strcmpi(label,'Channel Noise') || strcmpi(label,'Eye')
            fprintf('%d. %s - Wavelet tresholding... ',i,label);
            
            WLT ='coif5';
            [thresh,sorh,~] = ddencmp('den','wv',IC(ICsArtifact(i),:)); % get automatic threshold value
            swc = swt(IC(ICsArtifact(i),:),LVL,WLT);                    % use stationary wavelet transform (SWT) to wavelet transform the ICs
            Y = wthresh(swc,sorh,thresh);                               % threshold the wavelet to remove small values
            wIC(ICsArtifact(i),:) = iswt(Y,WLT);                        % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
            
            fprintf('Used treshold %1.2f\n',thresh);

            % elseif strcmpi(label,'Eye')
            %     fprintf('%d. %s - Wavelet tresholding... ',i,label);
            %
            %     WLT ='bior3.1';
            %     [thresh,sorh,~] = ddencmp('den','wv',IC(ICsArtifact(i),:)); % get automatic threshold value
            %     swc = swt(IC(ICsArtifact(i),:),LVL,WLT);                    % use stationary wavelet transform (SWT) to wavelet transform the ICs
            %     Y = wthresh(swc,sorh,thresh);                               % threshold the wavelet to remove small values
            %     wIC(ICsArtifact(i),:) = iswt(Y,WLT);                        % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
            %
            %     fprintf('Used treshold %1.2f\n',thresh);
            %
            %     % VEOG and HEOG movements
            %     fprintf('%d. %s - Wavelet tresholding...\n',i,label);
            %     [rho0v, p0v] = corr(IC(ICsArtifact(i),:)',dataveog);
            %     [rho0h, p0h] = corr(IC(ICsArtifact(i),:)',dataheog);
            %
            %     % What if this fails?
            %     if abs(rho0v)>=abs(rho0h)
            %         fprintf('This eye IC is similar to VEOG: r = %1.2f, p = %1.2f\n',rho0v,p0v);
            %         dataeog = dataveog;
            %     else
            %         fprintf('This eye IC is similar to HEOG: r = %1.2f, p = %1.2f\n',rho0h,p0h);
            %         dataeog = dataheog;
            %     end
            %
            %     [thresh,sorh,~] = ddencmp('den','wv',IC(ICsArtifact(i),:));
            %     swc = swt(IC(ICsArtifact(i),:),LVL,WLT);
            %
            %     r = NaN(NK,1);
            %     p = NaN(NK,1);
            %     for j = 1:NK
            %         thresh2 = thresh*K(j);
            %         Y = wthresh(swc,sorh,thresh2);
            %         tmp = iswt(Y,WLT);
            %         % tmp = tmp(1:end-numel(extra));
            %         [r(j), p(j)] = corr(tmp',dataeog);
            %     end
            %     % figure; plot(K,abs(r));
            %     [a,b] = max(abs(r));
            %     thresh = thresh*K(b);
            %     fprintf('Max{r} : r = %1.2f, p = %1.2f\n',a,p(b));
            %     fprintf('Final treshold %1.2f; used factor %1.2f\n',thresh,K(b));
            %
            %     Y = wthresh(swc,sorh,thresh);
            %     wIC(ICsArtifact(i),:) = iswt(Y,WLT);
            %     % clearvars y thresh sorh swc
        end

        % figure; multisignalplot(swc,EEG.srate,'r');
        % figure; multisignalplot(Y,EEG.srate,'r');
        % % figure; multisignalplot([sig; wIC(ICsArtifact(i),:)],EEG.srate,'r');
        tmp = [];
        tmp.data  = wIC(ICsArtifact(i),:);
        tmp.srate = EEG.srate;
        [pow, freq] = checkpowerspectrum(tmp,1,[]);
    end

    % Remove extra padding
    if ~isempty(extra)
        wIC = wIC(:,1:end-numel(extra));
    end

    % % Plot the ICs vs. wICs
    % figure;
    % subplot(3,1,1);
    % multisignalplot(IC,EEG.srate,'r');
    % title('ICs');
    % subplot(3,1,2);
    % multisignalplot(wIC,EEG.srate,'r');
    % title('wICs')
    % subplot(3,1,3);
    % multisignalplot(IC-wIC,EEG.srate,'r');
    % title('Difference (IC - wIC)');

    % Reconstruct artifacts and subtract
    artifacts = EEG.icawinv*wIC;
    artifacts = reshape(artifacts,size(artifacts,1),EEG.pnts,EEG.trials);

    chaneeg   = strcmp({EEG.chanlocs.type},'EEG');
    EEG.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;

    % % Visualise
    % % vis_artifacts(EEGNEW,EEG);
    %
    NOISE = EEG;
    NOISE.data(chaneeg,:,:) = reshape(artifacts,size(artifacts,1),EEG.pnts,EEG.trials);
    vis_artifacts(EEG,NOISE);
    % % [pow, freq] = checkpowerspectrum(NOISE,1:5,[]);
else
    fprintf('Skipping wavelet thresholding since there ar no obvious bad ICs detected...\n');
    NBIC = 0;
end

%% Log
label = EEG.ALSUTRECHT.ica.ICLabel_clss(EEG.ALSUTRECHT.ica.ICLabel_cvec(ICsArtifact));
EEG.ALSUTRECHT.ica.proportionArtifactICsReducedbywICA    = NBIC./length(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAmuscle  = sum(strcmpi(label,'Muscle'));
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAchannel = sum(strcmpi(label,'Channel Noise'));
EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAeye     = sum(strcmpi(label,'Eye'));

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA cleaning\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of bad ICs:    %d\n',NBIC);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of muscle IC:  %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAmuscle);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of channel IC: %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAchannel);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of eye IC:     %d\n',EEG.ALSUTRECHT.ica.numberArtifactICsReducedbywICAeye);

end