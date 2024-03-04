function EEGNEW = do_CEEMDAN_ICA(EEG)
%
% See: RELAX_wICA_on_ICLabel_artifacts
%
% SDukic, February 2024
%

% Make sure that 'icaaact' is not empty
if isempty(EEG.icaact)==1
    EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);
end
IC = reshape(EEG.icaact, size(EEG.icaact,1), []);

% labels = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};
ICsMuscle = find(EEG.reject.gcompreject);
ICsMuscle = ICsMuscle(strcmpi(EEG.ALSUTRECHT.ica.ICLabel_clss(EEG.ALSUTRECHT.ica.ICLabel_cvec(ICsMuscle)),'Muscle'));

% Only do CEEMDAN on 'muscle' ICs identified by ICLabel
fprintf('Performing CEEMDAN thresholding on %d muscle ICs...\n',length(ICsMuscle));

% Log
fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'CEEMDAN on muscle IC\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Muscle IC: %d\n',length(ICsMuscle));

wIC = zeros(size(IC));
% wIC = IC;
for i = 1:length(ICsMuscle)
    % tic
    sig = IC(ICsMuscle(i),:);
    mxInputData      = gpuArray(single(sig)); % input1
    mxInputDataIndex = gpuArray(int32(0:length(sig)-1)); % input2
    noiseStrength    = single(0.2); % input4
    numNoise         = int32(30);   % input5
    max_iter         = int32(500);  % input6
    num_IMFs         = int32(10);   % input7
    mxOutputData = gpuArray(single(zeros(int32(length(sig)), num_IMFs))); % input3

    IMFs  = iceemdanCUDA_MexFun(mxInputData, mxInputDataIndex, mxOutputData, noiseStrength, numNoise, max_iter, num_IMFs);
    modes = gather(IMFs);

    K = mean(sum(modes,2)./sig');
    modes = modes./K;

    % % sig = IC(ICsMuscle(i),:);
    % % [modes,ort,fvs,iterNum] = emd_sssc(sig,EEG.srate,'display',1);
    % % resid = modes(end,:); modes(end,:) = [];
    % % modes = modes';
    % % % modes = flip(modes,2);
    % % num_IMFs = size(modes,2);
    %
    % % % Plot
    % % % figure; multisignalplot([IC(ICsMuscle(i),:); modes'],EEG.srate,'r');
    % % figure; multisignalplot(modes',EEG.srate,'r');
    % modes3 = reshape(modes',num_IMFs,EEG.pnts,[]);
    % modes3 = mean(modes3,3);
    % sig2 = reshape(sig,EEG.pnts,[]);
    % sig2 = mean(sig2,2);
    % % figure; multisignalplot(modes3,EEG.srate,'r');
    %
    % figure; hold on;
    % plot(EEG.times,sum(modes3,1),'LineWidth',1.2,'Color',[0.8 0.2 0.2]);
    % plot(EEG.times,sig2,'LineWidth',1.2,'Color',[0 0 0]);
    % legend('EMD','Original');
    %
    % K = mean(sum(modes3,1)./sig2');
    % figure; hold on;
    % plot(EEG.times,sum(modes3,1)./K,'LineWidth',1.2,'Color',[0.8 0.2 0.2]);
    % plot(EEG.times,sig2,'LineWidth',1.2,'Color',[0 0 0]);
    % legend('EMD','Original');
    %
    % figure; hold on;
    % plot(EEG.times,sum(modes3(2:end,:)),'LineWidth',1.2);
    % plot(EEG.times,sum(modes3(3:end,:)),'LineWidth',1.2);
    % plot(EEG.times,sum(modes3(1,:),1),'LineWidth',1.2);
    % legend('2:end','3:end','1');
    % xlim([EEG.times([1 end])]);
    %
    % tmp = [];
    % tmp.data  = [sig; modes'];
    % tmp.srate = EEG.srate;
    % [pow, freq] = checkpowerspectrum(tmp,1:10,[]);

    % Power
    [pow, freq] = powerspectrum(modes',EEG.srate);
    % figure; plot(freq,zscore(pow')');

    % Zscore and see which components have a lot of >30Hz oscillations
    pow = zscore(pow');
    powdiff = max(pow(freq>50,:))-max(pow(freq<50,:));
    mask = powdiff>0;

    fprintf(EEG.ALSUTRECHT.subject.fid,'%d: muscle artifact is estimated using the first %d modes.\n',i,sum(mask));
    wIC(ICsMuscle(i),:) = sum(modes(:,mask),2);

    % toc
end

%% Plot the ICs vs. wICs
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

%% NWB added section to remove wICA artifact and reconstruct data within this function (rather than in the main script):
artifacts = EEG.icawinv*wIC;

NOISE = EEG;
NOISE.data = artifacts;
vis_artifacts(EEG,NOISE);

% Subtract out wavelet artifact signal from EEG signal
EEGNEW = EEG;
EEGNEW.data = reshape(EEG.data,EEG.nbchan,[])-artifacts;

% Visualise
vis_artifacts(EEGNEW,EEG);

%%

EEG.ALSUTRECHT.ica.proportionArtifactICsReducedbywICA = mean(EEG.reject.gcompreject);

% =========================================================================

function [powout, freq] = powerspectrum(data,fs)
NCHN = size(data,1);

L = 4*fs;
NTRL = floor(size(data,2)/(L));
data = reshape(data(:,1:NTRL*L),NCHN,L,NTRL);

NFRQ = L/2+1;
freq = fs*(0:(L/2))/L;

winhan = hanning(L)';

powout = NaN(NCHN,NFRQ);
for i = 1:NCHN
    pow = NaN(NTRL,NFRQ);
    for j = 1:NTRL
        tmp = fft(winhan.*squeeze(data(i,:,j)));
        pow(j,:) = abs(tmp(1:NFRQ)).^2;
    end
    powout(i,:) = mean(pow,1);
end
