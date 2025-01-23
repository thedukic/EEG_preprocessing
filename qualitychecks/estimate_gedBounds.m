function EEG = estimate_gedBounds(EEG)
% gedBounds method to obtain empirical frequency boundaries based on spatial correlations of eigenvectors.
%
% IMPORTANT NOTE! This is designed to work for empirical EEG data in eeglab
% format. It is highly likely that you will need to make at least minor
% modifications to the code for it to work on your data.
%

%% frequency parameters
% select data, set remaining params
eegchans = strcmp({EEG.chanlocs.type},'EEG');
dataTmp  = double(EEG.data(eegchans,:));
[nbchan, pnts] = size(dataTmp);

% frequency resolution and range
numfrex  = 100;
lowfreq  = 2;  % Hz
highfreq = 80; % Hz
frex = logspace(log10(lowfreq),log10(highfreq),numfrex);

% standard deviations for the Gaussian filtering
stds = linspace(2,5,numfrex);

% onset times for epoching resting-state data
% onsets = EEG.srate*2:2*EEG.srate:pnts-EEG.srate*4;
snipn  = 2*EEG.srate;
onsets = 1:snipn:pnts;
ntrls  = length(onsets);

%% create R covariance matrix
% full R
R = zeros(ntrls,nbchan,nbchan);
for segi = 1:ntrls
    snipdat = dataTmp(:,onsets(segi):onsets(segi)+snipn-1);
    snipdat = bsxfun(@minus,snipdat,mean(snipdat,2));
    R(segi,:,:) = snipdat*snipdat' / snipn;
end

% clean R
meanR = reshape(mean(R),[],1);
dists = zeros(1,nbchan);
parfor segi = 1:nbchan
    r = R(segi,:,:);
    dists(segi) = sqrt( sum((r(:)-meanR).^2) );
end
R = squeeze(mean( R(zscore(dists)<3,:,:) ,1));

% regularized R
gamma = .01;
Rr = R*(1-gamma) + eye(nbchan)*gamma*mean(eig(R));

%% loop over frequencies
Sall = NaN(nbchan,nbchan,numfrex);
for fi = 1:numfrex
    % Filter data
    fdat = filterFGx(dataTmp,EEG.srate,frex(fi),stds(fi));

    % Compute S
    % Full S
    S = zeros(ntrls,nbchan,nbchan);
    parfor segi = 1:ntrls
        snipdat = fdat(:,onsets(segi):onsets(segi)+snipn-1);
        snipdat = bsxfun(@minus,snipdat,mean(snipdat,2));
        S(segi,:,:) = snipdat*snipdat' / snipn;
    end

    % Clean S
    meanS = reshape(mean(S),[],1);
    dists = zeros(1,size(S,1));
    parfor segi = 1:size(S,1)
        s = S(segi,:,:);
        dists(segi) = sqrt( sum((s(:)-meanS).^2) );
    end
    S = squeeze(mean( S(zscore(dists)<3,:,:) ,1));

    % Global variance normalise
    S = S / (std(S(:))/std(R(:)));
    Sall(:,:,fi) = S;
end

%% Output
EEG.ALSUTRECHT.gedBounds.Rr   = Rr;
EEG.ALSUTRECHT.gedBounds.Sall = Sall;

end