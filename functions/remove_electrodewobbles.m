function EEG = remove_electrodewobbles(EEG,ICAtype)

% Use only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
EEGICA = pop_select(EEG,'channel',{EEG.chanlocs(chaneeg).labels});

% Rereference (robust)
EEGICA.data = EEGICA.data - trimmean(EEGICA.data,10,'round',1);

% Data rank
assert(getrank(EEGICA.data)==EEGICA.nbchan);

% ICA
% EEGICA = pop_runica(EEGICA,'icatype','cudaica','pca',70);
EEGICA = pop_runica(EEGICA,'icatype',ICAtype);
EEGICA = eeg_checkset(EEGICA,'ica');

% Make sure ICA activations are estimated
if isempty(EEGICA.icaact)
    EEGICA.icaact = (EEGICA.icaweights*EEGICA.icasphere)*EEGICA.data(EEGICA.icachansind,:);
    EEGICA.icaact = reshape(EEGICA.icaact, size(EEGICA.icaact,1), EEGICA.pnts, EEGICA.trials);
end
IC = reshape(EEGICA.icaact, size(EEGICA.icaact,1), []);

% Wavelet padding, 2^level
WLT ='coif5';
LVL = 5;
modulus = mod(size(IC,2),2^LVL);
if modulus~=0
    extra = zeros(1,(2^LVL)-modulus);
    IC = [IC, repmat(extra,size(IC,1),1)];
else
    extra = [];
end

% ICLabel
EEGICA = iclabel(EEGICA);

% Interim data saving
% EEGICA = pop_saveset(EEGICA,'filename','TMP.set','filepath',EEGICA.ALSUTRECHT.subject.preproc);

% Get the cassification labels
[ICLabel_pvec, ICLabel_cvec] = max(EEGICA.etc.ic_classification.ICLabel.classifications,[],2);

% (6) channel noise
badChanICs = find(ICLabel_cvec==6);
% badChanICs = find(ICLabel_cvec==6 & ICLabel_pvec>0.5);
badChanICs = unique([badChanICs; find(EEGICA.etc.ic_classification.ICLabel.classifications(:,6)>=0.2)]);
% badChanICs = setdiff(badChanICs,find(EEGICA.etc.ic_classification.ICLabel.classifications(:,1)>=0.2));

% The idea is that a true "bad channel" IC should have only one
% electrode with a higher weight, eg >3STD
icawinv = abs(zscore(EEGICA.icawinv(:,badChanICs),0,1));
mask    = icawinv>5;
% goodic  = sum(mask)==1;
% mask    = any(mask(:,goodic),2);

[a,b] = find(mask);
badChans = {EEGICA.chanlocs(unique(a)).labels};

% % Report
% myCmap = brewermap([],'*RdBu');
% iclabels = EEGICA.etc.ic_classification.ICLabel.classes;
%
% fh = figure; tiledlayout('flow');
% for i = 1:length(badChanICs)
%     nexttile;
%     topoplot(EEGICA.icawinv(:,badChanICs(i)),EEGICA.chanlocs,'maplimits',0.9*max(abs(EEGICA.icawinv(:,badChanICs(i))))*[-1 1],'headrad','rim','whitebk','on','style','map','electrodes','off');
%     title({['ICA' num2str(badChanICs(i))], [iclabels{ICLabel_cvec(badChanICs(i))} ', P = ' num2str(round(ICLabel_pvec(badChanICs(i)),2))]}); axis tight; colormap(myCmap);
% end
% % EEGICA = pop_saveset(EEGICA,'filename','TMP.set','filepath','C:\Users\Stefan\OneDrive\Desktop');

% Get the weights and determine the channel that is problematic
if ~isempty(badChanICs)
    wIC = zeros(size(IC));
    for i = 1:length(badChanICs)
        [thresh,sorh,~] = ddencmp('den','wv',IC(badChanICs(i),:));
        thresh = thresh*1.2;
        swc = swt(IC(badChanICs(i),:),LVL,WLT);
        Y   = wthresh(swc,sorh,thresh);
        wIC(badChanICs(i),:) = iswt(Y,WLT);

        % tmp = [];
        % tmp.data  = [sig; wIC(badChanICs(i),:)];
        % tmp.srate = EEG.srate;
        % [pow, freq] = checkpowerspectrum(tmp,1:2,[]);
    end

    % Remove extra padding
    if ~isempty(extra)
        wIC = wIC(:,1:end-numel(extra));
    end

    % Channel-level artifacts
    artifacts = EEGICA.icawinv*wIC;

    % Subtract out wavelet artifact signal from EEG signal
    EEG.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;

    % EEGNEW = EEG;
    % EEGNEW.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;
    % vis_artifacts(EEGNEW,EEG);
end

% Log
EEG.ALSUTRECHT.badchaninfo.wica.ics   = badChanICs;
EEG.ALSUTRECHT.badchaninfo.wica.pvec  = ICLabel_pvec;
EEG.ALSUTRECHT.badchaninfo.wica.cvec  = ICLabel_cvec;
EEG.ALSUTRECHT.badchaninfo.wica.fixed = badChans;

fprintf('wICA bad electrode ICs: %d\n', length(badChanICs));
str = arrayfun(@(x) num2str(x,'%1.1f'),ICLabel_pvec(badChanICs),'uni',0);
str = strjoin(str,', ');
fprintf('wICA P-values: %s\n', str);
fprintf('wICA likely electrodes: %d\n', length(badChans));

end