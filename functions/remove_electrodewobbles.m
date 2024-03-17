function EEG = remove_electrodewobbles(EEG,ICAtype)

% Use only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
EEGICA = pop_select(EEG,'channel',{EEG.chanlocs(chaneeg).labels});

% Rereference (robust)
EEGICA.data = EEGICA.data - trimmean(EEGICA.data,10,'round',1);

% Data rank
assert(get_rank(EEGICA.data)==EEGICA.nbchan);

% No need to estimate all ICs
[~,~,~,~,explained] = pca(EEGICA.data');
explained = cumsum(explained./sum(explained)); % figure; bar(explained);
NICA = find(explained>0.99,1);
% NICA = min(NICA,100);
% NICA = max(NICA,50);

% ICA
% EEGICA = pop_runica(EEGICA,'icatype',ICAtype);
EEGICA = pop_runica(EEGICA,'icatype',ICAtype,'extended',1,'pca',NICA);
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
% badChanICs = unique([badChanICs; find(EEGICA.etc.ic_classification.ICLabel.classifications(:,6)>=0.2)]); % BETTER DO NOT DO THIS
% badChanICs = setdiff(badChanICs,find(EEGICA.etc.ic_classification.ICLabel.classifications(:,1)>=0.2));

% The idea is that a true "bad channel" IC should usually have only one
% electrode with a higher weight, eg >3STD
% Tho this might not be the case as 2-3 electrodes can jump at the same time
icawinv = abs(zscore(EEGICA.icawinv(:,badChanICs),0,1));
mask    = icawinv>5;
% goodic  = sum(mask)==1;
% mask    = any(mask(:,goodic),2);

[a,b] = find(mask);
badChans = {EEGICA.chanlocs(unique(a)).labels};

% Report
myCmap = brewermap([],'*RdBu');
icLabels = EEGICA.etc.ic_classification.ICLabel.classes;

fh = figure; th = tiledlayout('flow');
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:length(badChanICs)
    nexttile;
    topoplot(EEGICA.icawinv(:,badChanICs(i)),EEGICA.chanlocs,'maplimits',max(abs(EEGICA.icawinv(:,badChanICs(i))))*[-1 1],'headrad','rim','whitebk','on','style','map','electrodes','off');
    title({['ICA' num2str(badChanICs(i))], [icLabels{ICLabel_cvec(badChanICs(i))} ', P = ' num2str(round(ICLabel_pvec(badChanICs(i)),2))]}); axis tight; colormap(myCmap);
end

% Save
plotX=35; plotY=20;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_ICwobbles']),'-dtiff','-r200');
close(fh);

% EEGICA = pop_saveset(EEGICA,'filename','TMP.set','filepath','C:\Users\Stefan\OneDrive\Desktop');

% Get the weights and determine the channel that is problematic
if ~isempty(badChanICs)
    wIC = zeros(size(IC));
    for i = 1:length(badChanICs)
        [thresh,sorh,~] = ddencmp('den','wv',IC(badChanICs(i),:));

        % Adjusting the treshold
        % Luckily channel ICs are mostly with low power/variance in general
        % (e.g. > 20th IC)
        % So even if the trheshold is too low
        % (i.e. higher trheshold excludes less data, by keeping higher freq)
        % then we wont lose too much brain/good data
        if ICLabel_pvec(badChanICs(i))>0.5
            % Playing safe and trying to avoid losing too much good/brain data
            thresh = thresh*1.2;
        else
            % If ICLabel is less certain, probably there is more
            % brain/good data here, so further increase the treshold
            thresh = thresh*1.5;
        end
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
EEG.ALSUTRECHT.badchaninfo.wica.icmax = NICA;
EEG.ALSUTRECHT.badchaninfo.wica.ics   = badChanICs;
EEG.ALSUTRECHT.badchaninfo.wica.pvec  = ICLabel_pvec;
EEG.ALSUTRECHT.badchaninfo.wica.cvec  = ICLabel_cvec;
EEG.ALSUTRECHT.badchaninfo.wica.fixed = badChans;

% Report
fprintf('\nwICA was just now used to fix electrode wobbles and jumps...\n');
fprintf('Total number of estimated ICs: %d\n', NICA);
fprintf('Detected number of bad-electrode ICs: %d\n', length(badChanICs));
str = arrayfun(@(x) num2str(x,'%1.1f'),ICLabel_pvec(badChanICs),'uni',0);
str = strjoin(str,', ');
fprintf('Their ICLabel P-values: %s\n', str);
fprintf('The estimated electrodes that are affected: %d\n', length(badChans));

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA bad electrodes\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Total number of estimated ICs: %d\n', NICA);
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected number of bad-electrode ICs: %d\n', length(badChanICs));
fprintf(EEG.ALSUTRECHT.subject.fid,'Their ICLabel P-values: %s\n', str);
fprintf(EEG.ALSUTRECHT.subject.fid,'The estimated electrodes that are affected: %d\n', length(badChans));

end