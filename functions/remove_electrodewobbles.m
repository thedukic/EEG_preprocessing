function EEG = remove_electrodewobbles(EEG,ICAtype)
%
% Removes jumps and wobbles using ICA, ICLabel and wavelet tresholding
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

% We use 2 values:
% 1st - when ICLabel is confident (P>0.5) that the it is a bad electrode IC
% 2nd - when ICLabel is unsure (P<0.5) that the it is a bad electrode IC
% Maybe good to be even more strict, to preserve more brain signals
% But this might mean that we put back some noise
% K = [1.5 2.5];
W ='coif5';
L = 5;

% Use only EEG
chaneeg = strcmp({EEG.chanlocs.type},'EEG');
EEGICA = pop_select(EEG,'channel',{EEG.chanlocs(chaneeg).labels});

% Rereference (robust)
EEGICA = do_reref(EEGICA,'aRobust');

% Data rank
assert(get_rank(EEGICA.data)==EEGICA.nbchan);

% No need to estimate all 128 ICs (?)
thisTask = EEG(1).ALSUTRECHT.subject.task;
switch thisTask
    case {'RS','EO','EC'}
        % 6 or 12 min
        NICA = 70;
    case 'SART'
        % 3 or 4 *5 = 15-20min
        NICA = 80;
    case 'MMN'
        % 3*7 = 20 min
        NICA = 80;
    case 'MT'
        % 3*7 = 20 min
        NICA = 80;
end
% NICA ^2*30 /256/60

% Check variance explained
if NICA ~= EEGICA.nbchan
    [~,~,~,~,explained] = pca(EEGICA.data(:,:)');
    explained = explained./sum(explained);
    VarRank = sum(explained(1:NICA));
else
    VarRank = 1;
end

% ICA
if NICA == EEGICA.nbchan
    EEGICA = pop_runica(EEGICA,'icatype',ICAtype,'extended',1,'lrate',1e-4,'maxsteps',2000);
else
    EEGICA = pop_runica(EEGICA,'icatype',ICAtype,'extended',1,'pca',NICA,'lrate',1e-4,'maxsteps',2000);
end
EEGICA = eeg_checkset(EEGICA,'ica');

% Make sure ICA activations are estimated
if isempty(EEGICA.icaact)
    EEGICA.icaact = (EEGICA.icaweights*EEGICA.icasphere)*EEGICA.data(EEGICA.icachansind,:);
    EEGICA.icaact = reshape(EEGICA.icaact, size(EEGICA.icaact,1), EEGICA.pnts, EEGICA.trials);
end
IC = reshape(EEGICA.icaact, size(EEGICA.icaact,1), []);

% Wavelet padding, 2^level
modulus = mod(size(IC,2),2^L);
if modulus~=0
    extra = zeros(1,(2^L)-modulus);
    IC = [IC, repmat(extra,size(IC,1),1)];
else
    extra = [];
end

% ICLabel
EEGICA = iclabel(EEGICA);

% Get the cassification labels
[ICLabel_pvec, ICLabel_cvec] = max(EEGICA.etc.ic_classification.ICLabel.classifications,[],2);

% Bad electrode ICs:
% 1. ICLabel
% badIC1 = find(ICLabel_cvec==6);
badIC1 = find(ICLabel_cvec==6 & ICLabel_pvec>0.5);
% badIC1 = unique([badIC1; find(EEGICA.etc.ic_classification.ICLabel.classifications(:,6)>=0.2)]); % BETTER DO NOT DO THIS
% badIC1 = setdiff(badIC1, find(EEGICA.etc.ic_classification.ICLabel.classifications(:,1)>=0.2));

% 2. Spatial smoothnes estimates
[spatialSmoothness, badIC2] = estimate_spatialsmoothnes(EEGICA);
muscleIC = find(ICLabel_cvec==2 & ICLabel_pvec>0.5);
badIC2 = setdiff(badIC2(:),muscleIC(:));

% Combine
badChanICs = unique([badIC1(:); badIC2(:)]);

if ~isempty(badChanICs)
    % Estimate bad electrodes in these ICs
    % The idea is that a true "bad channel" IC should usually has only one electrode with a higher weight, eg >3STD
    % Tho this might not be the case as 2-3 electrodes can jump at the same time
    % icawinv = abs(zscore(EEGICA.icawinv(:,badChanICs),0,1));
    icawinv = abs(robust_zscore(EEGICA.icawinv(:,badChanICs)));
    mask    = icawinv>5;
    % goodic  = sum(mask)==1;
    % mask    = any(mask(:,goodic),2);

    [a,b] = find(mask);
    badChans = {EEGICA.chanlocs(unique(a)).labels};

    % Report
    myCmap = brewermap(128,'*RdBu');
    icLabels = EEGICA.etc.ic_classification.ICLabel.classes;

    fh = figure; th = tiledlayout('flow');
    th.TileSpacing = 'compact'; th.Padding = 'compact';

    NBICS = length(badChanICs);
    myYlabel = cell(1,NBICS);
    for i = 1:NBICS
        nexttile;
        topoplot(EEGICA.icawinv(:,badChanICs(i)),EEGICA.chanlocs,'maplimits',max(abs(EEGICA.icawinv(:,badChanICs(i))))*[-1 1],'headrad','rim','whitebk','on','style','map','electrodes','on','shading','interp');
        myYlabel{i} = {['ICA' num2str(badChanICs(i))], [icLabels{ICLabel_cvec(badChanICs(i))} ', P = ' num2str(round(ICLabel_pvec(badChanICs(i)),2))]};
        title(myYlabel{i}); axis tight; colormap(myCmap);
    end

    % Save
    plotX=35; plotY=20;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_wobbleIC1']),'-dtiff','-r200');
    close(fh);

    % fh = figure;
    % th = tiledlayout(NBICS,4); th.TileSpacing = 'compact'; th.Padding = 'compact';
    % colormap(fh,myCmap);
    % LEpoch = EEG.srate;
    % NEpoch = floor(size(IC,2)/LEpoch);

    % Get the weights and determine the channel that is problematic
    wIC = zeros(size(IC));
    for i = 1:NBICS
        [thresh,sorh,~] = ddencmp('den','wv',IC(badChanICs(i),:));

        % Adjusting the treshold
        % Bad electrode ICs are mostly with low power/variance, e.g. > 20th IC
        % So even if the trheshold is too low
        % (i.e. higher trheshold excludes less data, by keeping higher freq)
        % then we wont lose too much brain/good data
        % if ICLabel_pvec(badChanICs(i))>0.5
        %     % Playing safe and trying to avoid losing too much good/brain data
        %     thresh = thresh.*K(1);
        % else
        %     % If ICLabel is less certain, probably there is more
        %     % brain/good data here, so further increase the treshold
        %     thresh = thresh.*K(2);
        % end
        % thresh = 0.9*thresh;

        swc = swt(IC(badChanICs(i),:),L,W);
        Y   = wthresh(swc,sorh,thresh);
        wIC(badChanICs(i),:) = iswt(Y,W);

        % tmp = [];
        % tmp.data  = wIC(badChanICs(i),:);
        % tmp.srate = EEG.srate;
        % [pow, freq] = checkpowerspectrum(tmp,1,[]);

        % % Apply only on targeted segments?
        % ICsz = abs(robust_zscore(IC(badChanICs(i),:)))<3;
        % ICsz = reshape(ICsz,LEpoch,NEpoch)';
        % ICs_all    = reshape(IC(badChanICs(i),:),LEpoch,NEpoch)';
        % ICs_noise  = reshape(wIC(badChanICs(i),:),LEpoch,NEpoch)';
        % ICs_signal = ICs_all-ICs_noise;
        % myClim = 0.9*max(abs(ICs_signal(:)));

        % nexttile; topoplot(EEGICA.icawinv(:,badChanICs(i)),EEGICA.chanlocs,'maplimits',max(abs(EEGICA.icawinv(:,badChanICs(i))))*[-1 1],'headrad','rim','whitebk','on','style','map','electrodes','on','shading','interp');
        % axis tight; colormap(myCmap);
        % th = nexttile; imagesc(th,ICs_all); axis off; clim(th,myClim*[-1 1]);
        % % text(th.Position(1),mean(th.Position([2 4])),myYlabel{i},'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10);
        % ICs_noise(ICsz) = 0;
        % th = nexttile; imagesc(th,ICs_noise); axis off; clim(th,myClim*[-1 1]);
        % th = nexttile; imagesc(th,ICs_signal); axis off; clim(th,myClim*[-1 1]);
    end

    % figure; hold on; plot(ICs0(2,:)); plot(ICs1(2,:));

    % Remove extra padding
    if ~isempty(extra)
        wIC = wIC(:,1:end-numel(extra));
    end

    % Channel-level artifacts
    artifacts = EEGICA.icawinv*wIC;

    % EEGNEW = EEG;
    % EEGNEW.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;
    % vis_artifacts(EEGNEW,EEG);

    % Subtract out wavelet artifact signal from EEG signal
    EEG.data(chaneeg,:,:) = EEG.data(chaneeg,:,:)-artifacts;

    % % Save
    % plotX=35; plotY=55;
    % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    % print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_wobbleIC2']),'-dtiff','-r300');
    % close(fh);
else
    NBICS = 0;
    badChans = {};
end

% % Visual check
% EEGNEW = pop_select(EEG,'channel',{EEG.chanlocs(chaneeg).labels});
% vis_artifacts(EEGNEW,EEGICA);

% Remove IC signals
EEG.icaact = [];

% Log
EEG.ALSUTRECHT.badchaninfo.wica.icmax   = NICA;
EEG.ALSUTRECHT.badchaninfo.wica.VarRank = VarRank;
EEG.ALSUTRECHT.badchaninfo.wica.ics     = badChanICs;
EEG.ALSUTRECHT.badchaninfo.wica.pvec    = ICLabel_pvec;
EEG.ALSUTRECHT.badchaninfo.wica.cvec    = ICLabel_cvec;
EEG.ALSUTRECHT.badchaninfo.wica.fixed   = badChans;

% Report
fprintf('\nwICA was just now used to fix electrode wobbles and jumps...\n');
fprintf('Number of estimated ICs: %d\n', NICA);
fprintf('Used varaince: %1.2f\n', VarRank);
fprintf('Detected number of bad-electrode ICs: %d\n', NBICS);
str = arrayfun(@(x) num2str(x,'%1.1f'),ICLabel_pvec(badChanICs),'uni',0);
str = strjoin(str,', ');
fprintf('Their ICLabel P-values: %s\n', str);
fprintf('The estimated electrodes that are affected: %d\n', length(badChans));

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'wICA bad electrodes\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Total number of estimated ICs: %d\n', NICA);
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected number of bad-electrode ICs: %d\n', NBICS);
fprintf(EEG.ALSUTRECHT.subject.fid,'Their ICLabel P-values: %s\n', str);
fprintf(EEG.ALSUTRECHT.subject.fid,'The estimated electrodes that are affected: %d\n', length(badChans));

end