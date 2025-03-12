function EEG = do_ICA(EEG,cfg)

fprintf('\n================================\n');
fprintf('ICA\n');
fprintf('================================\n');

% Check the rank
rankICA  = select_nICs(EEG,cfg.ica.icMax);
rankData = get_rank(EEG.data(:,:));

if rankICA > rankData
    warning('Many electrodes were already excluded... Probably low quality data.');
    printf('Data rank %d < %d\n',rankData,rankICA);
    rankICA = rankData;
end

% Check variance explained
[U,S,V] = svd(double(EEG.data(:,:)) * double(EEG.data(:,:))','econ');
S = diag(S);
explained = S ./ sum(S);

% 95% variance
explainedTmp = cumsum(explained); % figure; bar(explained);
NPCA95 = find(explainedTmp >= 0.95,1);

% Variance used for ICA
VarRank = 100 * sum(explained(1:rankICA));

% ICA
fprintf('\n--------------------------------\n');
fprintf('Doing ICA (%s)\n',cfg.ica.type2);
fprintf('--------------------------------\n');

if strcmpi(cfg.ica.type2,'AMICA')
    % AMICA folder
    prevpath = pwd;
    if exist(EEG.ALSUTRECHT.subject.icadata,'dir')~=7, mkdir(EEG.ALSUTRECHT.subject.icadata); end
    cd(EEG.ALSUTRECHT.subject.icadata);

    % Run AMICA
    % runamica15(EEG.data,'num_chans',EEG.nbchan,'outdir',EEG.ALSUTRECHT.subject.icadata,'pcakeep',rankICA,'num_models',1,'numprocs',1,'max_threads',10,'do_reject',1,'numrej',15,'rejsig',3,'rejint',1);
    runamica15(EEG.data,'num_chans',EEG.nbchan,'outdir',EEG.ALSUTRECHT.subject.icadata,'pcakeep',rankICA,'num_models',1,'numprocs',1,'max_threads',10,'max_iter',1000);

    EEG.etc.amica   = loadmodout15(EEG.ALSUTRECHT.subject.icadata);
    EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs,:);
    EEG.icaweights  = EEG.etc.amica.W;
    EEG.icasphere   = EEG.etc.amica.S;

    % % Run multimode AMICA for RS (EO+EC) data - MT too???
    % runamica15(EEG.data,'num_chans',EEG.nbchan,'outdir',EEG.ALSUTRECHT.subject.icadata,'pcakeep',rankICA,'num_models',2,'numprocs',1,'max_threads',10,'max_iter',1000);
    %
    % model_index = 2;
    % EEG.etc.amica   = loadmodout15(EEG.ALSUTRECHT.subject.icadata);
    % EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs,:);
    % EEG.icawinv     = EEG.etc.amica.A(:,:,model_index);
    % EEG.icaweights  = EEG.etc.amica.W(:,:,model_index);
    % EEG.icasphere   = EEG.etc.amica.S;
    % EEG = eeg_checkset(EEG,'ica');
    %
    % pop_topoplot(EEG,0);
    %
    % figure; hold on;
    % plot(EEG.etc.amica.LL);
    % plot([0 1000],EEG.etc.amica.LL(800)*[1 1]);
    % plot([0 1000],EEG.etc.amica.LL(1000)*[1 1]);
    %
    % % Compute model probability
    % model_prob = 10.^EEG.etc.amica.v;
    % figure; imagesc(model_prob);
    %
    % M = movmean(model_prob,5*EEG.srate,2);
    % figure; plot((1:size(model_prob,2))/EEG.srate/60,M);

    % Previous folder
    cd(prevpath);

    % If you want to know which data points were rejected, you can check EEG.etc.amica.Lht.
    % If any datapoint shows 0, it means the datapoints were rejected by AMICA.
else
    EEG = pop_runica(EEG,'icatype',lower(cfg.ica.type2),'extended',1,'pca',rankICA,'lrate',1e-4,'maxsteps',2000);
end

% Double-check
EEG = eeg_checkset(EEG,'ica');

% Make sure IC activations are present
EEG.icaact = (EEG.icaweights*EEG.icasphere) * EEG.data(EEG.icachansind,:);
EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);

% Calculate the variances/power
NICA = length(EEG.reject.gcompreject);
varICs = NaN(NICA,1);
for i = 1:NICA
    [~, varICs(i)] = compvar(EEG.data,EEG.icaact,EEG.icawinv,i);
end

% Remove (not needed)
EEG.icaact = [];

% Log
EEG.ALSUTRECHT.ica.icmax0  = cfg.ica.icMax;
EEG.ALSUTRECHT.ica.icmax1  = rankICA;
EEG.ALSUTRECHT.ica.icmax2  = NPCA95;
EEG.ALSUTRECHT.ica.VarRank = VarRank;
EEG.ALSUTRECHT.ica.varICs  = varICs;

% Report
fprintf('\nICA is done.\n');
fprintf('Number of estimated ICs (max %d): %d\n', rankICA,max(cfg.ica.icMax));
fprintf('Used varaince: %1.2f\n',VarRank);
fprintf('95%% PCA variance: %d\n', NPCA95);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'ICA\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of estimated ICs (max %d): %d (Var = %1.2f)\n', rankICA,max(cfg.ica.icMax),VarRank);
fprintf(EEG.ALSUTRECHT.subject.fid,'95% ICs: %d\n', NPCA95);

end