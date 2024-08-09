function EEG = do_ICA(EEG,cfg)

% Check the rank
rankData = get_rank(EEG.data(:,:));
rankICA  = cfg.ica.icMax;

if rankICA>rankData
    warning('Many electrodes were already excluded.. Probably low quality data.');
    printf('Data rank %d < %d\n',rankData,cfg.ica.icMax);
    rankICA = rankData;
end
% assert(get_rank(EEG.data(:,:))>cfg.ica.icMax);

% Check variance explained
[~,~,~,~,explained] = pca(EEG.data(:,:)');
explained = cumsum(explained./sum(explained)); % figure; bar(explained);
NICA95 = find(explained>0.95,1);

% ICA
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

% Make sure IC activations are present
EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1), EEG.pnts, EEG.trials);

% Calculate the variances
NICA = length(EEG.reject.gcompreject);
varICs = NaN(NICA,1);
for i = 1:NICA
    [~, varICs(i)] = compvar(EEG.data,EEG.icaact,EEG.icawinv,i);
end

% Double-check
EEG = eeg_checkset(EEG,'ica');

% Log
EEG.ALSUTRECHT.ica.icmax0 = cfg.ica.icMax;
EEG.ALSUTRECHT.ica.icmax1 = rankICA;
EEG.ALSUTRECHT.ica.icmax2 = NICA95;
EEG.ALSUTRECHT.ica.varICs = varICs;

% Report
fprintf('\nICA has just finished...\n');
fprintf('Total number of estimated ICs (max %d): %d\n', rankICA,cfg.ica.icMax);
fprintf('Suggested number of useful ICs: %d\n', NICA95);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'ICA\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Total number of estimated ICs (max %d): %d\n', rankICA,cfg.ica.icMax);
fprintf(EEG.ALSUTRECHT.subject.fid,'Suggested number of useful ICs: %d\n', NICA95);

end