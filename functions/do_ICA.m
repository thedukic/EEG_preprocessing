function EEG = do_ICA(EEG,cfg)

% Check the rank
assert(get_rank(EEG.data(:,:))>cfg.ica.icMax);

% Check variance explained
[~,~,~,~,explained] = pca(EEG.data(:,:)');
explained = cumsum(explained./sum(explained)); % figure; bar(explained);
NICA = find(explained>0.95,1);

% ICA
EEG = pop_runica(EEG,'icatype',lower(cfg.ica.type),'extended',1,'pca',cfg.ica.icMax,'maxsteps',1000); % ,'lrate',1e-5 //// 0.00065/log(cfgica.icMax)
EEG = eeg_checkset(EEG,'ica');

% Log
EEG.ALSUTRECHT.ica.icmax  = cfg.ica.icMax;
EEG.ALSUTRECHT.ica.icmax2 = NICA;

% Report
fprintf('\nICA has just finished...\n');
fprintf('Total number of estimated ICs: %d\n', cfg.ica.icMax);
fprintf('Suggested number of useful ICs: %d\n', NICA);

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'ICA\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Total number of estimated ICs: %d\n', cfg.ica.icMax);
fprintf(EEG.ALSUTRECHT.subject.fid,'Suggested number of useful ICs: %d\n', NICA);

end