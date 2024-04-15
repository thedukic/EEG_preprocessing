function EEG = detect_elecpopsanddrifts(EEG)

fprintf('\nUsing HEAR method to fix electrode pops and drifts...\n');

chaneeg = strcmp({EEG.chanlocs.type},'EEG');
chanlocs = EEG.chanlocs(chaneeg);
D = utl_chaninterpmatrix(chanlocs,6);

% fit HEAR to the calibration data
hear_mdl = HEAR(EEG.srate,true,[],[],D);
hear_mdl.train(EEG.data(chaneeg,:));

% detect and remove transient, high variance artifacts
[p_art, p_confidence, data2] = hear_mdl.apply(EEG.data(chaneeg,:));

EEG0 = EEG;
EEG.data(chaneeg,:) = data2;
vis_artifacts(EEG0,EEG);

fprintf('Done.\n');