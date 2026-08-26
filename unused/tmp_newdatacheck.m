
subj_bulk = list_participants("\\Ds.umcutrecht.nl\data\HER\onderzoeksarchief\19-462_ALS-Electrode_BS\E_ResearchData\2_ResearchData\2_PREPROCESSED\RS\ALS\T1",[]);
subj_bulk(ismember(subj_bulk,'reports')) = [];

subj_pc = list_participants("C:\DATA\MATLAB\EEG\1_EEG_DATA\ALS\T1",[]);

subj_bulk(~ismember(subj_bulk,subj_pc))