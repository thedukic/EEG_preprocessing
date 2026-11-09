function EEG_filt = do_filtering_method(EEG, cut_off, cfg)
% Filter data just for the method (GEDAI / ICA)

fprintf('Preparing temporary %.2f Hz filtered dataset...\n', cut_off);

% Ensure that this is applied to all possible tasks that could be processed
% Can be improved by changing all avialble fields under .filt 
cfg_tmp = cfg;
cfg_tmp.flt.erp.hp = [cut_off, 4];
cfg_tmp.flt.rs.hp  = [cut_off, 4];
cfg_tmp.flt.mt.hp  = [cut_off, 4];

% Filter
EEG_filt = do_filtering(EEG, 'highpass', 'eeg', cfg_tmp);

end