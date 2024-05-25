function EEG = detect_badICs(EEG,EXT,cfg)

% ICLabel
EEG = iclabel(EEG);

% Flag artifact ICs
EEG = pop_icflag(EEG,cfg.ica.iclabel);

% Log IClabel info
EEG.ALSUTRECHT.ica.ICLabel_bics = find(EEG.reject.gcompreject);
EEG.ALSUTRECHT.ica.ICLabel_clss = EEG.etc.ic_classification.ICLabel.classes;
[EEG.ALSUTRECHT.ica.ICLabel_pvec, EEG.ALSUTRECHT.ica.ICLabel_cvec] = max(EEG.etc.ic_classification.ICLabel.classifications,[],2);

% Manually add heart ICs
K = find(contains(EEG.ALSUTRECHT.ica.ICLabel_clss,'Heart'));
ICHeart = EEG.ALSUTRECHT.ica.ICLabel_cvec==K;
if any(ICHeart)
    if EEG.ALSUTRECHT.ica.ICLabel_pvec(ICHeart)>0.5
        EEG.reject.gcompreject(ICHeart) = true;
    end
end

% Update the log
EEG.ALSUTRECHT.ica.ICLabel_bics = find(EEG.reject.gcompreject);

end