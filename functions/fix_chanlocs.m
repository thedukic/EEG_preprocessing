function EEG = fix_chanlocs(EEG,chanlocs)
% The assumption is that the first N (= 128) electrodes are EEG
% While the rest are EMG and/or EXT electrodes
% [EEG.chanlocs] = deal(chanlocs(ismember({chanlocs.labels},{EEG(1).chanlocs.labels})));

NBLK = length(EEG);
NCHN = length(chanlocs);
for i = 1:NBLK
    for j = 1:NCHN
        assert(strcmpi(EEG(i).chanlocs(j).labels,chanlocs(j).labels));
        % EEG(i).chanlocs(j).labels = chanlocs(j).labels;
        EEG(i).chanlocs(j).type   = chanlocs(j).type;
        EEG(i).chanlocs(j).theta  = chanlocs(j).theta;
        EEG(i).chanlocs(j).radius = chanlocs(j).radius;
        EEG(i).chanlocs(j).X      = chanlocs(j).X;
        EEG(i).chanlocs(j).Y      = chanlocs(j).Y;
        EEG(i).chanlocs(j).Z      = chanlocs(j).Z;
        EEG(i).chanlocs(j).sph_theta      = chanlocs(j).sph_theta;
        EEG(i).chanlocs(j).sph_phi        = chanlocs(j).sph_phi;
        EEG(i).chanlocs(j).sph_radius     = chanlocs(j).sph_radius;
        EEG(i).chanlocs(j).sph_theta_besa = chanlocs(j).sph_theta_besa;
        EEG(i).chanlocs(j).sph_phi_besa   = chanlocs(j).sph_phi_besa;
    end
end

if EEG(1).nbchan==136
    % MMN/SART/RS
    for i = 1:NBLK
        for j = NCHN+1:EEG(1).nbchan
            EEG(i).chanlocs(j).type = 'EXT';
        end
    end
elseif EEG(1).nbchan>136
    % MT
    for i = 1:NBLK
        for j = NCHN+1:EEG(1).nbchan
            if isempty(regexp(EEG(i).chanlocs(j).labels,'E\d{1}','once'))
                EEG(i).chanlocs(j).type = 'EXT';
            else
                EEG(i).chanlocs(j).type = 'EMG';
            end
        end
    end
end

% Check
EEG = eeg_checkset(EEG);

for i = 1:NBLK
    EEG(i).allchans = EEG(i).chanlocs;
end

end