function EEG = fix_chanlocs(EEG,chanlocs)
% The assumption is that the first N (= 128) electrodes are EEG
% While the rest are EMG and/or EXT electrodes
% [EEG.chanlocs] = deal(chanlocs(ismember({chanlocs.labels},{EEG(1).chanlocs.labels})));

% Define possible EXT channel labels
chansEXT = {'LEL','REL','VEOGS','VEOGI','HEOGL','HEOGR','LM','RM'};
NBLK = length(EEG);
NCHN = length(chanlocs);

for i_blk = 1:NBLK
    for i_chn = 1:NCHN
        assert(strcmpi(EEG(i_blk).chanlocs(i_chn).labels, chanlocs(i_chn).labels));
        % EEG(i_blk).chanlocs(i_chn).labels = chanlocs(i_chn).labels;
        EEG(i_blk).chanlocs(i_chn).type   = chanlocs(i_chn).type;
        EEG(i_blk).chanlocs(i_chn).theta  = chanlocs(i_chn).theta;
        EEG(i_blk).chanlocs(i_chn).radius = chanlocs(i_chn).radius;
        EEG(i_blk).chanlocs(i_chn).X      = chanlocs(i_chn).X;
        EEG(i_blk).chanlocs(i_chn).Y      = chanlocs(i_chn).Y;
        EEG(i_blk).chanlocs(i_chn).Z      = chanlocs(i_chn).Z;
        EEG(i_blk).chanlocs(i_chn).sph_theta      = chanlocs(i_chn).sph_theta;
        EEG(i_blk).chanlocs(i_chn).sph_phi        = chanlocs(i_chn).sph_phi;
        EEG(i_blk).chanlocs(i_chn).sph_radius     = chanlocs(i_chn).sph_radius;
        EEG(i_blk).chanlocs(i_chn).sph_theta_besa = chanlocs(i_chn).sph_theta_besa;
        EEG(i_blk).chanlocs(i_chn).sph_phi_besa   = chanlocs(i_chn).sph_phi_besa;
    end
end

if EEG(1).nbchan == 136
    % MMN/SART/RS
    for i_blk = 1:NBLK
        cntEXT = 0;

        for i_chn = NCHN+1:EEG(1).nbchan
            EEG(i_blk).chanlocs(i_chn).type = 'EXT';

            % DUB data sometimes has undefined EXT channel labels
            if length(EEG(i_blk).chanlocs(i_chn).labels) > 2
                % -> Not LM / RM
                if strcmpi(EEG(i_blk).chanlocs(i_chn).labels(1:3), 'EXG')
                    cntEXT = cntEXT + 1;
                    if cntEXT == 1 && i_blk == 1
                        fprintf('DUB data? Fixing external channel labels...\n');
                    end
                    if i_blk == 1
                        fprintf('%s -> ',EEG(i_blk).chanlocs(i_chn).labels);
                    end

                    % Rename
                    EEG(i_blk).chanlocs(i_chn).labels = chansEXT{cntEXT};

                    if i_blk == 1
                        fprintf('%s\n',EEG(i_blk).chanlocs(i_chn).labels);
                    end
                end
            end
        end
    end

elseif EEG(1).nbchan > 136
    % MT
    for i_blk = 1:NBLK
        for i_chn = NCHN+1:EEG(1).nbchan
            if isempty(regexp(EEG(i_blk).chanlocs(i_chn).labels,'E\d{1}','once'))
                EEG(i_blk).chanlocs(i_chn).type = 'EXT';
            else
                EEG(i_blk).chanlocs(i_chn).type = 'EMG';
            end
        end
    end
end

% Check
EEG = eeg_checkset(EEG);

% Store info of all channels
for i_blk = 1:NBLK
    EEG(i_blk).allchans = EEG(i_blk).chanlocs;
end

end