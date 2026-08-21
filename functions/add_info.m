function EEG = add_info(EEG, subject, cfg)

fprintf('\n================================\n');
fprintf('Adding participant / channel info\n');
fprintf('================================\n');

% Load standard 128-channel BioSemi montage
chanlocs = readlocs('biosemi128_eeglab.ced');

% Define external and muscle channel labels
labels_ext = {'VEOGS', 'VEOGI', 'HEOGL', 'HEOGR', 'ECGL', 'ECGR', 'LM', 'RM'};
labels_emg = {'APB', 'FDI', 'FPB', 'EPB', 'EDC', 'FDS'};

% Fix chanlocs field
EEG = fix_chanlocs(EEG, chanlocs, labels_ext, labels_emg);

% Attach experimental configuration and metadata
num_blocks = length(EEG);
for i_block = 1:num_blocks
    EEG(i_block).ALSUTRECHT.subject = subject;
    EEG(i_block).ALSUTRECHT.cfg     = cfg;
end

fprintf('Channel metadata and subject info successfully added.\n');

end


function EEG = fix_chanlocs(EEG, chanlocs, labels_ext, labels_emg)
% FIX_CHANLOCS
% Assigns 3D coordinates to channels 1:128 (EEG) and maps generic BioSemi
% external channels (EXG1..EXG8) and EMG leads (EXG9..EXG20 / E1..E12)
% to standardised anatomical labels and channel types.

num_blocks   = length(EEG);
num_eeg_chan = length(chanlocs); % Typically 128

for i_block = 1:num_blocks
    total_chan = EEG(i_block).nbchan;

    % ---------------------------------------------------------------------
    % 1. Map 128 Scalp EEG Channels
    % ---------------------------------------------------------------------
    for i_channel = 1:num_eeg_chan
        % Validate label alignment if labels already exist
        if ~isempty(EEG(i_block).chanlocs(i_channel).labels) && ...
                ~strncmpi(EEG(i_block).chanlocs(i_channel).labels, 'EXG', 3)
            assert(strcmpi(EEG(i_block).chanlocs(i_channel).labels, chanlocs(i_channel).labels), ...
                'Mismatch between EEG dataset label (%s) and CED file (%s) at channel %d.', ...
                EEG(i_block).chanlocs(i_channel).labels, chanlocs(i_channel).labels, i_channel);
        end

        % Copy coordinates and set type
        EEG(i_block).chanlocs(i_channel).labels         = chanlocs(i_channel).labels;
        EEG(i_block).chanlocs(i_channel).type           = 'EEG';
        EEG(i_block).chanlocs(i_channel).theta          = chanlocs(i_channel).theta;
        EEG(i_block).chanlocs(i_channel).radius         = chanlocs(i_channel).radius;
        EEG(i_block).chanlocs(i_channel).X              = chanlocs(i_channel).X;
        EEG(i_block).chanlocs(i_channel).Y              = chanlocs(i_channel).Y;
        EEG(i_block).chanlocs(i_channel).Z              = chanlocs(i_channel).Z;
        EEG(i_block).chanlocs(i_channel).sph_theta      = chanlocs(i_channel).sph_theta;
        EEG(i_block).chanlocs(i_channel).sph_phi        = chanlocs(i_channel).sph_phi;
        EEG(i_block).chanlocs(i_channel).sph_radius     = chanlocs(i_channel).sph_radius;
        EEG(i_block).chanlocs(i_channel).sph_theta_besa = chanlocs(i_channel).sph_theta_besa;
        EEG(i_block).chanlocs(i_channel).sph_phi_besa   = chanlocs(i_channel).sph_phi_besa;
    end

    % ---------------------------------------------------------------------
    % 2. Map External Channels (129:136) -> EOG / ECG / Mastoids
    % ---------------------------------------------------------------------
    if total_chan >= 136
        cnt_ext = 0;
        for i_channel = (num_eeg_chan + 1):min(total_chan, num_eeg_chan + 8)
            cnt_ext = cnt_ext + 1;

            % Assign standard name if channel has a generic hardware label
            cur_label = EEG(i_block).chanlocs(i_channel).labels;
            if isempty(cur_label) || strncmpi(cur_label, 'EXG', 3)
                EEG(i_block).chanlocs(i_channel).labels = labels_ext{cnt_ext};
            end

            % Assign specific types based on assigned label
            EEG(i_block).chanlocs(i_channel).type = 'EXT';
            % cur_label = upper(cur_label);
            % if contains(cur_label, 'EOG')
            %     EEG(i_block).chanlocs(i_channel).type = 'EOG';
            % elseif contains(cur_label, 'ECG')
            %     EEG(i_block).chanlocs(i_channel).type = 'ECG';
            % else
            %     EEG(i_block).chanlocs(i_channel).type = 'EXT';
            % end
        end
    end

    % ---------------------------------------------------------------------
    % 3. Map EMG Channels (137:148) -> 6 Muscle Pairs (12 Leads)
    % ---------------------------------------------------------------------
    if total_chan > 136
        num_leads = total_chan - 136;
        num_pairs = num_leads / 2;
        assert(length(labels_emg) == num_pairs);

        cnt_emg = 0;
        for i_channel = 1:num_pairs
            cnt_emg = cnt_emg + 1;
            cur_label = labels_emg{cnt_emg};

            lead1_idx = 136 + 2*i_channel - 1;
            lead2_idx = 136 + 2*i_channel;

            % Assign distinctive monopolar labels
            EEG(i_block).chanlocs(lead1_idx).labels = [cur_label, '1'];
            EEG(i_block).chanlocs(lead1_idx).type   = 'EMG';

            EEG(i_block).chanlocs(lead2_idx).labels = [cur_label, '2'];
            EEG(i_block).chanlocs(lead2_idx).type   = 'EMG';
        end
    end

    % Validate and cache full channel backup
    EEG(i_block) = eeg_checkset(EEG(i_block));
    EEG(i_block).allchans = EEG(i_block).chanlocs;
end

end