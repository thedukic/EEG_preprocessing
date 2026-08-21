function DATA = make_extbipolar(DATA)
% MAKE_EXTBIPOLAR
% Converts paired external monopolar electrodes (ECG, VEOG, HEOG, EMG)
% into bipolar channels and removes redundant reference leads.
%
% Strictly requires VEOG (VEOGS/VEOGI) and HEOG (HEOGL/HEOGR).

do_bipolar_emg = false; % Set to true to derive bipolar EMG pairs

fprintf('\n==================================================\n');
fprintf('Making External & EMG channels bipolar\n');
fprintf('==================================================\n');

% if do_bipolar_emg
%     fprintf('Mode: EMG -> Bipolar derivation\n');
% else
%     fprintf('Mode: EMG -> Retain monopolar\n');
% end

num_block = length(DATA);

for i_block = 1:num_block
    if num_block > 1
        fprintf('\n--- Block %d of %d ---\n', i_block, num_block);
    end

    eleclabels = {DATA(i_block).chanlocs.labels};
    nbchan_old = DATA(i_block).nbchan;
    assert(length(eleclabels) == nbchan_old);

    % ---------------------------------------------------------------------
    % 1. Bipolar ECG Derivation (Dedicated ECGL/ECGR or Earlobe LEL/REL)
    % ---------------------------------------------------------------------
    eleclabels = {DATA(i_block).chanlocs.labels};

    idx_ecgl = find(strcmpi(eleclabels, 'ECGL'), 1);
    idx_ecgr = find(strcmpi(eleclabels, 'ECGR'), 1);
    idx_lel  = find(strcmpi(eleclabels, 'LEL'), 1);
    idx_rel  = find(strcmpi(eleclabels, 'REL'), 1);

    if ~isempty(idx_ecgl) && ~isempty(idx_ecgr)
        % Dedicated bipolar ECG
        DATA(i_block).data(idx_ecgl, :, :) = DATA(i_block).data(idx_ecgl, :, :) - DATA(i_block).data(idx_ecgr, :, :);
        DATA(i_block).chanlocs(idx_ecgl).labels = 'ECG';
        DATA(i_block).chanlocs(idx_ecgl).type   = 'EXT';

        % Remove redundant reference channel
        DATA(i_block).data(idx_ecgr, :, :) = [];
        DATA(i_block).chanlocs(idx_ecgr)   = [];
        eleclabels = {DATA(i_block).chanlocs.labels};

        % Metadata logging
        DATA(i_block).ALSUTRECHT.subject.ecg = 'recorded';

        fprintf('  [ECG]  True bipolar ECG derived (ECGL - ECGR -> ECG).\n');

    elseif ~isempty(idx_lel) && ~isempty(idx_rel)
        % Extract 2-channel earlobe data across all continuous samples
        lel_data = squeeze(DATA(i_block).data(idx_lel, :, :));
        rel_data = squeeze(DATA(i_block).data(idx_rel, :, :));

        % Flatten if data is epoched [samples x 2]
        ear_mat = [lel_data(:), rel_data(:)];

        % Extract 1st Principal Component (maximises shared cardiac variance)
        [~, score] = pca(ear_mat);
        ecg_approx = score(:, 1);

        % Enforce positive R-peak polarity (R-peaks generate positive skewness)
        if skewness(ecg_approx) < 0
            ecg_approx = -ecg_approx;
        end

        % Reshape back to original dimensions
        DATA(i_block).data(idx_lel, :, :) = reshape(ecg_approx, size(DATA(i_block).data(idx_lel, :, :)));
        DATA(i_block).chanlocs(idx_lel).labels = 'ECG';
        DATA(i_block).chanlocs(idx_lel).type   = 'EXT';

        % Remove redundant reference channel
        DATA(i_block).data(idx_rel, :, :) = [];
        DATA(i_block).chanlocs(idx_rel)   = [];
        eleclabels = {DATA(i_block).chanlocs.labels};

        % Metadata logging
        DATA(i_block).ALSUTRECHT.subject.ecg = 'approximated';

        fprintf('  [ECG]  Approximated bipolar ECG derived from earlobes (LEL - REL -> ECG).\n');

    else
        error('Something went wrong.');
    end

    % ---------------------------------------------------------------------
    % 2. Bipolar VEOG (Mandatory: VEOGS - VEOGI)
    % ---------------------------------------------------------------------
    m1_veog = find(strcmpi(eleclabels, 'VEOGS'));
    m2_veog = find(strcmpi(eleclabels, 'VEOGI'));

    if isempty(m1_veog) || isempty(m2_veog)
        error('Block %d: Missing mandatory VEOG channels. Found VEOGS: %d, VEOGI: %d', ...
            i_block, ~isempty(m1_veog), ~isempty(m2_veog));
    end

    DATA(i_block).data(m1_veog, :, :) = DATA(i_block).data(m1_veog, :, :) - DATA(i_block).data(m2_veog, :, :);
    DATA(i_block).chanlocs(m1_veog).labels = 'VEOG';
    DATA(i_block).chanlocs(m1_veog).type   = 'EXT';

    DATA(i_block).data(m2_veog, :, :) = [];
    DATA(i_block).chanlocs(m2_veog)   = [];

    fprintf('  [VEOG] Derived bipolar VEOG (VEOGS - VEOGI -> VEOG).\n');

    % Refresh label list after deletion
    eleclabels = {DATA(i_block).chanlocs.labels};

    % ---------------------------------------------------------------------
    % 3. Bipolar HEOG (Mandatory: HEOGL - HEOGR)
    % ---------------------------------------------------------------------
    m1_heog = find(strcmpi(eleclabels, 'HEOGL'));
    m2_heog = find(strcmpi(eleclabels, 'HEOGR'));

    if isempty(m1_heog) || isempty(m2_heog)
        error('Block %d: Missing mandatory HEOG channels. Found HEOGL: %d, HEOGR: %d', ...
            i_block, ~isempty(m1_heog), ~isempty(m2_heog));
    end

    DATA(i_block).data(m1_heog, :, :) = DATA(i_block).data(m1_heog, :, :) - DATA(i_block).data(m2_heog, :, :);
    DATA(i_block).chanlocs(m1_heog).labels = 'HEOG';
    DATA(i_block).chanlocs(m1_heog).type   = 'EXT';

    DATA(i_block).data(m2_heog, :, :) = [];
    DATA(i_block).chanlocs(m2_heog)   = [];

    fprintf('  [HEOG] Derived bipolar HEOG (HEOGL - HEOGR -> HEOG).\n');

    % ---------------------------------------------------------------------
    % 4. EMG Channel Handling
    % ---------------------------------------------------------------------
    labels = {DATA(i_block).chanlocs.labels};

    % Safe extraction of types (handles empty type fields without dimension mismatch)
    types = cell(size(labels));
    for ch = 1:length(DATA(i_block).chanlocs)
        if isfield(DATA(i_block).chanlocs(ch), 'type') && ~isempty(DATA(i_block).chanlocs(ch).type)
            types{ch} = char(DATA(i_block).chanlocs(ch).type);
        else
            types{ch} = '';
        end
    end

    % Extract unique muscle base names (e.g. 'APB1' & 'APB2' -> 'APB')
    emg_mask   = strcmpi(types, 'EMG');
    emg_labels = labels(emg_mask);
    muscles    = unique(regexprep(emg_labels, '\d+$', ''), 'stable');

    if isempty(muscles)
        fprintf('  [EMG]  No EMG channels detected. Skipped.\n');
    elseif ~do_bipolar_emg
        fprintf('  [EMG]  Detected %d monopolar EMG channels (%d muscles: %s). Retained monopolar.\n', ...
            length(emg_labels), length(muscles), strjoin(muscles, ', '));
    else
        fprintf('  [EMG]  Deriving bipolar EMG for %d muscles:\n', length(muscles));
        for m = 1:length(muscles)
            m_name = muscles{m};
            idx1   = find(strcmpi(labels, [m_name, '1']));
            idx2   = find(strcmpi(labels, [m_name, '2']));

            if ~isempty(idx1) && ~isempty(idx2)
                % Compute bipolar signal (Lead 1 - Lead 2)
                DATA(i_block).data(idx1, :, :) = DATA(i_block).data(idx1, :, :) - DATA(i_block).data(idx2, :, :);
                DATA(i_block).chanlocs(idx1).labels = m_name;

                % Drop the redundant reference lead
                DATA(i_block).data(idx2, :, :) = [];
                DATA(i_block).chanlocs(idx2)   = [];

                fprintf('         Derived %s (%s1 - %s2). Removed %s2.\n', m_name, m_name, m_name, m_name);

                % Refresh labels after deletion for subsequent iterations
                labels = {DATA(i_block).chanlocs.labels};
            else
                fprintf('         Warning: Incomplete pair for %s (Lead 1: %d, Lead 2: %d). Skipped.\n', ...
                    m_name, ~isempty(idx1), ~isempty(idx2));
            end
        end
    end

    % ---------------------------------------------------------------------
    % 5. Finalise and Validate Structure
    % ---------------------------------------------------------------------
    DATA(i_block).nbchan = size(DATA(i_block).data, 1);
    DATA(i_block) = eeg_checkset(DATA(i_block));

    fprintf('  [INFO] Total remaining channels: %d/%d\n', DATA(i_block).nbchan, nbchan_old);
end

fprintf('\nDone!\n');

end