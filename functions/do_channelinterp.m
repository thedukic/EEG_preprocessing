function EEG = do_channelinterp(EEG, thisType)

if nargin < 2 || isempty(thisType)
    thisType = 'spherical';
end

fprintf('\n================================\n');
fprintf('Interpolating EEG electrodes (%s)\n', thisType);
fprintf('================================\n');

if EEG.nbchan < 128
    chanlocs_target = readlocs('biosemi128_eeglab.ced');
    chanlocs_old    = EEG.chanlocs;

    % Check if ICA mixing matrix exists and needs interpolation
    has_icawinv = isfield(EEG, 'icawinv') && ~isempty(EEG.icawinv);
    do_interp_ica = has_icawinv && (size(EEG.icawinv, 1) ~= 128);

    if do_interp_ica
        fprintf('Interpolating ICA mixing matrix from %d to 128 channels...\n', size(EEG.icawinv, 1));

        % Copy struct to retain chaninfo, srate, etc.
        EEG_temp             = EEG;
        EEG_temp.data        = double(EEG.icawinv); % [Channels x ICs]
        EEG_temp.chanlocs    = chanlocs_old;
        EEG_temp.nbchan      = length(chanlocs_old);
        EEG_temp.pnts        = size(EEG.icawinv, 2);
        EEG_temp.trials      = 1;

        % Explicitly wipe ICA fields to prevent eeg_interp from reindexing icachansind
        EEG_temp.icaweights  = [];
        EEG_temp.icasphere   = [];
        EEG_temp.icawinv     = [];
        EEG_temp.icachansind = [];
        EEG_temp.icaact      = [];

        % Interpolate missing channel rows across all component topographies
        EEG_temp = pop_interp(EEG_temp, chanlocs_target, thisType);

        % Store interpolated mixing matrix
        new_icawinv = EEG_temp.data;
    end

    % 1. Interpolate continuous scalp data
    EEG = pop_interp(EEG, chanlocs_target, thisType);

    % 2. Assign interpolated ICA mixing weights back to main EEG struct
    if do_interp_ica
        EEG.icawinv       = new_icawinv;
        EEG.icachansind   = 1:EEG.nbchan;
        EEG.icasphere     = eye(EEG.nbchan);
        EEG.icaweights    = pinv(EEG.icawinv);
        EEG.icaact        = [];

        if isfield(EEG, 'ALSUTRECHT') && isfield(EEG.ALSUTRECHT, 'ica')
            EEG.ALSUTRECHT.ica.icawinv     = EEG.icawinv;
            EEG.ALSUTRECHT.ica.icaweights  = EEG.icaweights;
            EEG.ALSUTRECHT.ica.icasphere   = EEG.icasphere;
            EEG.ALSUTRECHT.ica.icachansind = EEG.icachansind;
        end
    end
else
    fprintf('Skipping! All electrodes are already in the data.\n');
end

end