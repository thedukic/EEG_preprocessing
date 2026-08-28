function [psdspectra, freq, chaneeg, chanemg] = estimate_power(EEG, this_script)
% ESTIMATE_POWER Helper function for estimating power spectra.
%
% Syntax:
%   [psdspectra, freq, chaneeg, chanemg] = estimate_power(EEG, thisScript)
%
% Inputs:
%   EEG        - EEGLAB data structure or struct array
%   thisScript - Mode selector: 'preproc2', 'freport', or 'speaks'
%
% Outputs:
%   psdspectra - Power spectral density matrix (frequencies x channels)
%   freq       - Frequency vector corresponding to PSD rows
%   chaneeg    - Logical index vector for EEG channels
%   chanemg    - Logical index vector for EMG channels

% Channel type masks (case-insensitive)
if isfield(EEG(1).chanlocs, 'type')
    chaneeg = strcmpi({EEG(1).chanlocs.type}, 'EEG');
    chanemg = strcmpi({EEG(1).chanlocs.type}, 'EMG');
else
    chaneeg = true(1, EEG(1).nbchan);
    chanemg = false(1, EEG(1).nbchan);
end

% Fallback if no channels are explicitly labelled 'EEG'
if ~any(chaneeg)
    chaneeg = true(1, size(EEG(1).data, 1));
end

fs = EEG(1).srate;

switch lower(this_script)
    case 'preproc2'
        % @preproc_cleaning2 (@generate_finalplots)
        assert(ndims(EEG.data) == 3, 'Input EEG.data must be 3D for preproc2 mode.');

        is_MT = isfield(EEG, 'ALSUTRECHT') && ...
            isfield(EEG.ALSUTRECHT, 'subject') && ...
            isfield(EEG.ALSUTRECHT.subject, 'task') && ...
            strcmpi(EEG.ALSUTRECHT.subject.task, 'MT');

        if is_MT
            dataeeg = EEG.data;
        else
            dataeeg = EEG.data(chaneeg, :, :);
        end

    case 'speaks'
        % @preproc_cleaning1 (@reduce_spectrapeaks)
        winSizeCompleteSpectrum = 10; % [s]

        % Concatenate across struct array blocks if multiple files provided
        data_all = cat(2, EEG(:).data);
        data_sub = data_all(chaneeg, :);
        [nchn, npts_total] = size(data_sub);

        % Guard against zero-division for short recordings
        min_segments = 8;
        if winSizeCompleteSpectrum * fs > npts_total / min_segments
            winSizeCompleteSpectrum = max(1, floor(npts_total / min_segments / fs));
            warning('Dataset is short. Adjusted window size to %d s.', winSizeCompleteSpectrum);
        end

        npts = winSizeCompleteSpectrum * fs;
        ntrl = floor(npts_total / npts);

        if ntrl < 1
            error('Dataset too short (%d samples) for spectral calculation.', npts_total);
        end
        dataeeg = reshape(data_sub(:, 1:ntrl * npts), nchn, npts, ntrl);

        % case 'freport'
        %     % @preproc_cleaning1 (@leftover report) and @report_final
        %     if ndims(EEG.data) == 3
        %         dataeeg = EEG.data(chaneeg, :, :);
        %     else
        %         npts = 2 * fs; % 2-second segments
        %         data_sub = EEG.data(chaneeg, :);
        %         [nchn, npts_total] = size(data_sub);
        %         ntrl = floor(npts_total / npts);
        %
        %         if ntrl < 1
        %             error('Dataset too short (< %d samples) to construct 2s epochs.', npts);
        %         end
        %         dataeeg = reshape(data_sub(:, 1:ntrl * npts), nchn, npts, ntrl);
        %     end
        %
    otherwise
        error('Unrecognised mode: ''%s''. Valid options: ''preproc2'', ''freport'', ''speaks''.', this_script);
end

% Compute power spectra
[nchn, npts, ntrl] = size(dataeeg);
assert(ntrl >= 1, 'No valid trials available for PSD computation.');

% Baseline-correct each epoch and orient to [samples x channels x trials]
dataeeg = double(dataeeg);
dataeeg = dataeeg - mean(dataeeg, 2);
dataeeg = permute(dataeeg, [2 1 3]);

psdspectra = NaN(floor(npts / 2 + 1), nchn, ntrl);

for i = 1:ntrl
    [psdspectra(:, :, i), freq] = pwelch(dataeeg(:, :, i), npts, 0, npts, fs);
end

% Average across epochs/trials
psdspectra = mean(psdspectra, 3);

end