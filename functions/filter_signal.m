function EEG = filter_signal(EEG,lp,hp,chansfilt,type)
% FILTER_SIGNAL: Filters EEG data (EEGLAB or FieldTrip format).
% Arguments:
%   EEG       - EEG data structure (EEGLAB/FieldTrip format).
%   lp        - Low-pass filter parameters: [Cutoff (Hz), Order].
%   hp        - High-pass filter parameters: [Cutoff (Hz), Order].
%   chansfilt - Channels to filter (indices).
%   type      - Data type: 'eeglab' or 'fieldtrip'.

% =========================================================================
% !!! IMPORTANT !!!
% Because filtfilt performs a forward-reverse filtering, it doubles the filter order!
% Filters should be applied to the continuous (rather than segmented) data.
%
% Order 1 ~ 6 dB/oct
% Order 1 becomes order 2  cos of two-pass and this is then 12 dB/oct
% Order 2 becomes order 4  cos of two-pass and this is then 24 dB/oct
% Order 4 becomes order 8  cos of two-pass and this is then 48 dB/oct
% Order 5 becomes order 10 cos of two-pass and this is then 60 dB/oct (?)
%
% Phase shift is frequency dependent and it is nonlinear in IIR.
% Increasing the order, makes the pahse shift more nonlinear
%
% Roisin SART: Highpass 0.3 Hz 5th order Butterworth
% https://socialsci.libretexts.org/Bookshelves/Psychology/Book%3A_Applied_Event-Related_Potential_Data_Analysis_(Luck)/14%3A_Appendix_3%3A_Example_Processing_Pipeline
% =========================================================================

% Check Inputs
% Validate 'lp' and 'hp' as either empty or a two-element vector of positive numbers
if ~isempty(lp)
    validateattributes(lp, {'numeric'}, {'numel', 2, 'positive'}, mfilename, 'lp');
else
    lp = []; % Ensure empty input is preserved
end
if ~isempty(hp)
    validateattributes(hp, {'numeric'}, {'numel', 2, 'positive'}, mfilename, 'hp');
else
    hp = [];
end
assert(ismember(lower(type), {'fieldtrip', 'eeglab'}), 'Unsupported data type.');

% Define variables
if strcmpi(type,'fieldtrip')
    NTRL = length(EEG.trial);
    FNYQ = EEG.fsample / 2;
    NCHN = size(EEG.trial{1},1);
    warning('Filtering ALL given signals across %d trials...\n',NTRL);

elseif  strcmpi(type,'eeglab')
    NBLK = length(EEG);
    FNYQ = EEG(1).srate / 2;
    NCHN = length(chansfilt);
    fprintf('Filtering signals across %d blocks...\n',NBLK);
end

% Calculate filter coefficients
typeFilter = [~isempty(lp), ~isempty(hp)];

% Check
if ~any(typeFilter)
    fprintf('Skipping the filtering as the parameters are not defined.\n'); return;
end

% Filter design
if typeFilter(1)
    [bLow, aLow] = butter(lp(2)/2, lp(1)/FNYQ, 'low');
    assert(isstable(bLow, aLow), 'Low-pass filter unstable.');
    
    fprintf('Using Butterworth lowpass [%1.2f Hz, order %d] filter.\n',lp(1),lp(2));
    if lp(1)<20, warning('N.B. Very low (<20 Hz) lowpass filter settings.'); end
end
if typeFilter(2)
    [bHigh, aHigh] = butter(hp(2)/2, hp(1)/FNYQ, 'high');
    assert(isstable(bHigh, aHigh), 'High-pass filter unstable.');
    
    fprintf('Using Butterworth highpass [%1.2f Hz, order %d] filter.\n',hp(1),hp(2));
    if hp(1)>1.5, warning('N.B. Very high (>1.5 Hz) highpass filter settings.'); end
end

% Filtering
if  strcmpi(type,'fieldtrip')
    for i = 1:NTRL
        assert(size(EEG.trial{i},1) == NCHN);

        % Remove DC (big offsets can cause artifacts)
        EEG.trial{i} = remove_dcsignal(double(EEG.trial{i}), FNYQ);

        % Bandpass
        % First lowpass and then highpass
        if all(typeFilter)
            EEG.trial{i} = filtfilt(bLow,aLow, EEG.trial{i}');
            EEG.trial{i} = filtfilt(bHigh,aHigh, EEG.trial{i})';
        end
        % Lowpass
        if typeFilter(1) && ~typeFilter(2)
            EEG.trial{i} = filtfilt(bLow,aLow, EEG.trial{i}')';
        end
        % Highpass
        if ~typeFilter(1) && typeFilter(2)
            EEG.trial{i} = filtfilt(bHigh,aHigh, EEG.trial{i}')';
        end

        assert(size(EEG.trial{i},1) == NCHN);
    end
else
    for i = 1:NBLK
        assert(size(EEG(i).data(chansfilt,:),1) == NCHN);

        if ~isempty(EEG(i).event)
            eventLabels = {EEG(i).event(:).type};
            eventTimes  = [EEG(i).event(:).latency];
            mask = contains(eventLabels,'boundary');
            boundaries = unique(floor([0 eventTimes(mask) size(EEG(i).data,2)]));
        else
            boundaries = [0 size(EEG(i).data,2)];
        end

        NChunk = length(boundaries)-1;
        if NChunk>1
            fprintf('Block %d: Boundaries detected! Each chunk (N = %d) will be filtered separately.\n',i,NChunk);
        else
            fprintf('Block %d: Boundaries not detected.\n',i);
        end

        for j = 1:NChunk
            % Indices
            % dataInd = [boundaries(j)+1 boundaries(j+1)];
            % dataInd = dataInd(1):dataInd(2);
            dataInd = (boundaries(j)+1):boundaries(j+1);

            if length(dataInd) > FNYQ*2
                % Remove DC (big offsets can cause artifacts)
                dataTmp = double(EEG(i).data(chansfilt,dataInd));
                dataTmp = remove_dcsignal(dataTmp, FNYQ);

                % Bandpass
                if all(typeFilter)
                    % First lowpass and then highpass
                    dataTmp = filtfilt(bLow,aLow, dataTmp')';
                    dataTmp = filtfilt(bHigh,aHigh, dataTmp')';
                end
                % Lowpass
                if typeFilter(1) && ~typeFilter(2)
                    dataTmp = filtfilt(bLow,aLow, dataTmp')';
                end
                % Highpass
                if ~typeFilter(1) && typeFilter(2)
                    dataTmp = filtfilt(bHigh,aHigh, dataTmp')';
                end

                % Place back
                EEG(i).data(chansfilt,dataInd) = dataTmp;
            else
                warning('Chunk %d: Too small (<= %d) for filtering (L = %d samples)!',j,FNYQ*2,length(dataInd));
            end

        end

        assert(size(EEG(i).data(chansfilt,:),1) == NCHN);
    end
end

fprintf('Done!\n');

end