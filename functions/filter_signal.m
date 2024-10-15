function EEG = filter_signal(EEG,lp,hp,chansfilt,type)
% % Remove DC offsets and apply a high-pass filter (non-causal Butterworth impulse response function, 0.1 Hz half-amplitude cut-off, 12 dB/oct roll-off)
% EEG  = pop_basicfilter( EEG,  1:33 , 'Boundary', 'boundary', 'Cutoff',  0.1, 'Design', 'butter', 'Filter', 'highpass', 'Order',  2, 'RemoveDC', 'on' );
% % Apply a low-pass filter (non-causal Butterworth impulse response function, 20 Hz half-amplitude cut-off, 48 dB/oct roll-off) to the ERP waveforms
% ERP = pop_filterp( ERP,  1:35 , 'Cutoff',  20, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  8 );
% =========================================================================
% !!! IMPORTANT !!!
% Because filtfilt performs a forward-reverse filtering, it doubles the filter order!!!
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

% disp('==================================================================');

if strcmpi(type,'fieldtrip') % length(EEG)==1 % iscell(EEG)
    % Fieldtrip
    NTRL = length(EEG.trial);
    FNYQ = EEG.fsample/2;
    NCHN = size(EEG.trial{1},1);
    fprintf('Filtering ALL given signals across %d trials...\n',NTRL);

elseif  strcmpi(type,'eeglab')
    % EEGLAB
    NBLK = length(EEG);
    FNYQ = EEG(1).srate/2;

    % NCHN = size(EEG(1).data,1);
    % chaneeg = 1:NCHN;
    % This will break if FT data strct if given; will fix later
    % chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
    NCHN = length(chansfilt);

    fprintf('Filtering signals across %d blocks...\n',NBLK);
end

% Calculate filter coefficients
filterflag = false(2,1);
if ~isempty(lp)
    [bl, al] = butter(lp(2)./2,lp(1)/FNYQ,'low'); % Input: [ORDER, CUTOFF]
    assert(isstable(bl,al));
    filterflag(1) = true;

    fprintf('Using Butterworth lowpass [%1.2f Hz, order %d] filter.\n',lp(1),lp(2));
    if lp(1)<20, warning('Very low (<20 Hz) lowpass filter settings!'); end
end
if ~isempty(hp)
    [bh, ah] = butter(hp(2)./2,hp(1)/FNYQ,'high'); % Input: [ORDER, CUTOFF]
    assert(isstable(bh,ah));
    filterflag(2) = true;

    fprintf('Using Butterworth highpass [%1.2f Hz, order %d] filter.\n',hp(1),hp(2));
    if hp(1)>1.5, warning('Very high (>1.5 Hz) highpass filter settings!'); end
end
% [bs, as] = butter(1,[48 52]/FNYQ,'stop');

% Filtering
if  strcmpi(type,'fieldtrip')
    % Fieldtrip
    for i = 1:NTRL
        assert(size(EEG.trial{i},1)==NCHN);

        % Remove DC (big offsets can cause artifacts)
        EEG.trial{i} = remove_dcsignal(double(EEG.trial{i}), FNYQ);

        % Bandpass
        % First lowpass and then highpass
        if all(filterflag)
            EEG.trial{i} = filtfilt(bl,al, EEG.trial{i}');
            EEG.trial{i} = filtfilt(bh,ah, EEG.trial{i})';
        end
        % Lowpass
        if filterflag(1) && ~filterflag(2)
            EEG.trial{i} = filtfilt(bl,al, EEG.trial{i}')';
        end
        % Highpass
        if ~filterflag(1) && filterflag(2)
            EEG.trial{i} = filtfilt(bh,ah, EEG.trial{i}')';
        end

        assert(size(EEG.trial{i},1)==NCHN);
    end
else
    % EEGLAB
    for i = 1:NBLK
        assert(size(EEG(i).data(chansfilt,:),1)==NCHN);

        if ~isempty(EEG(1).event)
            eventLabels = {EEG(i).event(:).type};
            eventTimes  = [EEG(i).event(:).latency];
            mask = contains(eventLabels,'boundary');

            jumpsBad = floor([0 eventTimes(mask) size(EEG(i).data,2)]);
            jumpsBad = unique(jumpsBad);
        else
            jumpsBad = [0 size(EEG(i).data,2)];
        end

        NChunk = length(jumpsBad)-1;
        if NChunk>1
            fprintf('Block %d: Boundaries detected! Each chunk (N = %d) will be filtered separately.\n',i,NChunk);
        else
            fprintf('Block %d: Boundaries not detected.\n',i);
        end

        for j = 1:NChunk
            % Indices
            dataInd = [jumpsBad(j)+1 jumpsBad(j+1)];
            dataInd = dataInd(1):dataInd(2);

            if length(dataInd)>FNYQ*2
                dataTmp = double(EEG(i).data(chansfilt,dataInd));

                % Remove DC (big offsets can cause artifacts)
                dataTmp = remove_dcsignal(dataTmp, FNYQ);

                % Bandpass
                % First lowpass and then highpass
                if all(filterflag)
                    dataTmp = filtfilt(bl,al, dataTmp')';
                    dataTmp = filtfilt(bh,ah, dataTmp')';
                end
                % Lowpass
                if filterflag(1) && ~filterflag(2)
                    dataTmp = filtfilt(bl,al, dataTmp')';
                end
                % Highpass
                if ~filterflag(1) && filterflag(2)
                    dataTmp = filtfilt(bh,ah, dataTmp')';
                end
                % eeg(i).data = filtfilt(bs,as, eeg(i).data)';

                % Place back
                EEG(i).data(chansfilt,dataInd) = dataTmp;
            else
                warning('Chunk %d: Data too small (<=%d) for filtering (L = %d samples)!',j,FNYQ*2,length(dataInd));
            end

        end

        assert(size(EEG(i).data(chansfilt,:),1)==NCHN);
    end
end

fprintf('Done!\n');
% disp('==================================================================');

end