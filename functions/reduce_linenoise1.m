function DATA = reduce_linenoise1(DATA, cfg)
% More or less just a wrapper function for the Zapline function
% Zapline is slightly modifed to accept p.fig1 / p.fig2 inputs

method = 'Zaplineplus';
freq_noise = 50;

fprintf('\n================================\n');
fprintf('Removing %d Hz noise (%s)\n', freq_noise, method);
fprintf('================================\n');

% ECG was not recorded at the beginning of the study
% but if VEOG/HEOG are ok, then ECG should also be OK
fprintf('Making sure that external channels are not bipolar...\n');
chanheog1 = contains({DATA(1).chanlocs.labels}, 'HEOG');
chanheog2 = contains({DATA(1).chanlocs.labels}, 'VEOG');
assert(sum(chanheog1) == 2 & sum(chanheog2) == 2);

num_blocks = length(DATA);

% Method
switch method
    case 'Zaplineplus'
        % Is it better to merge block and do it in one go?
        % If using RS blocks, window size should be <= 120/8s = 15s
        fprintf('Cleaning all blocks at once....\n');

        % Clean all blocks at once
        [data_clean, ~, analyticsResults, fh] = clean_data_with_zapline_plus(double(cat(2, DATA(:).data)), DATA(1).srate, ...
            'noisefreqs', freq_noise, ...
            'winSizeCompleteSpectrum', 20, ...
            'plotResults', true, ...
            'plotVisible', cfg.figure.visible);

        % Split back into blocks
        mask_block = make_blockmasks(DATA);
        mask_block = mask_block(1).ALSUTRECHT.blockinfo.rs_mask;
        assert(size(data_clean, 2) == size(mask_block,2));

        for i_block = 1:num_blocks
            DATA(i_block).data = data_clean(:, mask_block(i_block,:));
            assert(size(DATA(i_block).data, 2) == DATA(i_block).pnts);
        end

        % Check
        DATA = eeg_checkset(DATA);

        % Save
        save_figure(fh, DATA(1).ALSUTRECHT.subject.figures, [DATA(1).ALSUTRECHT.subject.id '_linenoise_1'], [35 20]);

        % Log
        for i_block = 1:num_blocks
            DATA(i_block).ALSUTRECHT.LineNoiseCleaning1 = analyticsResults;
        end

        % % Old code: does cleaninging per block
        % for i = 1:NBLK
        %     fprintf('\nCleaning block: %1d\n', i);
        %     [EEG(i).data, ~, analyticsResults(i), fh] = clean_data_with_zapline_plus(EEG(i).data,EEG(i).srate,'noisefreqs',50,'winSizeCompleteSpectrum',20,'plotResults',true);
        %
        %     % LineNoiseCleaning.ratioNoiseClean(i,:) = analyticsResults.ratioNoiseClean;
        %     % LineNoiseCleaning.proportionRemovedBelowNoise(i,:) = analyticsResults.proportionRemovedBelowNoise;
        %
        %     % % Save
        %     % plotX=35; plotY=20;
        %     % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
        %     % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
        %     % print(fh, fullfile(EEG(1).ALSUTRECHT.subject.figures, [EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval_' num2str(i)]),'-dtiff','-r400');
        %     % close(fh);
        % end

    case 'Zapline'
        % fh = figure('visible','off');
        % th = tiledlayout(NBLK,2);
        % th.TileSpacing = 'compact'; th.Padding = 'compact';
        %
        % % Normalised line frequnecy
        % freqLine = 50/EEG(1).srate;
        %
        % % Remove line noise only from EEG channels (?)
        % % eegchan = strcmp({EEG(1).chanlocs.type},'EEG');
        % eegchan = 1:EEG(1).nbchan;
        %
        % % Function parameters
        % params = [];
        % params.nfft        = 4*EEG(1).srate; % default: 1024
        % params.nkeep       = [];             % [] == use all PCs
        % params.niterations = 1;
        %
        % % Remove line noise
        % % Improve: Check automaticaly how many components should be removed
        % for i = 1:NBLK
        %     fprintf('Cleaning block: %1d\n', i);
        %
        %     params.fig1 = nexttile(2*i-1);
        %     params.fig2 = nexttile(2*i);
        %
        %     % nt_zapline(EEG(i).data',fline);
        %     EEG(i).data(eegchan,:) = nt_zapline(EEG(i).data(eegchan,:)',freqLine,2,params,true)';
        % end
        %
        % % Save
        % plotX=25; plotY=25;
        % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
        % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
        % print(fh, fullfile(EEG(1).ALSUTRECHT.subject.figures, [EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval']),'-dtiff','-r400');
        % close(fh);
end

end