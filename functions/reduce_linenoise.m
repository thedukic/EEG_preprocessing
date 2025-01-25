function EEG = reduce_linenoise(EEG)
%
% More or less just a wrapper function for the Zapline function
% Zapline is slightly modifed to accept p.fig1 / p.fig2 inputs
% SDukic, May 2024
%

thisMethod = 'Zaplineplus';
fprintf('\n%s 50Hz noise cleaning...\n',thisMethod);
fprintf('Make sure that none of the channels are bipolar!\n');

NBLK = length(EEG);

switch thisMethod
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
        % print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval']),'-dtiff','-r400');
        % close(fh);

    case 'Zaplineplus'
        % Is it better to merge block and do it in one go?
        % If using RS blocks, window size should be <= 120/8s = 15s
        fprintf('Cleaning all blocks at once....\n');

        % Clean all blocks at once
        [dataeeg, ~, analyticsResults, fh] = clean_data_with_zapline_plus(cat(2,EEG(:).data),EEG(1).srate,'noisefreqs',50,'winSizeCompleteSpectrum',20,'plotResults',true);

        % Split back the blocks
        EEG2 = make_rsmasks(EEG);
        assert(size(dataeeg,2) == size(EEG2(1).ALSUTRECHT.blockinfo.rs_mask,2));

        for i = 1:NBLK
            EEG(i).data = dataeeg(:,EEG2(i).ALSUTRECHT.blockinfo.rs_mask(i,:));
            assert(size(EEG(i).data,2) == EEG(i).pnts);
        end
        EEG = eeg_checkset(EEG);

        % Save
        plotX=35; plotY=20;
        set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
        set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
        print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval_1']),'-dtiff','-r300');
        close(fh);

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
        %     % print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval_' num2str(i)]),'-dtiff','-r400');
        %     % close(fh);
        % end

        % Log / Report
        for i = 1:NBLK
            EEG(i).ALSUTRECHT.LineNoiseCleaning1 = analyticsResults;
        end
end


end