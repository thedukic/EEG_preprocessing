function EEG = reduce_linenoise(EEG)
%
% More or less just a wrapper function for the Zapline function
% Zapline is slightly modifed to accept p.fig1 / p.fig2 inputs
% SDukic, Feb 2024
%

thisMethod = 'Zaplineplus';
fprintf('\n%s 50Hz noise cleaning...\n',thisMethod);
NBLK = length(EEG);

switch thisMethod
    case 'Zapline'
        fh = figure('visible','off');
        th = tiledlayout(NBLK,2);
        th.TileSpacing = 'compact'; th.Padding = 'compact';

        % Normalised line frequnecy
        fline = 50/EEG(1).srate;

        % Remove line noise only from EEG channels (?)
        % eegchan = strcmp({EEG(1).chanlocs.type},'EEG');
        eegchan = 1:EEG(1).nbchan;

        if length(eegchan)==EEG(1).nbchan
            fprintf('Make sure that none of the channels are bipolar!\n');
        end

        % Function parameters
        p = [];
        p.nfft        = 4*EEG(1).srate; % default: 1024
        p.nkeep       = [];             % [] == use all PCs
        p.niterations = 1;

        % Remove line noise
        % Improve: Check automaticaly how many components should be removed
        for i = 1:NBLK
            fprintf('Cleaning block: %1d\n', i);

            p.fig1 = nexttile(2*i-1);
            p.fig2 = nexttile(2*i);

            % nt_zapline(EEG(i).data',fline);
            EEG(i).data(eegchan,:) = nt_zapline(EEG(i).data(eegchan,:)',fline,2,p,true)';
        end

        % Save
        plotX=25; plotY=25;
        set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
        set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
        print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval']),'-dtiff','-r400');
        close(fh);

    case 'Zaplineplus'
        % Is it better to merge block and do it in one go?
        % If using RS blocks, window size should be <= 120/8s = 15s
        for i = 1:NBLK
            fprintf('\nCleaning block: %1d\n', i);
            [EEG(i).data, ~, analyticsResults(i), fh] = clean_data_with_zapline_plus(EEG(i).data,EEG(i).srate,'noisefreqs',50,'winSizeCompleteSpectrum',20,'plotResults',true);

            % LineNoiseCleaning.ratioNoiseClean(i,:) = analyticsResults.ratioNoiseClean;
            % LineNoiseCleaning.proportionRemovedBelowNoise(i,:) = analyticsResults.proportionRemovedBelowNoise;

            % Save
            plotX=35; plotY=20;
            set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
            set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
            print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval_' num2str(i)]),'-dtiff','-r400');
            close(fh);
        end

        % Log / Report
        for i = 1:NBLK
            EEG(i).ALSUTRECHT.LineNoiseCleaning = analyticsResults;
        end
        % if all(LineNoiseCleaning.ratioNoiseClean(:,1)<1.2 & LineNoiseCleaning.ratioNoiseClean(:,1)>0.8)
        %     disp('50 Hz noise cleanining seems good!');
        % else
        %     disp('50 Hz noise cleanining might be suboptimal!');
        % end
end


end