if strcmpi(myPaths.task,'RS')
    for sub = 1:NSUB
        load(fullfile(myPaths.preproc,subjects{sub},[subjects{sub} '_' myPaths.visit '_' myPaths.task '_cleandata_' rnum 'b.mat']),'EEG');

        % Select only EEG
        chaneeg = strcmp({EEG.chanlocs.type},'EEG');
        dataeeg = EEG.data(chaneeg,:,:);

        % Compute power spectra
        [NCHN,NPTS,NTRL]= size(dataeeg);
        psdspectra = NaN(floor(NPTS/2+1),NCHN,NTRL);

        for trial = 1:NTRL
            [psdspectra(:,:,trial), freq] = pwelch(dataeeg(:,:,trial)',NPTS,0,NPTS,EEG.srate);
        end

        % Average the spectra
        psdspectra = mean(psdspectra,3);

        % Calculate measure of spread
        SelectFreqBins = freq >= 40 & freq <= 70;
        psdMean = mean(psdspectra(SelectFreqBins,:),1);
        psdSpread = std(psdMean);
        cutoff = 0.04; % adjust if needed
        
        % Plot and output if it exceeds the cut-off:
        if psdSpread > cutoff
            fprintf('Subject: %s - The spread of power across electrodes (%.2f) exceeds the cutoff value of %.2f.\n', subjects{sub}, psdSpread, cutoff);

            dataCmap = brewermap(sum(chaneeg),'BrBG');

            fh = figure;
            th = tiledlayout(1,2);
            th.TileSpacing = 'compact'; th.Padding = 'compact';

            th = nexttile;
            hold on; box off;
            plot(freq,psdspectra,'LineWidth',1.1);

            colororder(th,dataCmap); xlim([freq(1), 40]);
            title([subjects{sub},', group: ',myPaths.group, ', visit: ', myPaths.visit, ', task: ', myPaths.task,  ' - Power Spread: ', num2str(psdSpread, '%.2f')]);
            pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('Power (a.u.)');

            th = nexttile;
            hold on; box off;

            dataTmp = log10(psdspectra);
            dataClim = 1.1*[min(dataTmp(:)), max(dataTmp(:))];
            dataClim(1) = floor(dataClim(1));
            dataClim(2) = ceil(dataClim(2));
            plot(freq,dataTmp,'LineWidth',1.1);

            colororder(th,dataCmap); xlim([freq(1), 70]);
            title([subjects{sub},', group: ',myPaths.group, ', visit: ', myPaths.visit, ', task: ', myPaths.task, ' - Power Spread: ', num2str(psdSpread, '%.2f')]);
            pbaspect([1.618 1 1]); xlabel('Frequency (Hz)'); ylabel('log_{10}(Power) (a.u.)');

            % plotX=20; plotY=8;
            % set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
            % set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
            % print(fh,fullfile(subject.preproc,[subjects{sub} '_pspectra_final']),'-dtiff','-r300');
            % close(fh);
        end
    end
end