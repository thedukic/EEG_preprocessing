function EEG = reduce_linenoiseleftovers(EEG)
%
% TODO:
%   1. topoplots of pvalues masked using pvals
%   2. final spectra using a log-scale
%
% SDukic, October 2024

% Pvalue trashold
Ptrsh   = 1e-4; % P ~ 0.05/136
Lepoch0 = 4;

fprintf('\n================================\n');
fprintf('Checking if there are any 50 Hz noise leftovers\n');
fprintf('================================\n');

NBLK    = length(EEG);
chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
myCmap1 = brewermap(256,'BuPu');
myCmap2 = brewermap(length(chaneeg),'BrBG');

fh = figure;
th = tiledlayout(NBLK,3);
th.TileSpacing = 'compact'; th.Padding = 'compact';

badelec = cell(NBLK,1);
pval    = cell(NBLK,1);

for i_blk = 1:NBLK
    % Epoch into 1s (or maybe better into 1s with 0.5 overlap, but OK)
    Lepoch  = Lepoch0 * EEG(i_blk).srate;
    Nepoch  = floor(size(EEG(i_blk).data,2) / Lepoch);

    if Nepoch > 10
        dataeeg = reshape(EEG(i_blk).data(:,1:Nepoch*Lepoch),EEG(i_blk).nbchan,Lepoch,Nepoch);

        % Estimate power
        % Windowing function
        winfunc    = hann(Lepoch);
        label      = {EEG(i_blk).chanlocs.labels};
        psdspectra = NaN(Lepoch/2+1,EEG(i_blk).nbchan,Nepoch);
        for i_trl = 1:Nepoch
            [psdspectra(:,:,i_trl), freq] = pwelch(dataeeg(:,:,i_trl)', winfunc, 0, Lepoch, EEG(i_blk).srate);
        end

        % Log transform
        psdspectra = 10 * log10(psdspectra);

        % Averge across trials
        psdspectra2a = mean(psdspectra,3)'; % better to catch the extreme values
        % psdspectra2a = trimmean(psdspectra,10,'round',3)';

        % % Normalise
        % freqsel = true(size(freq));
        % freqsel(freq>49 & freq<51) = false;
        % psdspectra = trimmean(psdspectra,10,'round',3)';
        % psdspectra = psdspectra ./ sum(psdspectra(:,freqsel),2);

        % % Log transform
        % psdspectra2a = 10*log10(psdspectra2a);

        % % Plot
        % fh = figure;
        % plot(freq,psdspectra2,'LineWidth',1.2);
        % xlim([1 70]); ylim([0 0.15]);
        %
        % % Plot
        % fh = figure;
        % plot(freq,psdspectra2(129:end,:),'LineWidth',1.2);
        % xlim([1 70]); ylim([0 0.15]); legend(label(129:end));

        %%
        % Frequency masks
        if i_blk == 1
            freqsel50 = false(size(freq));
            % if length(lineBandwith) == 2
            %     freqsel50(freq>=lineBandwith(1) & freq<=lineBandwith(2)) = true;
            % elseif length(lineBandwith) == 1
            %     freqsel50(freq == lineBandwith) = true;
            % end
            freqsel50(freq == 50) = true;

            freqselrest = false(size(freq));
            % freqselrest(freq>=48 & freq<=49) = true;
            % freqselrest(freq>=51 & freq<=52) = true;
            freqselrest(freq == 49) = true;
            freqselrest(freq == 51) = true;

            % freq(freqsel50)
            % freq(freqselrest)
        end

        % X = mean(psdspectra2(:,freqsel50),2);
        % Y = mean(psdspectra2(:,freqselrest),2);
        %
        % fh = figure;
        % plot(X-Y,'LineWidth',1.2);

        X = squeeze(mean(psdspectra(freqsel50,:,:),1))';
        Y = squeeze(mean(psdspectra(freqselrest,:,:),1))';

        % Only test if 50 Hz is higher than the neighbouring freqs
        [h, pval{i_blk}] = ttest2(X,Y,'tail','right');
        % [h,pval{i}] = ttest(X,Y,'tail','right');

        [pvalsort,b] = sort(pval{i_blk},2,"descend");

        % Determine channels with 50 Hz
        badelec{i_blk} = find(pval{i_blk} < Ptrsh);

        % nexttile; plot(-log10(pval{i}),'LineWidth',1.2);
        % ylim([0 4]); xlim([1 EEG(i).nbchan]); ylabel('-log10(p)');
        % title(['Block ' num2str(i)]);

        nexttile((i_blk-1)*3+1);
        topoplot(-log10(pval{i_blk}(chaneeg)),EEG(i_blk).chanlocs,'maplimits',[0 -log(Ptrsh)],'headrad','rim','whitebk','on','style','map','electrodes','on','emarker2',{badelec{i_blk},'d','k',10,1},'shading','interp');
        title(['Block ' num2str(i_blk)]); axis tight; colormap(myCmap1);
        hcb = colorbar;
        hcb.Title.String = "-log_{10}(P)";

        if ~isempty(badelec{i_blk})
            NCHN = length(badelec{i_blk});
            fprintf('Block %d: Leftover line noise found in %d electrode(s). Fixing them now...\n',i_blk,NCHN);
            % disp(pval{i});

            % % Notch filter
            % [z1, p1] = butter(2,[47 53]./(EEG(i).srate/2),'stop');
            % for j = 1:NCHN
            %     % Pay attention to boundary events!
            %     EEG(i).data(badelec{i}(j),:) = filtfilt(z1,p1,EEG(i).data(badelec{i}(j),:));
            % end

            % Spectrum interpolation
            EEG(i_blk).data(badelec{i_blk},:) = my_dftfilter(EEG(i_blk).data(badelec{i_blk},:),EEG(i_blk).srate,50,'neighbour',1,2);

            % Recompute
            dataeeg = reshape(EEG(i_blk).data(:,1:Nepoch*Lepoch),EEG(i_blk).nbchan,Lepoch,Nepoch);
            psdspectra = NaN(Lepoch/2+1,EEG(i_blk).nbchan,Nepoch);
            for i_trl = 1:Nepoch
                [psdspectra(:,:,i_trl),freq] = pwelch(dataeeg(:,:,i_trl)',winfunc,0,Lepoch,EEG(i_blk).srate);
            end

            % Average and log-transform
            psdspectra2b = mean(psdspectra,3)';
            % psdspectra2b = trimmean(psdspectra,10,'round',3)';
            % psdspectra2b = psdspectra2b ./ sum(psdspectra2(b:,freqsel),2);
            psdspectra2b = 10*log10(psdspectra2b);

            % Mean P-value
            Pmean = mean(pval{i_blk}(badelec{i_blk}));
        else
            fprintf('Block %d: Nice! No leftover line noise is found.\n',i_blk);
            NCHN = 0;
            Pmean = NaN;
        end

        % Sort spectra
        psdspectra2a = psdspectra2a(b,:);

        % Plot 2
        nexttile((i_blk-1)*3+2); hold on;
        plot([50 50],[-100 50],'LineWidth',1.2,'Color',0.6*ones(1,3));
        plot([100 100],[-100 50],'LineWidth',1.2,'Color',0.6*ones(1,3));
        plot(freq,psdspectra2a','LineWidth',1.2); colororder(myCmap2);
        xlim([0 128]); ylim([-100 50]); pbaspect([1.618 1 1]);
        xlabel('Frequency (Hz)'); ylabel('10log_{10}(power)');

        % Plot 3
        nexttile((i_blk-1)*3+3); hold on;
        if ~isempty(badelec{i_blk})
            % Sort spectra
            psdspectra2b = psdspectra2b(b,:);
            % badelectmp   = b(badelec{i});
            badelectmp   = size(psdspectra2a,1):-1:size(psdspectra2a,1)-length(badelec{i_blk})+1;
            A = psdspectra2a(badelectmp,:)-mean(psdspectra2a(badelectmp,freqselrest),2);
            B = psdspectra2b(badelectmp,:)-mean(psdspectra2b(badelectmp,freqselrest),2);
            plot(freq,A,'LineWidth',1.2,'Color',0.7*ones(1,3));
            plot(freq,B,'LineWidth',1.2,'Color',0.3*ones(1,3));
        else
            A = psdspectra2a-mean(psdspectra2a(:,freqselrest),2);
            plot(freq,A,'LineWidth',1.2); colororder(myCmap2);
        end
        xlim([45 55]); ylim([-20 20]); pbaspect([1.618 1 1]);
        xlabel('Frequency (Hz)'); ylabel('10log_{10}(power)');
        if isnan(Pmean)
            title(['N = ' num2str(NCHN)]);
        else
            title(['N = ' num2str(NCHN) ', Pm = ' num2str(Pmean)]);
        end
        
    else
        warning('Block %d: Too little data (N = %ds) in this recoridng block.', i_blk, Nepoch*Lepoch0);
    end
end

% Save
if NBLK == 2
    plotX=20; plotY=10;
elseif NBLK==3
    plotX=20; plotY=15;
elseif NBLK==4
    plotX=20; plotY=20;
else
    plotX=25; plotY=30;
end
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_linenoiseremoval_2']),'-dtiff','-r300');
close(fh);

% Fix for reporting
pval{i_blk} = pval{i_blk}(badelec{i_blk});
for i_blk = 1:NBLK
    badelec{i_blk} = label(badelec{i_blk});
end

% Log / Report
for i_blk = 1:NBLK
    EEG(i_blk).ALSUTRECHT.LineNoiseCleaning2.badelec = badelec;
    EEG(i_blk).ALSUTRECHT.LineNoiseCleaning2.pval    = pval;
end

fprintf(EEG(1).ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG(1).ALSUTRECHT.subject.fid,'Electrodes with leftover 50 Hz noise fixed\n');
fprintf(EEG(1).ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');

badelec = [badelec{:}];
pval = [pval{:}];

if isempty(badelec)
    str = 'None';
else
    str = strjoin(badelec,', ');
end
fprintf(EEG(1).ALSUTRECHT.subject.fid,'Electrodes: %s\n', str);

if isempty(badelec)
    str = 'None';
else
    str = arrayfun(@(x) num2str(x,'%1.3f'),pval,'uni',0);
    str = strjoin(str,', ');
end
fprintf(EEG(1).ALSUTRECHT.subject.fid,'P-values:   %s\n', str);

end