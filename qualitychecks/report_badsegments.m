function report_badsegments(EEG,maskNoise,figTag)
% Plot excluded chunks of EEG data
% EEG: EEGLAB structure with sampling rate and block length
% maskNoise: cell array with samples being excluded (each row represents a chunk)

NBLK = length(EEG);
assert(length(maskNoise) == NBLK);

% Create a figure
fh = figure;
th = tiledlayout(NBLK,1);
th.TileSpacing = 'compact'; th.Padding = 'compact';

for i = 1:NBLK
    % Time axis in minutes
    timeFactor = EEG(i).srate * 60;
    timeVector = (0:EEG(i).pnts - 1) / timeFactor;

    nexttile; hold on;
    set(gca, 'YLim', [0 1], 'YTick', []);
    rectangle('Position', [0 0 timeVector(end) 1], 'FaceColor', [0.6 0.8 1], 'EdgeColor', 'none');

    % % Identify the excluded chunks
    % excludedStarts = find(diff([false; excludedMask]) == 1); % Start indices of excluded chunks
    % excludedStops = find(diff([excludedMask; false]) == -1); % Stop indices of excluded chunks

    % Plot the excluded regions with red stripes
    for j = 1:size(maskNoise{i},1)
        startSamp = maskNoise{i}(j,1);
        stopSamp  = maskNoise{i}(j,end);
        xPos = [startSamp, stopSamp] / timeFactor;
        rectangle('Position', [xPos(1), 0, diff(xPos), 1], 'FaceColor', [1 0 0 0.5], 'EdgeColor', 'none');
    end

    % Final touches
    if i == 1
        title('Excluded chunks of EEG data');
    end
    if i == NBLK
        xlabel('Time (min)');
    end
    ylabel(['Block ' num2str(i)]);
    grid on; hold off; axis tight;
end

% Save
plotX=20; plotY=NBLK*3;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG(1).ALSUTRECHT.subject.preproc,[EEG(1).ALSUTRECHT.subject.id '_badchunks_' figTag]),'-dtiff','-r300');
close(fh);

end