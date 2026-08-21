function report_badsegments(EEG, mask_noise, tag_figure, opt_figure)
% Plot excluded chunks of EEG data
% EEG: EEGLAB structure with sampling rate and block length
% mask_noise: cell array with samples being excluded (each row represents a chunk)

NBLK = length(EEG);
assert(length(mask_noise) == NBLK);

% Create a figure
fh = figure('Visible', opt_figure);
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
    for j = 1:size(mask_noise{i},1)
        startSamp = mask_noise{i}(j,1);
        stopSamp  = mask_noise{i}(j,end);
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
save_figure(fh, EEG(1).ALSUTRECHT.subject.figures, [EEG(1).ALSUTRECHT.subject.id '_detected_' tag_figure], [15 NBLK*3]);

end