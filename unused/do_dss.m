function DATA = do_dss(DATA, artifactType)

fprintf('\n================================\n');
fprintf('Removing %s artifact\n', artifactType);
fprintf('================================\n');

NDSS = 0.5;
NBLK = length(DATA);

% Load individual/template ICA templates
templatesICA = load_ictempalteweights(DATA(1));

% Make a block mask
block_mask = make_blockmask([DATA(:).pnts]);

fh = figure;
tiledlayout(NBLK, 3, "TileSpacing", "compact", "Padding", "compact");

for i_blk = 1:NBLK
    DATA_block = DATA(i_blk);
    fprintf('\nBlock %d: Removing %s from %d electrodes at once...\n', i_blk, artifactType, DATA_block.nbchan);

    % -------------------------------------------------------------------------
    % % Make external electrodes bipolar
    % [DATA_tmp, isRecordedECG] = make_bipolar0(DATA_block);

    % Detect artifacts
    switch artifactType
        case 'ecg'
            % Detect QRS
            R = 0.80;
            Artifact_templates = {templatesICA.Heartweights0, templatesICA.Heartweights1, templatesICA.Heartweights2};
            Artifact_mask       = templatesICA.Heartmask;
            % Artifact_epochs   = templatesICA.Heartepochs;

            % winHeart = 30;
            % [Artifact_mask, Artifact_epochs, Artifact_latency, Artifact_data, pulsEstimate] = detect_ecg(DATA_tmp, winHeart, isRecordedECG);

        case 'veog'
            % Detect VEOG
            R = 0.80;
            Artifact_templates = {templatesICA.Blinkweights0, templatesICA.Blinkweights1, templatesICA.Blinkweights2};
            Artifact_mask      = templatesICA.Blinkmask;

            % winBlink = 150;
            % [Artifact_mask, Artifact_epochs, Artifact_maxLatency, Artifact_data, Artifact_treshold] = detect_veog(DATA_tmp, winBlink);

        case 'heog'
            % Detect HEOG
            R = 0.80;
            Artifact_templates = templatesICA.Saccadeweights0;
            Artifact_mask      = templatesICA.SaccademaskL | templatesICA.SaccademaskR;
            % Artifact_epochs   = [templatesICA.Saccadeepochs1; templatesICA.Saccadeepochs2];

            % winSaccade = 150;
            % [Artifact_mask, Artifact_epochs, Artifact_maxLatency, Artifact_data, Artifact_treshold] = detect_heog(DATA_tmp, winSaccade);

            % Artifact_epochs = [Artifact_epochs{1}; Artifact_epochs{2}];
            % Artifact_mask   = Artifact_mask{1} | Artifact_mask{2};
    end

    Artifact_mask = Artifact_mask(block_mask == i_blk);
    assert(size(DATA_block.data,2) == length(Artifact_mask));

    % -------------------------------------------------------------------------
    % Estimate data/noise covariances
    % [C1, A1] = estimate_covariance(DATA_block, Artifact_mask);
    C1 = nt_cov(bsxfun(@times, DATA_block.data', double(Artifact_mask')));
    C0 = nt_cov(bsxfun(@times, DATA_block.data', double(~Artifact_mask')));

    % Clean
    [DATA_block_new, p1, indxRemove] = do_dss0(DATA_block, C0, C1, NDSS, Artifact_templates, R);
    % disp(indxRemove);

    % Double-check
    % [C1, A1old] = estimate_covariance(DATA_block, Artifact_epochs);
    % [C1, A1new] = estimate_covariance(DATA_block_new, Artifact_epochs);
    A1old = mean(bsxfun(@times, DATA_block.data', double(Artifact_mask')),1);
    A1new = mean(bsxfun(@times, DATA_block_new.data', double(Artifact_mask')),1);

    % Plot 1
    % nexttile; plot(p1, '.-'); ylabel('score'); xlabel('component'); title (['DSS: Block ' num2str(i_blk)]); pbaspect([1.618 1 1]);
    components = 1:length(p1);
    figure(fh); nexttile;

    % Plot all component scores (the original line)
    plot(components, p1, 'k.-', 'MarkerSize', 10, 'LineWidth', 1);
    hold on;

    % Mark the components to be removed
    % We index both the components (X-axis) and p1 (Y-axis) using indxRemove.
    plot(components(indxRemove), p1(indxRemove), 'ro', ...
        'MarkerSize', 10, ... % Make the marker large and visible
        'LineWidth', 2, ...   % Make the outline stand out
        'DisplayName', 'Removed Components');

    % Add labels and aesthetics
    hold off; pbaspect([1.618 1 1]);
    ylabel('score'); xlabel('component'); title(['DSS: Block ' num2str(i_blk)]);

    mytopoplot(A1old(1:128),[],'Data cov: Before',nexttile); axis tight;
    mytopoplot(A1new(1:128),[],'Data cov: After',nexttile); axis tight;

    % % Plot 2
    % DATA_block_new.etc = rmfield(DATA_block_new.etc,'clean_sample_mask');
    % DATA_block.etc = rmfield(DATA_block.etc,'clean_sample_mask');
    % vis_artifacts(DATA_block_new, DATA_block);

    % Return
    if strcmpi(artifactType, 'ecg')
        % All electrodes (EEG + EXT + EMG) are cleaned, except ECG
        channel_mask1 = ~contains({DATA(i_blk).chanlocs.labels},'ECG');
        channel_mask2 = ~contains({DATA_block_new.chanlocs.labels},'ECG');
    elseif strcmpi(artifactType, 'veog')
        % Only EEG electrodes are cleaned
        channel_mask1 = strcmp({DATA(i_blk).chanlocs.type},'EEG');
        channel_mask2 = strcmp({DATA_block_new.chanlocs.type},'EEG');
    end
    DATA(i_blk).data(channel_mask1,:) = DATA_block_new.data(channel_mask2,:);
end

% Save
if NBLK == 2
    plotX=20; plotY=10;
elseif NBLK == 3
    plotX=18; plotY=15;
elseif NBLK == 4
    plotX=18; plotY=20;
else
    plotX=20; plotY=30;
end

set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(DATA(1).ALSUTRECHT.subject.preproc, [DATA(1).ALSUTRECHT.subject.id '_removal_' artifactType]),'-dtiff','-r400');
close(fh);

end

% =========================================================================
% Helper fuctions
% =========================================================================
% function [C1, A1] = estimate_covariance(DATA, epochs)
% N = size(epochs,1);
% C1 = NaN(DATA.nbchan,DATA.nbchan,N);
% A1 = NaN(DATA.nbchan,N);
%
% for i = 1:N
%     y = DATA.data(:, epochs(i,1):epochs(i,2));
%     % y = y - mean(y,1);
%     % y = y - mean(y,2);
%     C1(:,:,i) = y * y';
%     A1(:,i) = mean(y,2);
% end
%
% % [u, s, v] = svd(A1(1:128,:),'econ');
% % figure;
% % for i = 1:10
% % mytopoplot(u(:,i), [], '', nexttile);
% % end
%
% C1 = mean(C1,3);
% A1 = mean(A1,2);
%
% end

function block_mask = make_blockmask(Ns)
block_mask = [];
for i = 1:length(Ns)
    block_mask = [block_mask, i*ones(1, Ns(i))];
end
end

function [DATAclean, p1, indxRemove] = do_dss0(DATA, C0, C1, NDSS, TEMPLATE, R)
% Number of components to estimate
Nestim = round(NDSS * DATA.nbchan);
fprintf('DSS components estimated: %d\n', Nestim);

% DSS
[todss, pwr0, pwr1] = nt_dss0(C0,C1,Nestim,[]);

% Explained
p1 = pwr1 ./ pwr0;
p1 = p1 / sum(p1);

% nexttile; plot(p1, '.-');
% ylabel('score'); xlabel('component'); title ('DSS'); pbaspect([1.618 1 1]);

% DSS components
fromdss = pinv(todss);

% Detect ECG components
Ntemplate = length(TEMPLATE);
Ncheck = 10;
r = NaN(Ncheck,Ntemplate);
for i = 1:Ntemplate
    if ~isnan(TEMPLATE{i})
        [r(:,i), p] = corr(fromdss(1:Ncheck,1:128)', TEMPLATE{i}(1:128));
    end
end

% Remove NaNs
r = r(:,~isnan(r(1,:)));

r = abs(r);
indxRemove = find(any(r > R, 2));
% indxRemove = indxRemove(indxRemove < 6);

% Plot
figure; tiledlayout(1, Ncheck+1, "TileSpacing", "compact", "Padding", "compact");
for i_cmp = 1:Ncheck
    if any(indxRemove == i_cmp-1)
        strRempve = '*';
    else
        strRempve = '';
    end
    if i_cmp == 1
        mytopoplot(TEMPLATE{1}(1:128), [], 'Data covariance', nexttile);
    else
        mytopoplot(fromdss(i_cmp-1,1:128), [], [strRempve 'S = ' num2str(round(p1(i_cmp-1),3)) ', R = ' num2str(round(mean(r(i_cmp-1,:)),3)) strRempve], nexttile);
    end
end

if isempty(indxRemove)
    fprintf('DSS cleaning will be skipped: max{R} = %1.2f\n', max(r));
    DATAclean = DATA;
else
    fprintf('DSS cleaning: max{R} = %1.2f\n', max(r));
    fprintf('DSS components removed: %d\n', length(indxRemove));
    % disp(indxRemove);

    % Remove
    noise_components = todss' * DATA.data;
    cleandata = nt_tsr(DATA.data', noise_components(indxRemove,:)')';

    % Return
    DATAclean = DATA;
    DATAclean.data = cleandata;
end

end