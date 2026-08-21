function DATA = do_dss_all(DATA, artifactType)

fprintf('\n================================\n');
fprintf('Removing %s artifact\n', artifactType);
fprintf('================================\n');

NDSS = 0.2;

% Load individual/template ICA templates
templatesICA = load_ictempalteweights(DATA(1));

% Make a block mask
block_mask = make_blockmask([DATA(:).pnts]);

NBLK = length(DATA);
DATAALL = [DATA(:).data];
fprintf('\nRemoving %s from %d electrodes from all %d blocks at once...\n', artifactType, DATA(1).nbchan, NBLK);

% Detect artifacts
switch artifactType
    case 'ecg'
        % % Detect QRS
        % R = 0.85;
        % Artifact_templates = templatesICA.Heartweights0;
        % Artifact_mask = templatesICA.Heartmask;
    case 'veog'
        % Detect VEOG
        R = 0.9;
        Artifact_templates = {templatesICA.Blinkweights0, templatesICA.Blinkweights1, templatesICA.Blinkweights2};
        Artifact_mask = templatesICA.Blinkmask;

    case 'heog'
        % Detect HEOG
        R = 0.9;
        Artifact_templates = {...
            templatesICA.Saccadeweights0, ...
            templatesICA.Saccadeweights1L, templatesICA.Saccadeweights1R, ...
            templatesICA.Saccadeweights2L, templatesICA.Saccadeweights2R ...
            };
        Artifact_mask = templatesICA.SaccademaskL | templatesICA.SaccademaskR;

end

assert(size(DATAALL,2) == length(Artifact_mask));

fh = figure;
tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");

% -------------------------------------------------------------------------
% Estimate data/noise covariances
C1 = nt_cov(bsxfun(@times, DATAALL', double(Artifact_mask')));
C0 = nt_cov(bsxfun(@times, DATAALL', double(~Artifact_mask')));

% Clean
[DATAALL_new, p1, indxRemove] = do_dss0(DATAALL, C0, C1, NDSS, Artifact_templates, R);
% disp(indxRemove);

% Double-check
A1_old = mean(bsxfun(@times, DATAALL', double(Artifact_mask')),1);
A1_new = mean(bsxfun(@times, DATAALL_new', double(Artifact_mask')),1);

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
ylabel('score'); xlabel('component'); title(['DSS: All ' num2str(NBLK) ' blocks']);

mytopoplot(A1_old(1:128),[],'Data cov: Before',nexttile);
mytopoplot(A1_new(1:128),[],'Data cov: After',nexttile);

% % Plot 2
% DATA_block_new.etc = rmfield(DATA_block_new.etc,'clean_sample_mask');
% DATA_block.etc = rmfield(DATA_block.etc,'clean_sample_mask');
% vis_artifacts(DATA_block_new, DATA_block);

% Return
if strcmpi(artifactType, 'ecg')
    % All electrodes (EEG + EXT + EMG) are cleaned, except ECG
    channel_mask = ~contains({DATA(1).chanlocs.labels}, 'ECG');
elseif ismember(artifactType, {'veog','heog'})
    % Only EEG electrodes are cleaned
    channel_mask = strcmp({DATA(1).chanlocs.type}, 'EEG');
end
DATA = store_back(DATA, DATAALL_new, channel_mask, block_mask);

% Save
plotX=20; plotY=10;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(DATA(1).ALSUTRECHT.subject.preproc, [DATA(1).ALSUTRECHT.subject.id '_removal_' artifactType]),'-dtiff','-r400');
close(fh);

end

% =========================================================================
% Helper fuctions
% =========================================================================

function [DATAclean, p1, indxRemove] = do_dss0(DATA, C0, C1, NDSS, TEMPLATE, R)
% Number of components to estimate
Nestim = round(NDSS * size(DATA,1));
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
    [r(:,i), p] = corr(fromdss(1:Ncheck,1:128)', TEMPLATE{i}(1:128));
end

r = abs(r);
indxRemove = find(any(r > R, 2));
% indxRemove = indxRemove(indxRemove < 6);

% Plot
figure; tiledlayout(1, Ncheck+1, "TileSpacing", "compact", "Padding", "compact");
for i_cmp = 1:Ncheck+1
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
    noise_components = todss' * DATA;
    DATAclean = nt_tsr(DATA', noise_components(indxRemove,:)')';
end
end

function block_mask = make_blockmask(Ns)
block_mask = [];
for i = 1:length(Ns)
    block_mask = [block_mask, i*ones(1, Ns(i))];
end
end

function DATA = store_back(DATA, DATAALL_new, channel_mask, block_mask)
for i_blk = 1:length(DATA)
    block_mask_tmp = block_mask == i_blk;
    assert(size(DATA(i_blk).data,2) == sum(block_mask_tmp));
    DATA(i_blk).data(channel_mask,:) = DATAALL_new(channel_mask, block_mask_tmp);
end
DATA = eeg_checkset(DATA);
end
