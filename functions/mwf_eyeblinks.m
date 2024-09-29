function [EEG, noiseMask] = mwf_eyeblinks(EEG)
%
% Test whether only eye blinks are present using bipolar VEOG
% SDukic, March 2023
%
fprintf('\nMWF (VEOG) eye blinks...\n');

% =========================================================================
% Detect eye blinks
[noiseMask, eyeBlinksEpochs, ~, dataeog] = detect_eog(EEG,400,true); % 400 ms
noiseMask = double(noiseMask);

% Account for very bad epochs affected across all channels
if any(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1)
    assert(length(noiseMask)==length(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1));
    noiseMask(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1) = NaN;
end

% =========================================================================
if ~isempty(eyeBlinksEpochs)

    NEOG = size(eyeBlinksEpochs,1);
    badIndx = find(EEG.ALSUTRECHT.extremeNoise.extremeNoiseEpochs1);
    T = linspace(-0.4,0.4,length(eyeBlinksEpochs(1,1):eyeBlinksEpochs(1,2)));
    cnt = 0;

    fh = figure; hold on;
    for i = 1:NEOG
        if ~any(ismember(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2),badIndx))
            y = dataeog(eyeBlinksEpochs(i,1):eyeBlinksEpochs(i,2));

            % absDiff_100ms   = abs(T - (-0.2));
            % minDiff_100ms   = min(absDiff_100ms(:));
            % [~, col_n100ms] = find(absDiff_100ms == minDiff_100ms);
            % absDiff_100ms   = abs(T - 0.2);
            % minDiff_100ms   = min(absDiff_100ms(:));
            % [~, col_p100ms] = find(absDiff_100ms == minDiff_100ms);
            % y = y - mean(y([1:col_n100ms, col_p100ms:end]));
            % % y = y - mean(y(1:col_m150ms));

            plot(T,y,'LineWidth',1.2);
        else
            cnt = cnt+1;
        end
    end
    NEOG = NEOG-cnt;

    set(gca,'ColorOrder',brewermap(NEOG,'BuGn'));
    title(['N = ' num2str(NEOG)]);
    xlabel('Time a.u.'); ylabel('VEOG amplitude');

    % Save
    plotX=35; plotY=20;
    set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
    set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
    print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_mwf_blinks']),'-dtiff','-r400');
    close(fh);

else
    noiseMask = zeros(1,EEG.pnts);
    NEOG = 0;
end

%% Multi-channel Wiener filter
% Noise mask:
% NaN - ignored samples (== very bad data)
% 0   - good samples
% 1   - bad samples     (== bad data that will be corrected)

% Double-check
assert(length(noiseMask)==size(EEG.data,2));

% Log info
EEG.ALSUTRECHT.MWF.R2.badElectrodes          = NaN;
EEG.ALSUTRECHT.MWF.R2.noiseMask              = noiseMask;
EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF = mean(noiseMask,'omitnan');

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'MWF (VEOG) eye blinks\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
% fprintf(EEG.ALSUTRECHT.subject.fid,'Eye blink amplitude threshold: %1.2f\n',treshold);
fprintf(EEG.ALSUTRECHT.subject.fid,'Number of eye blinks detected: %d\n',NEOG);
% fprintf(EEG.ALSUTRECHT.subject.fid,'Average eye blink duration for MWF: %1.2f\n',TEOG);
fprintf(EEG.ALSUTRECHT.subject.fid,'Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF);

fprintf('Bad data for MWF: %1.2f\n',EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF);

if EEG.ALSUTRECHT.MWF.R2.proportionMarkedForMWF>0.05
    % Parameters
    delayNew = round((10/1000)*EEG.srate);
    delaySpacingNew = round((16/1000)*EEG.srate);
    params  = mwf_params('delay',delayNew,'delay_spacing',delaySpacingNew);

    % MWF EEG only
    chaneeg = strcmp({EEG.chanlocs.type},'EEG');
    [cleanEEG, d, W, SER, ARR] = mwf_process(EEG.data(chaneeg,:),noiseMask,params);

    % Check if there were any problems
    if contains(lastwarn,"eigenvectors")
        warning('The MWF delay is too long?'); SER = Inf; ARR = Inf;
    end
    if isnan(SER) || isnan(ARR)
        warning('MWF did not fail but the MWF quality measures (SER/ARR) are NaN. The bad data might be too short.');
    end

    % % Visual inspection
    % EEG0 = EEG;
    % EEG0.data(chaneeg,:) = cleanEEG;
    % vis_artifacts(EEG0,EEG);

    % Return the clean data
    EEG.data(chaneeg,:) = cleanEEG;

    % Log
    % EEG.ALSUTRECHT.MWF.R2.estimatedArtifactInEachChannel = d;
    % EEG.ALSUTRECHT.MWF.R2.matrixUsedToEstimateArtifacts  = W;
    EEG.ALSUTRECHT.MWF.R2.status                 = 1;
    EEG.ALSUTRECHT.MWF.R2.signalToErrorRatio     = SER;
    EEG.ALSUTRECHT.MWF.R2.artifactToResidueRatio = ARR;

    fprintf('Signal-to-error ratio:     %1.2f\n',SER);
    fprintf('Artifact-to-residue ratio: %1.2f\n',ARR);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Signal-to-error ratio:     %1.2f\n',SER);
    fprintf(EEG.ALSUTRECHT.subject.fid,'Artifact-to-residue ratio: %1.2f\n',ARR);
else
    EEG.ALSUTRECHT.MWF.R2.status                 = 0;
    EEG.ALSUTRECHT.MWF.R2.signalToErrorRatio     = NaN;
    EEG.ALSUTRECHT.MWF.R2.artifactToResidueRatio = NaN;

    fprintf(EEG.ALSUTRECHT.subject.fid,'MWF is not done. Too little data.\n');
    fprintf('MWF is not done. Too little data.\n');
end

end