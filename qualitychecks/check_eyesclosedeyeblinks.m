function EEG = check_eyesclosedeyeblinks(EEG)
%
% SDukic, January 2025

fprintf('\n================================\n');
fprintf('Detecting eye blinks in the eyes closed resting-state recording\n');
fprintf('================================\n');

% Select these
chaneeg  = strcmp({EEG.chanlocs.type},'EEG'); % & contains({EEG.chanlocs.labels},'C');
chaneog  = strcmp({EEG.chanlocs.labels},'VEOG'); assert(sum(chaneog) == 1);
dataeeg  = EEG.data(chaneeg,:);
dataeog  = EEG.data(chaneog,:);
chanlocs = EEG.chanlocs(chaneeg);

% Only EC
ec_mask = EEG.ALSUTRECHT.blockinfo.rs_mask(~EEG.ALSUTRECHT.blockinfo.eo_mask,:);
ec_mask = any(ec_mask,1);
assert(length(ec_mask) == size(dataeeg,2));

dataeeg = dataeeg(:,ec_mask);
dataeogeo = dataeog(~ec_mask)';
dataeogec = dataeog(ec_mask)';
times     = EEG.times(ec_mask)/1000;
times     = times - times(1);

%% Strong EMG
% Temporarily filter EOG
[bl, al] = butter(2,15/(EEG.srate/2),'low');
[bh, ah] = butter(2,1/(EEG.srate/2),'high');
assert(isstable(bl,al), 'Lowpass filter unstable.');
assert(isstable(bh,ah), 'Highpass filter unstable.');

dataeogeo = filtfilt(bl,al,dataeogeo);
dataeogec = filtfilt(bl,al,dataeogec);
dataeogeo = filtfilt(bh,ah,dataeogeo);
dataeogec = filtfilt(bh,ah,dataeogec);

% Baseline correct
dataeogeo = dataeogeo - trimmean(dataeogeo,10);
dataeogec = dataeogec - trimmean(dataeogec,10);

% =========================================================================
% Fit t-location-scale distribution for data1
params1 = mle(dataeogeo, 'distribution', 'tlocationscale'); % Location, scale, df
fitted_pdf1 = @(x) tpdf((x - params1(1)) / params1(2), params1(3)) / params1(2);

% Fit t-location-scale distribution for data2
params2 = mle(dataeogec, 'distribution', 'tlocationscale');
fitted_pdf2 = @(x) tpdf((x - params2(1)) / params2(2), params2(3)) / params2(2);

% Goodness-of-fit testing
% Calculate log-likelihood
loglik1 = sum(log(fitted_pdf1(dataeogeo)));
loglik2 = sum(log(fitted_pdf2(dataeogec)));

% Number of parameters in t-location-scale (mu, sigma, nu)
num_params = 3;

% Compute BIC
bicEO = -2 * loglik1 + num_params * log(length(dataeogeo));
bicEC = -2 * loglik2 + num_params * log(length(dataeogec));

fprintf('BIC for EO VEOG histogram: %.4f\n', bicEO);
fprintf('BIC for EC VEOG histogram: %.4f\n', bicEC);

if bicEC < bicEO
    % -> better fit!
    flagGoodFit = true;
    fprintf('Great! The BIC for EC fit should be lower.\n');
else
    flagGoodFit = false;
    warning('Unexpected. The BIC for EC fit should be lower.');
end

% Plot histograms and fitted PDFs
fh = figure;
th = tiledlayout(2,2);
th.TileSpacing = 'compact'; th.Padding = 'compact';

nexttile; hold on;
hb = histogram(dataeogeo, 'Normalization', 'pdf');
hb.FaceColor = 0.7*ones(1,3);
myXlim = [min(dataeogeo) max(dataeogeo)];
fplot(fitted_pdf1, myXlim, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
title(['EO VEOG histogram, BIC = ' num2str(round(bicEO))]);
xlabel('Amplitude (uV)'); ylabel('Probability Density'); xlim(myXlim);

nexttile; hold on;
hb = histogram(dataeogec, 'Normalization', 'pdf');
hb.FaceColor = 0.7*ones(1,3);
myXlim = [min(dataeogec) max(dataeogec)];
fplot(fitted_pdf2, myXlim, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
title(['EC VEOG histogram, BIC = ' num2str(round(bicEC))]);
xlabel('Amplitude (uV)'); ylabel('Probability Density'); xlim(myXlim);

% =========================================================================
% Estimate the treshold using the EO RS data
% -> use a very high cutoff
treshold = prctile(dataeogeo,75) + 5*iqr(dataeogeo);
treshold = max(treshold,100);
fprintf('Eye blinks detection using a treshold of %1.0f uV.\n',treshold);

% Find the blinks in EC using that (EO) treshold
dataDetect = dataeogec;
[qrspeaks, locs] = findpeaks(dataDetect,times,'MinPeakHeight',treshold);

% % Detect eye blinks
% EEG_EO = EEG;
% EEG_EO.data = EEG_EO.data(:,~ec_mask);
%
% eventLabels = {EEG_EO.event(:).type};
% eventTimes  = [EEG_EO.event(:).latency];
% mask = contains(eventLabels,'boundary');
% mask = eventTimes(mask) < size(EEG_EO.data,2);
% EEG_EO.event = EEG_EO.event(mask);
%
% [eyeBlinksMask, eyeBlinksEpochs, BlinkMaxLatency, eyeBlinkData, treshold]= detect_veog(EEG_EO,500);

NEOG = length(locs);
if NEOG > 0
    % figure; hold on;
    % plot(times(1:100*256),dataeog(1:100*256));
    % plot(locs(1:11),qrspeaks(1:11),'ro');

    badEpoch2 = round(locs * EEG.srate);
    badEpoch2 = [badEpoch2-80; badEpoch2+80]';
    N = length(badEpoch2(1,1):badEpoch2(1,2));

    badEpoch2(badEpoch2<1) = 1;
    badEpoch2(badEpoch2>length(dataDetect)) = length(dataDetect);

    EOG = NaN(NEOG,N);
    for i = 1:NEOG
        if length(badEpoch2(i,1):badEpoch2(i,2)) == N
            EOG(i,:) = dataDetect(badEpoch2(i,1):badEpoch2(i,2));
        end
    end
    EOG(isnan(EOG(:,1)),:) = [];
    mEOG = median(EOG,1);

    nexttile; hold on;
    F = (0:size(EOG,2)-1)./EEG.srate;
    plot(F,EOG,'LineWidth',1.2);
    set(gca, 'ColorOrder',brewermap(size(EOG,1),'Greys'));
    plot(F,mEOG,'Color',[0.8 0.1 0.1],'LineWidth',3);
    title(['N = ' num2str(NEOG) ', treshold = ' num2str(round(treshold)) ' uV']);
    xlabel('Time (s)'); ylabel('EC VEOG amplitude');
    % pbaspect([1.618 1 1]);

    EEGtmp = NaN(NEOG,sum(chaneeg));
    for i = 1:NEOG
        if length(badEpoch2(i,1):badEpoch2(i,2)) == N
            EEGtmp(i,:) = mean(dataeeg(:,badEpoch2(i,1):badEpoch2(i,2)).^2,2);
        end
    end
    EEGtmp(isnan(EEGtmp(:,1)),:) = [];
    mEEGtmp = median(EEGtmp,1);

    % myCmap = brewermap(128,'*RdBu');
    myCmap = brewermap(128,'Greys');
    nexttile;
    topoplot(mEEGtmp,chanlocs,'maplimits',max(abs(mEEGtmp))*[0 1],'headrad',0.5,'colormap',myCmap,'whitebk','on','electrodes','off','style','map','shading','interp');
    title('EC EEG timelocked to eye blinks');
    hcb = colorbar;
    hcb.Title.String = "uV^2";

    % Report
    warning('The participant has eye blinks detected during EC.');

else
    fprintf('Great! No eye blinks found during eyes closed resting-state recording.\n');
    mEEGtmp = NaN(1,128);
end

% Save
plotX=25; plotY=25;
set(fh,'InvertHardCopy','Off','Color',[1 1 1]);
set(fh,'PaperPositionMode','Manual','PaperUnits','Centimeters','PaperPosition',[0 0 plotX plotY],'PaperSize',[plotX plotY]);
print(fh,fullfile(EEG.ALSUTRECHT.subject.preproc,[EEG.ALSUTRECHT.subject.id '_EC_eyeblinks']),'-dtiff','-r200');
close(fh);

% Log
EEG.ALSUTRECHT.blockinfo.ec_NBlinks     = NEOG;
EEG.ALSUTRECHT.blockinfo.ec_flagGoodFit = flagGoodFit;
EEG.ALSUTRECHT.blockinfo.ec_BICFit      = [bicEC bicEO];
EEG.ALSUTRECHT.blockinfo.ec_blinksTopo  = mEEGtmp;

fprintf(EEG.ALSUTRECHT.subject.fid,'\n---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Eyes closed resting-state eye blinks detection\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'---------------------------------------------------------\n');
fprintf(EEG.ALSUTRECHT.subject.fid,'Detected: %d\n', NEOG);

end