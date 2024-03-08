function EEG = report_issues(EEG,thisTask)

if isfield(EEG,'data')
    % ALS number
    issues_to_check.aFileName = EEG.ALSUTRECHT.subject.id;
    issues_to_check.NumberTrials1 = sum([EEG.ALSUTRECHT.eventinfo{:,3}]);   % Total possible
    if strcmpi(thisTask,'RS')
        issues_to_check.NumberTrials2 = floor(length(EEG.times)/EEG.srate); % Left after preproc
    else
        issues_to_check.NumberTrials2 = size(EEG.data,3);
    end

    % Bad electrodes
    issues_to_check.FlatElectrodesDiscrepancy = EEG.ALSUTRECHT.badchaninfo.flatElectrodesDiscrepancy;
    if length(EEG.ALSUTRECHT.badchaninfo.badElectrodes)/128>0.2
        issues_to_check.RejectedTooManyElectrodes = length(EEG.ALSUTRECHT.badchaninfo.badElectrodes);
    else
        issues_to_check.RejectedTooManyElectrodes = 0;
    end

    % Extreme noise
    if EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier>0.2
        issues_to_check.HighProportionExcludedAsExtremeOutlier = EEG.ALSUTRECHT.extremeNoise.proportionExcludedForExtremeOutlier;
    else
        issues_to_check.HighProportionExcludedAsExtremeOutlier = 0;
    end
    if EEG.ALSUTRECHT.MWF.R1.ProportionOfDataShowingMuscleActivityTotal > 0.5
        issues_to_check.HighProportionOfEMG = EEG.ALSUTRECHT.MWF.R1.ProportionOfDataShowingMuscleActivityTotal;
    else
        issues_to_check.HighProportionOfEMG = 0;
    end

    % MWF
    MFWrounds = fields(EEG.ALSUTRECHT.MWF);
    for j = 1:length(MFWrounds)
        issues_to_check.(['MWF' MFWrounds{j} 'Status1']) =  EEG.ALSUTRECHT.MWF.(MFWrounds{j}).status;

        % If SER/ARR are NaN, then MWF failed
        if isnan(EEG.ALSUTRECHT.MWF.(MFWrounds{j}).signalToErrorRatio) || isnan(EEG.ALSUTRECHT.MWF.(MFWrounds{j}).artifactToResidueRatio)
            issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  0;
        else
            issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  1;
        end

        % Bad data should not be more than 50-ish%
        if EEG.ALSUTRECHT.MWF.(MFWrounds{j}).proportionMarkedForMWF>0.6
            issues_to_check.(['MWF' MFWrounds{j} 'BadData']) =  EEG.ALSUTRECHT.MWF.(MFWrounds{j}).proportionMarkedForMWF;
        else
            issues_to_check.(['MWF' MFWrounds{j} 'BadData']) =  0;
        end
    end
    % if ~isempty(tmpLabels)
    %     issues_to_check.HighProportionOfBadDataMWF = strjoin(aField(tmpLabels),', ');
    % else
    %     issues_to_check.HighProportionOfBadDataMWF = 0;
    % end

    % ICA/ICLabel/wICA
    if strcmpi(EEG.ALSUTRECHT.ica.DataLengthForValidICA,'OK')
        issues_to_check.DataTooShortForValidICA = 0;
    else
        issues_to_check.DataTooShortForValidICA = 1;
    end
    if EEG.ALSUTRECHT.ica.proportionArtifactICsReducedbywICA>0.3
        issues_to_check.HighProportionOfArtifactICs = EEG.ALSUTRECHT.ica.proportionArtifactICsReducedbywICA;
    else
        issues_to_check.HighProportionOfArtifactICs = 0;
    end

    % EC-RS eye blinks
    if strcmpi(thisTask,'RS')
        if EEG.ALSUTRECHT.blockinfo.ec_blinks>0
            issues_to_check.ECEyeBinksDetected = EEG.ALSUTRECHT.blockinfo.ec_blinks;
        else
            issues_to_check.ECEyeBinksDetected = 0;
        end
    end

    % Leftovers
    if EEG.ALSUTRECHT.leftovers.muscle>0.25
        issues_to_check.MuscleLeftovers = EEG.ALSUTRECHT.leftovers.muscle;
    else
        issues_to_check.MuscleLeftovers = 0;
    end
    if abs(mean(EEG.ALSUTRECHT.leftovers.blinks)-1)>0.05
        issues_to_check.EyeLeftovers = mean(EEG.ALSUTRECHT.leftovers.blinks)-1;
    else
        issues_to_check.EyeLeftovers    = 0;
    end

else
    % The participant was not processed
    issues_to_check.aFileName                              = EEG.ALSUTRECHT.subject.id;
    issues_to_check.NumberTrials1                          = NaN;
    issues_to_check.NumberTrials2                          = NaN;
    issues_to_check.FlatElectrodesDiscrepancy              = NaN;
    issues_to_check.RejectedTooManyElectrodes              = NaN;
    issues_to_check.HighProportionExcludedAsExtremeOutlier = NaN;
    issues_to_check.HighProportionOfEMG                    = NaN;

    % Note sure how to automate this part and not cheat
    MFWrounds = {'R1','R2','R3','R4'}; % Max 4 MWF rounds of cleaning
    for j = 1:length(MFWrounds)
            issues_to_check.(['MWF' MFWrounds{j} 'Status1']) =  NaN;
            issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  NaN;
            issues_to_check.(['MWF' MFWrounds{j} 'BadData']) =  NaN;
    end

    issues_to_check.DataTooShortForValidICA     = NaN;
    issues_to_check.HighProportionOfArtifactICs = NaN;

    if strcmpi(thisTask,'RS')
        issues_to_check.ECEyeBinksDetected = NaN;
    end

    issues_to_check.MuscleLeftovers = NaN;
    issues_to_check.EyeLeftovers    = NaN;
end

% Log
EEG.ALSUTRECHT.issues_to_check = issues_to_check;

end