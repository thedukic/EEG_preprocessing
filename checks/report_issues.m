function qa_data = report_issues(qa_data)

% ALS number
issues_to_check.aFileName     = qa_data.subject.id;
issues_to_check.NumberTrials1 = sum([qa_data.eventinfo{:, 3}]);           % Total possible
issues_to_check.NumberTrials2 = qa_data.epochRejections.initialEpochs;    % Left after preproc1
issues_to_check.NumberTrials3 = qa_data.epochRejections.remainingEpochs;  % Left after preproc2

issues_to_check.NumberIC2 = max(qa_data.ica.num_req);                     % ICA (request)
issues_to_check.NumberIC3 = qa_data.ica.num_done;                         % ICA (done)
issues_to_check.NumberIC4 = qa_data.ica.num_max;                          % ICA ()
% issues_to_check.DataTooShortForValidICA = 0;

% % Number of blocks
% issues_to_check.NumberOfBlocks = length(qa_data.subject.datablocks);

% CMS dropouts
issues_to_check.CMSDropout = sum(qa_data.cmsDropouts.NDropouts);

% Flat channelss
% issues_to_check.FlatElectrodesDiscrepancy = qa_data.badchaninfo.flatElectrodesDiscrepancy;

% % Bad electrodes
% if length(qa_data.badchaninfo.badElectrodes)/128 > 0.2
%     issues_to_check.RejectedTooManyElectrodes = length(qa_data.badchaninfo.badElectrodes);
% else
%     issues_to_check.RejectedTooManyElectrodes = 0;
% end

% Extreme noise that is exlcluded from the analysis
% if qa_data.extremeNoise.proportionExcludedForExtremeOutlier > 0.2
%     issues_to_check.HighProportionExcludedAsExtremeOutlier = qa_data.extremeNoise.proportionExcludedForExtremeOutlier;
% else
%     issues_to_check.HighProportionExcludedAsExtremeOutlier = 0;
% end
% issues_to_check.HighProportionExcludedAsExtremeOutlier = 0;

% if isfield(qa_data, 'MWF')
%     % EMG noise before cleaning
%     if qa_data.MWF.R1.ProportionOfDataShowingMuscleActivityTotal > 0.5
%         issues_to_check.HighProportionOfEMGInitial = qa_data.MWF.R1.ProportionOfDataShowingMuscleActivityTotal;
%     else
%         issues_to_check.HighProportionOfEMGInitial = 0;
%     end
%
%     % MWF
%     MFWrounds = fields(qa_data.MWF);
%     MFWrounds = MFWrounds(contains(MFWrounds,'R'));
%     for j = 1:length(MFWrounds)
%         issues_to_check.(['MWF' MFWrounds{j} 'Status1']) = qa_data.MWF.(MFWrounds{j}).status;
%
%         % If SER/ARR are NaN, then MWF failed
%         if isnan(qa_data.MWF.(MFWrounds{j}).signalToErrorRatio)
%             issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  NaN;
%         elseif isinf(qa_data.MWF.(MFWrounds{j}).signalToErrorRatio)
%             issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  Inf;
%         else
%             issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  1;
%         end
%     end
% end

% % ICA / ICLabel
% if strcmpi(qa_data.ica.DataLengthForValidICA, 'OK')
%     issues_to_check.DataTooShortForValidICA = 0;
% else
%     issues_to_check.DataTooShortForValidICA = 1;
% end

% if qa_data.ica.proportionArtifactICs > 0.3
%     issues_to_check.HighProportionOfArtifactICs = qa_data.ica.proportionArtifactICs;
% else
%     issues_to_check.HighProportionOfArtifactICs = 0;
% end
% issues_to_check.HighProportionOfArtifactICs = 0;

% % RS only
% if strcmpi(qa_data.subject.task,'RS') || strcmpi(qa_data.subject.task,'EO') || strcmpi(qa_data.subject.task,'EC')
%     % % EC eye blinks
%     % if qa_data.blockinfo.ec_blinks>0
%     %     issues_to_check.ECEyeBinksDetected = qa_data.blockinfo.ec_blinks;
%     % else
%     %     issues_to_check.ECEyeBinksDetected = 0;
%     % end
%
%     % Data lost due to epoching only
%     if qa_data.blockinfo.dataLostbyEpoching > 0.2
%         issues_to_check.RSdataLostbyEpoching = qa_data.blockinfo.dataLostbyEpoching;
%     else
%         issues_to_check.RSdataLostbyEpoching = 0;
%     end
% end

% % Leftovers
% if qa_data.leftovers.muscle2 > 0.25
%     issues_to_check.MuscleLeftovers = qa_data.leftovers.muscle2;
% else
%     issues_to_check.MuscleLeftovers = 0;
% end
% if abs(mean(qa_data.leftovers.blinks)-1)>0.05
%     issues_to_check.EyeLeftovers = mean(qa_data.leftovers.blinks)-1;
% else
%     issues_to_check.EyeLeftovers = 0;
% end

% % The participant will not be processed, likely cos data is missing
% issues_to_check.aFileName                              = qa_data.subject.id;
% issues_to_check.NumberTrials1                          = NaN;
% issues_to_check.NumberTrials2                          = NaN;
% issues_to_check.NumberIC1                              = NaN;
% issues_to_check.NumberIC2                              = NaN;
% issues_to_check.NumberIC3                              = NaN;
%
% % issues_to_check.badElectrodes                          = NaN;
% issues_to_check.FlatElectrodesDiscrepancy              = NaN;
% issues_to_check.RejectedTooManyElectrodes              = NaN;
% issues_to_check.HighProportionExcludedAsExtremeOutlier = NaN;
% issues_to_check.HighProportionOfEMG                    = NaN;
%
% % Note sure how to automate this part and not cheat
% MFWrounds = {'R1'}; % Max 4 MWF rounds of cleaning
% for j = 1:length(MFWrounds)
%     issues_to_check.(['MWF' MFWrounds{j} 'Status1']) =  NaN;
%     issues_to_check.(['MWF' MFWrounds{j} 'Status2']) =  NaN;
%     issues_to_check.(['MWF' MFWrounds{j} 'BadData']) =  NaN;
% end
%
% issues_to_check.DataTooShortForValidICA     = NaN;
% issues_to_check.HighProportionOfArtifactICs = NaN;
%
% if strcmpi(qa_data.subject.task,'RS') || strcmpi(qa_data.subject.task,'EO') || strcmpi(qa_data.subject.task,'EC')
%     % issues_to_check.ECEyeBinksDetected = NaN;
%     issues_to_check.RSdataLostbyEpoching = NaN;
% end
%
% issues_to_check.MuscleLeftovers1 = NaN;
% % issues_to_check.EyeLeftovers    = NaN;

% Log
qa_data.issues_to_check = issues_to_check;

end