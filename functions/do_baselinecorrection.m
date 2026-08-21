function EEG = do_baselinecorrection(EEG, typeMethod)

fprintf('\n================================\n');
fprintf('Baseline-correcting data (%s)\n', typeMethod);
fprintf('================================\n');

if strcmpi(typeMethod, 'traditional')
    % Traditional
    if matches(EEG.ALSUTRECHT.subject.task, {'RS', 'EO', 'EC'}, 'IgnoreCase', true)
        EEG = pop_rmbase(EEG, [], []);
    else
        EEG = pop_rmbase(EEG, [(EEG.xmin)*1000 0], []);
    end

elseif strcmpi(typeMethod, 'regression')
    % % Regression-based method
    % if strcmpi(EEG.ALSUTRECHT.subject.task,'MMN')
    %     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0);
    %
    % elseif strcmpi(EEG.ALSUTRECHT.subject.task,'SART')
    %     % SART wrt visual stimuli
    %     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0);
    %     % SART wrt response times
    %     condLabel2 = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0);
    %
    % elseif strcmpi(EEG.ALSUTRECHT.subject.task,'MT')
    %     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0);
    %
    % end
    % if strcmpi(EEG.ALSUTRECHT.subject.task,'MMN') || strcmpi(EEG.ALSUTRECHT.subject.task,'SART')
    %     EEGcell = correct_baseline(EEGcell,[(EEGcell.xmin)*1000 0],'Factor_1_Level_1',condLabel);
    % else
    %     error('Check this step for MT data.');
    %     % How to consider MT2 and MT3 as one condition and MT5 as the other
    %     % Make a copy and rename MT3 to MT2 in the .event struct?
    %     % EEGcell = correct_baseline(EEGcell,[(EEGcell.xmin)*1000 0]); % if only 1 stimulus condition present
    % end
else
    error('Method not defined');
end

end