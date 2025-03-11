function EEG = do_baselinecorrection(EEG,typeMethod,myPaths)

fprintf('\n================================\n');
fprintf('Baseline correction (%s)\n',typeMethod);
fprintf('================================\n');

if strcmpi(typeMethod,'traditional')
    % Baseline correction: Traditional
    if strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO') || strcmpi(myPaths.task,'EC')
        EEG = pop_rmbase(EEG,[],[]);
    else
        EEG = pop_rmbase(EEG,[(EEG.xmin)*1000 0],[]);
    end

elseif strcmpi(typeMethod,'regression')
    % Baseline correction: Regression-based method
    % if strcmpi(myPaths.task,'MMN')
    %     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mmn{1},'Uniformoutput',0);
    %
    % elseif strcmpi(myPaths.task,'SART')
    %     % SART wrt visual stimuli
    %     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart1{1},'Uniformoutput',0);
    %     % SART wrt response times
    %     condLabel2 = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.sart2{1},'Uniformoutput',0);
    %
    % elseif strcmpi(myPaths.task,'MT')
    %     condLabel = arrayfun(@(x) ['condition ' num2str(x)],cfg.trg.mt{1},'Uniformoutput',0);
    %
    % end
    % if strcmpi(myPaths.task,'MMN') || strcmpi(myPaths.task,'SART')
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