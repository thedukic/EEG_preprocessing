function EEG = do_manualfix(EEG,subject,myPaths)
% C48 is very strage - low quality data?

fprintf('\n================================\n');
fprintf('Manually fixing some datasets\n');
fprintf('================================\n');

% 1. Swapped C- and B-set
if strcmp(subject.id,'C50') && strcmpi(myPaths.task,'SART')
    warning([subject.id 'has swapped C- and B- set. Fixing that now...']);
    
    NBLK = length(EEG);
    for i = 1:NBLK
        EEG(i).data(33:96,:) = [EEG(i).data(65:96,:); EEG(i).data(33:64,:)];
    end
else
    fprintf('Nice, not needed!\n');
end

% Check 
EEG = eeg_checkset(EEG,'loaddata');

end