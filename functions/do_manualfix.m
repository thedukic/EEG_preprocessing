function EEG = do_manualfix(EEG,subject,myPaths)
% C48 (SART) is very strage - low quality data?
% Note:
% This part can be automated by checking the electric potential duing 
% eyeblinks (that can be automatically detected). These potentials are very 
% strong, (almost) always present and have a very well defined potential
% distribution (topolot).

fprintf('\n================================\n');
fprintf('Manually fixing some datasets\n');
fprintf('================================\n');

% Define
NBLK = length(EEG);

% Swapped sets
if strcmp(subject.id,'C50') && strcmpi(myPaths.task,'SART')
    warning([subject.id 'has swapped C- and B- set. Fixing that now...']);
    for i = 1:NBLK
        EEG(i).data(33:96,:) = [EEG(i).data(65:96,:); EEG(i).data(33:64,:)];
    end

elseif strcmp(subject.id,'P49') && (strcmpi(myPaths.task,'RS') || strcmpi(myPaths.task,'EO'))
    warning([subject.id 'has swapped C- and D- set. Fixing that now...']);
    for i = 1:NBLK
        EEG(i).data(65:128,:) = [EEG(i).data(97:128,:); EEG(i).data(65:96,:)];
    end

else
    fprintf('Nice, not needed!\n');
end

% Check
EEG = eeg_checkset(EEG,'loaddata');

end