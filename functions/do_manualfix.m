function EEG = do_manualfix(EEG, subject)
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

% Swapped electrode sets
if strcmp(subject.id, 'C50') && strcmpi(subject.task, 'SART')
    warning([subject.id 'has swapped C- and B- set. Fixing that now...']);
    for i_block = 1:NBLK
        EEG(i_block).data(33:96, :) = [EEG(i_block).data(65:96, :); EEG(i_block).data(33:64, :)];
    end

elseif strcmp(subject.id, 'P49') && (strcmpi(subject.task, 'RS') || strcmpi(subject.task, 'EO'))
    warning([subject.id 'has swapped C- and D- set. Fixing that now...']);
    for i_block = 1:NBLK
        EEG(i_block).data(65:128, :) = [EEG(i_block).data(97:128, :); EEG(i_block).data(65:96, :)];
    end

else
    fprintf('Nice, not needed!\n');
end

% Check
EEG = eeg_checkset(EEG, 'loaddata');

end