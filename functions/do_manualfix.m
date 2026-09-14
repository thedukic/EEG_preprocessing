function EEG = do_manualfix(EEG)
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
num_blocks = length(EEG);
id   = EEG(1).ALSUTRECHT.subject.id;
task = EEG(1).ALSUTRECHT.subject.task;

% Swapped electrode sets
if strcmp(id, 'C50') && strcmpi(task, 'SART')
    warning([id 'has swapped C- and B- set. Fixing that now...']);
    for i_block = 1:num_blocks
        EEG(i_block).data(33:96, :) = [EEG(i_block).data(65:96, :); EEG(i_block).data(33:64, :)];
    end

elseif strcmp(id, 'P49') && (strcmpi(task, 'RS') || strcmpi(task, 'EO'))
    warning([id 'has swapped C- and D- set. Fixing that now...']);
    for i_block = 1:num_blocks
        EEG(i_block).data(65:128, :) = [EEG(i_block).data(97:128, :); EEG(i_block).data(65:96, :)];
    end

else
    fprintf('Nice, not needed!\n');
end

% Check
EEG = eeg_checkset(EEG, 'loaddata');

end