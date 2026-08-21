function EEG = remove_datasetends(EEG)
%
% % EEGLAB data struct input
% SART events (1/11, 3/6, 40/30, 10/20)
% MMN  events (17 & 12)
%
fprintf('\n================================\n');
fprintf('Cutting the ends of each block\n');
fprintf('================================\n');
fprintf('Making sure that the data are divisible into 1s-long epochs.\n');

NBLK = length(EEG);
taskType = EEG(1).ALSUTRECHT.subject.task;

switch taskType
    case 'SART'
        for i_blk = 1:NBLK
            % lat = NaN(1,2);
            % event   = {EEG(i_blk).event.type};
            % latency = {EEG(i_blk).event.latency};
            %
            % % First event
            % indx = find(ismember(event,'condition 6'),1,'first');
            % lat(1) = latency{indx}./EEG(i_blk).srate - 1; % seconds
            %
            % % Last event
            % % if strcmpi(event{end},'40') || strcmpi(event{end},'condition 30') || strcmpi(event{end},'condition 10') || strcmpi(event{end},'condition 20')
            % %     lat(2) = latency{end}./EEG(i).srate + 1; % seconds
            % % else
            % %     error('Unexpeced!');
            % % end
            % indx1 = find(ismember(event,'condition 3'),1,'last');
            % indx2 = find(ismember(event,'condition 6'),1,'last');
            % indx = max(indx1,indx2);
            % lat(2) = latency{indx}./EEG(i_blk).srate + 2; % seconds
            %
            % % Just in case there is not enough buffer after the last trial
            % if lat(2) > EEG(i_blk).times(end)/1000-1
            %     lat(2) = EEG(i_blk).times(end)/1000-1;
            % end
            %
            % % Remove the ends
            % % lat = round(lat,1);
            % timedif = mod(lat,1);
            % if timedif(1) > timedif(2)
            %     lat(2) = lat(2) + abs(diff(timedif));
            % else
            %     lat(2) = floor(lat(2)) + 1 + timedif(1);
            % end
            % assert(mod(diff(lat),1) == 0);
            % lat(2) = lat(2) - 1/EEG(i_blk).srate;

            lat = extract_lats(EEG(i_blk), 'SART');
            EEG(i_blk) = pop_select(EEG(i_blk),'time',lat);
            check_datasize(EEG(i_blk));
        end

    case 'MMN'
        for i_blk = 1:NBLK

            % lat = NaN(1,2);
            % event   = {EEG(i_blk).event.type};
            % latency = {EEG(i_blk).event.latency};
            %
            % % First event
            % indx = find(ismember(event,'condition 12'),1,'first');
            % lat(1) = latency{indx}./EEG(i_blk).srate - 1; % seconds
            %
            % % Last event
            % indx1 = find(ismember(event,'condition 12'),1,'last');
            % indx2 = find(ismember(event,'condition 17'),1,'last');
            % indx = max(indx1,indx2);
            % lat(2) = latency{indx}./EEG(i_blk).srate + 2; % seconds
            %
            % % Just in case there is not enough buffer after the last trial
            % if lat(2) > EEG(i_blk).times(end)/1000-1
            %     lat(2) = EEG(i_blk).times(end)/1000-1;
            % end
            %
            % % Remove the ends
            % % lat = round(lat,1);
            % timedif = mod(lat,1);
            % if timedif(1) > timedif(2)
            %     lat(2) = lat(2) + abs(diff(timedif));
            % else
            %     lat(2) = floor(lat(2)) + 1 + timedif(1);
            % end
            % assert(mod(diff(lat),1) == 0);
            % lat(2) = lat(2)-1/EEG(i_blk).srate;

            lat = extract_lats(EEG(i_blk), 'MMN');
            EEG(i_blk) = pop_select(EEG(i_blk),'time',lat);
            check_datasize(EEG(i_blk));
        end

    case {'RS','EO','EC'}
        T = 2 * EEG(1).srate;
        % for i = 1:NTRL
        %     EEG.trial{i} = EEG.trial{i}(:,T+1:end-T);
        %     EEG.time{i}  = EEG.time{i}(1:end-2*T);
        % end

        % Make sure data is rounded to N*1s
        % Cuts T = 2s of data from both ends
        for i_blk = 1:NBLK
            N = floor(EEG(i_blk).pnts/EEG(i_blk).srate) * EEG(i_blk).srate - T;
            EEG(i_blk) = pop_select(EEG(i_blk),'point',[T+1 N]);
        end

    case 'MT'
        for i_blk = 1:NBLK

            % lat = NaN(1,2);
            % event   = {EEG(i_blk).event.type};
            % latency = {EEG(i_blk).event.latency};
            %
            % mask = ~ismember(event, 'boundary');
            % event = event(mask);
            % latency = event(latency);
            %
            % % First event
            % lat(1) = latency{1}./EEG(i_blk).srate - 1;   % seconds
            %
            % % Last event
            % lat(2) = latency{end}./EEG(i_blk).srate + 2; % seconds
            %
            % % MT2 often does not have buffer after the last trial.
            % % Enforce that lat(2) is never greater than the time limit (EEG end - 1s).
            % lat(2) = min(lat(2), EEG(i_blk).times(end)/1000 - 1);
            %
            %
            % % % Enforce that lat(1) is never greater than the time limit (EEG end - 1s).
            % % lat(1) = min(lat(1), 0);

            % % Remove the ends
            % % lat = round(lat,1);
            % timedif = mod(lat,1);
            % if timedif(1) > timedif(2)
            %     lat(2) = lat(2) + abs(diff(timedif));
            % else
            %     lat(2) = floor(lat(2)) + 1 + timedif(1);
            % end
            % assert(mod(diff(lat),1) == 0);
            % lat(2) = lat(2) - 1./EEG(i_blk).srate;
            %
            % % [a,b] = min(abs(EEG(i).times/1000-lat(1)));
            % % lat = lat + (EEG(i).times(b)/1000-lat(1));
            % % lat(1) = EEG(i).times(b)/1000;

            lat = extract_lats(EEG(i_blk), 'MT');
            EEG(i_blk) = pop_select(EEG(i_blk),'time',lat);
            check_datasize(EEG(i_blk));
        end
end

% Check
EEG = eeg_checkset(EEG);

fprintf('Done!\n');

end

% =========================================================================
% Helper functions
% =========================================================================
function lat = extract_lats(EEG, taskType)

% Extract
event   = {EEG.event.type};
latency = {EEG.event.latency};

% Remove boundary events
mask = ~ismember(event, 'boundary');
event = event(mask);
latency = latency(mask);

% Each task has their own event numbers
switch taskType
    case 'SART'
        indx1 = find(ismember(event,'condition 6'),1,'first');

        indx2a = find(ismember(event,'condition 3'),1,'last');
        indx2b = find(ismember(event,'condition 6'),1,'last');
        indx2 = max(indx2a,indx2b);

    case 'MMN'
        indx1 = find(ismember(event,'condition 12'),1,'first');

        indx2a = find(ismember(event,'condition 12'),1,'last');
        indx2b = find(ismember(event,'condition 17'),1,'last');
        indx2 = max(indx2a,indx2b);

    case 'MT'
        % Convert cell array of strings to a standard string array
        stringArray = string(event);

        % Extract only the number part (e.g., '50', '51', '52', '30')
        numberStrings = extractAfter(stringArray, "condition ");

        % Convert the string numbers to a numeric array
        conditionNumbers = str2double(numberStrings);

        % Calculate the condition type (the base 10 value)
        conditionType = floor(conditionNumbers / 10) * 10;

        % Take the most common one
        conditionType = mode(conditionType);

        % Double-check
        assert(ismember(conditionType,  [20, 30, 50]));

        indx1 = find(ismember(event,['condition ' num2str(conditionType)]),1,'first');
        indx2 = find(ismember(event,['condition ' num2str(conditionType + 2)]),1,'last');

end

% Allocate
lat = NaN(1,2);

% First event
lat(1) = latency{indx1}./EEG.srate - 1; % seconds

% Last event
lat(2) = latency{indx2}./EEG.srate + 2; % seconds

% MT2 often does not have buffer after the last trial.
% Enforce that lat(2) is never greater than the time limit (EEG end - 1s).
lat(2) = min(lat(2), EEG.times(end)/1000 - 1);

% Enforce that lat(1) is never greater than the time limit (EEG end - 1s).
lat(1) = max(lat(1), 0);

% Remove the ends
% lat = round(lat,1);
timedif = mod(lat,1);
if timedif(1) > timedif(2)
    lat(2) = lat(2) + abs(diff(timedif));
else
    lat(2) = floor(lat(2)) + 1 + timedif(1);
end
assert(mod(diff(lat),1) == 0);
lat(2) = lat(2) - 1 ./ EEG.srate;

end

function check_datasize(EEG)
L = EEG.srate;
N = floor(size(EEG.data,2) ./ L);
assert(size(EEG.data,2) == N*L);
end