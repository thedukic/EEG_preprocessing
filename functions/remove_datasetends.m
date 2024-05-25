function EEG = remove_datasetends(EEG)
%
% % EEGLAB data struct input
% SART events (1/11, 3/6, 40/30, 10/20)
% MMN  events (17 & 12)
%

fprintf('Cutting the ends of each block...\n');
fprintf('Making sure that the data is divisible into 1s epochs.\n');
NTRL = length(EEG);

switch EEG(1).ALSUTRECHT.subject.task
    case 'SART'
        for i = 1:NTRL
            lat = NaN(1,2);
            event   = {EEG(i).event.type};
            latency = {EEG(i).event.latency};

            % First event
            indx = find(ismember(event,'condition 6'),1,'first');
            lat(1) = latency{indx}./EEG(i).srate - 1; % seconds

            % Last event
            % if strcmpi(event{end},'40') || strcmpi(event{end},'condition 30') || strcmpi(event{end},'condition 10') || strcmpi(event{end},'condition 20')
            %     lat(2) = latency{end}./EEG(i).srate + 1; % seconds
            % else
            %     error('Unexpeced!');
            % end
            indx1 = find(ismember(event,'condition 3'),1,'last');
            indx2 = find(ismember(event,'condition 6'),1,'last');
            indx = max(indx1,indx2);
            lat(2) = latency{indx}./EEG(i).srate + 2; % seconds

            % Just in case there is not enough buffer after the last trial
            if lat(2)>EEG(i).times(end)/1000-1
                lat(2) = EEG(i).times(end)/1000-1;
            end

            % Remove the ends
            % lat = round(lat,1);
            timedif = mod(lat,1);
            if timedif(1)>timedif(2)
                lat(2) = lat(2)+abs(diff(timedif));
            else
                lat(2) = floor(lat(2))+1+(timedif(1));
            end
            assert(mod(diff(lat),1)==0);
            lat(2) = lat(2)-1/EEG(1).srate;
            EEG(i) = pop_select(EEG(i),'time',lat);

            L = EEG(i).srate;
            N = floor(size(EEG(i).data,2)/L);
            assert(N*L==size(EEG(i).data,2));
        end
    case 'MMN'
        for i = 1:NTRL
            lat = NaN(1,2);
            event   = {EEG(i).event.type};
            latency = {EEG(i).event.latency};

            % First event
            indx = find(ismember(event,'condition 12'),1,'first');
            lat(1) = latency{indx}./EEG(i).srate - 1; % seconds

            % Last event
            indx1 = find(ismember(event,'condition 12'),1,'last');
            indx2 = find(ismember(event,'condition 17'),1,'last');
            indx = max(indx1,indx2);
            lat(2) = latency{indx}./EEG(i).srate + 2; % seconds

            % Just in case there is not enough buffer after the last trial
            if lat(2)>EEG(i).times(end)/1000-1
                lat(2) = EEG(i).times(end)/1000-1;
            end

            % Remove the ends
            % lat = round(lat,1);
            timedif = mod(lat,1);
            if timedif(1)>timedif(2)
                lat(2) = lat(2)+abs(diff(timedif));
            else
                lat(2) = floor(lat(2))+1+(timedif(1));
            end
            assert(mod(diff(lat),1)==0);
            lat(2) = lat(2)-1/EEG(1).srate;
            EEG(i) = pop_select(EEG(i),'time',lat);

            L = EEG(i).srate;
            N = floor(size(EEG(i).data,2)/L);
            assert(N*L==size(EEG(i).data,2));
        end
    case {'RS','EO','EC'}
        T = 2*EEG(1).srate;
        % for i = 1:NTRL
        %     EEG.trial{i} = EEG.trial{i}(:,T+1:end-T);
        %     EEG.time{i}  = EEG.time{i}(1:end-2*T);
        % end

        % Make sure data is rounded to N*1s
        % Cuts T = 2s of data from both ends
        for i = 1:NTRL
            N = floor(EEG(i).pnts/EEG(i).srate)*EEG(i).srate-T;
            EEG(i) = pop_select(EEG(i),'point',[T+1 N]);
        end
    case 'MT'
        for i = 1:NTRL
            lat = NaN(1,2);
            % event   = {EEG(i).event.type};
            latency = {EEG(i).event.latency};

            % First event
            lat(1) = latency{1}./EEG(i).srate - 1;   % seconds

            % Last event
            lat(2) = latency{end}./EEG(i).srate + 2; % seconds

            % MT2 often does not have buffer after the last trial
            if lat(2)>EEG(i).times(end)/1000-1
                lat(2) = EEG(i).times(end)/1000-1;
            end

            % Remove the ends
            % lat = round(lat,1);
            timedif = mod(lat,1);
            if timedif(1)>timedif(2)
                lat(2) = lat(2)+abs(diff(timedif));
            else
                lat(2) = floor(lat(2))+1+(timedif(1));
            end
            assert(mod(diff(lat),1)==0);
            lat(2) = lat(2)-1/EEG(i).srate;

            % [a,b] = min(abs(EEG(i).times/1000-lat(1)));
            % lat = lat + (EEG(i).times(b)/1000-lat(1));
            % lat(1) = EEG(i).times(b)/1000;

            EEG(i) = pop_select(EEG(i),'time',lat);

            L = EEG(i).srate;
            N = floor(size(EEG(i).data,2)/L);
            assert(N*L==size(EEG(i).data,2));
        end
end

% Check
EEG = eeg_checkset(EEG);

end