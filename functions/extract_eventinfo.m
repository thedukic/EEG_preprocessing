function [EEG, eventinfo] = extract_eventinfo(EEG, cfg)

% Pradigm
task = EEG(1).ALSUTRECHT.subject.task;

fprintf('\n================================\n');
fprintf('Extracting %s event information\n', task);
fprintf('================================\n');

% Keep event info
num_block = length(EEG);
eventinfo = cell(num_block, 4);
block_length = NaN(num_block, 1);

switch task
    case 'SART'
        for i_block = 1:num_block
            % Block length
            block_length(i_block) = round((EEG(i_block).pnts / EEG(i_block).srate) / 60, 1);

            % Fix response events
            for i_event = 1:numel(EEG(i_block).event)
                % First add 'condition' prefix if needed
                if ~contains(EEG(i_block).event(i_event).type, 'condition') && ~strcmpi(EEG(i_block).event(i_event).type, 'boundary')
                    EEG(i_block).event(i_event).edftype   = str2double(EEG(i_block).event(i_event).type);
                    EEG(i_block).urevent(i_event).edftype = EEG(i_block).event(i_event).edftype;
                    EEG(i_block).event(i_event).type      = ['condition ' EEG(i_block).event(i_event).type];
                    EEG(i_block).urevent(i_event).type    = EEG(i_block).event(i_event).type;
                end

                % Then convert all 11/12/2s to 1s
                if strcmp(EEG(i_block).event(i_event).type, 'condition 11') || strcmp(EEG(i_block).event(i_event).type, 'condition 12') || strcmp(EEG(i_block).event(i_event).type, 'condition 2')
                    EEG(i_block).event(i_event).edftype   = 1;
                    EEG(i_block).urevent(i_event).edftype = 1;
                    EEG(i_block).event(i_event).type      = 'condition 1';
                    EEG(i_block).urevent(i_event).type    = 'condition 1';
                end
            end

            % Double-checks
            EEG(i_block) = eeg_checkset(EEG(i_block));

            labels1 = unique({EEG(i_block).event(:).type});
            labels1 = labels1(~contains(labels1, '255'));
            labels1 = labels1(~contains(labels1, 'boundary'));
            assert(all(contains(labels1, 'condition')));

            % labels1 = strjoin(labels1, ', ');
            % fprintf('SART%d unique events:\n', i_block); disp(labels1);
            print_events(task, i_block, labels1);

            % Extract events and latencies
            for i_event = 1:numel(EEG(i_block).event)
                if strcmp(EEG(i_block).event(i_event).type, 'condition 1') || strcmp(EEG(i_block).event(i_event).type, 'condition 3') || strcmp(EEG(i_block).event(i_event).type, 'condition 6')
                    eventinfo{i_block, 1} = [eventinfo{i_block, 1}, EEG(i_block).event(i_event).edftype];
                    eventinfo{i_block, 2} = [eventinfo{i_block, 2}, EEG(i_block).event(i_event).latency];
                end
            end
            assert(length(eventinfo{i_block, 1}) == length(eventinfo{i_block, 2}));

            eventinfo{i_block, 3} = sum(eventinfo{i_block, 1}==3 | eventinfo{i_block, 1}==6);
            eventinfo{i_block, 4} = EEG(i_block).srate;
            fprintf('SART%d has %d trials (Go+NoGo) and length of %0.1f min.\n', i_block, eventinfo{i_block, 3}, block_length(i_block));
        end

    case 'MMN'
        for i_block = 1:num_block
            % Block length
            block_length(i_block) = round((EEG(i_block).pnts / EEG(i_block).srate) / 60, 1);

            % Fix events if needed
            labels1 = {EEG(i_block).event(:).type};
            if any(~contains(labels1, 'condition'))
                for i_event = 1:numel(EEG(i_block).event)
                    if ~contains(EEG(i_block).event(i_event).type, 'condition') && ~strcmpi(EEG(i_block).event(i_event).type, 'boundary')
                        EEG(i_block).event(i_event).edftype   = str2double(EEG(i_block).event(i_event).type);
                        EEG(i_block).urevent(i_event).edftype = EEG(i_block).event(i_event).edftype;
                        EEG(i_block).event(i_event).type      = ['condition ' EEG(i_block).event(i_event).type];
                        EEG(i_block).urevent(i_event).type    = EEG(i_block).event(i_event).type;
                    end
                end
            end

            % Double-checks
            EEG(i_block) = eeg_checkset(EEG(i_block));

            labels1 = unique({EEG(i_block).event(:).type});
            labels1 = labels1(~contains(labels1, '255'));
            labels1 = labels1(~contains(labels1, 'boundary'));
            assert(all(contains(labels1, 'condition')));

            % 255 - start of MMN
            % 4   - not sure
            % labels1 = labels1(~(contains(labels1,'255') | contains(labels1,'condition 4')));

            % labels1 = strjoin(labels1, ', ');
            % fprintf('MMN%d unique events:\n', i_block); disp(labels1);
            print_events(task, i_block, labels1);

            % Cannot be asserted. At least one dataset has two different tasks (MMN3 + SART1)
            % labels2 = cell2mat(cellfun(@(x) str2double(x(end-1:end)), labels1, 'UniformOutput', false));
            % assert(length(labels1) == 2); % 12 / 17

            % Extract events and latencies
            for i_event = 1:numel(EEG(i_block).event)
                if strcmp(EEG(i_block).event(i_event).type, 'condition 12') || strcmp(EEG(i_block).event(i_event).type, 'condition 17')
                    eventinfo{i_block, 1} = [eventinfo{i_block, 1}, EEG(i_block).event(i_event).edftype];
                    eventinfo{i_block, 2} = [eventinfo{i_block, 2}, EEG(i_block).event(i_event).latency];
                end
            end
            assert(length(eventinfo{i_block, 1})==length(eventinfo{i_block, 2}));

            eventinfo{i_block, 3} = sum(eventinfo{i_block, 1}==12 | eventinfo{i_block, 1}==17);
            eventinfo{i_block, 4} = EEG(i_block).srate;
            fprintf('MMN%d has %d trials (standard+deviant) and length of %0.1f min.\n', i_block, eventinfo{i_block, 3}, block_length(i_block));
        end

    case 'MT'
        for i_block = 1:num_block
            % Block length
            block_length(i_block) = round((EEG(i_block).pnts / EEG(i_block).srate) / 60, 1);

            % Fix events if needed
            labels1 = {EEG(i_block).event(:).type};
            % labels1(~contains(labels1,'condition')) = cellfun(@(x) ['condition ' num2str(x)],labels1(~contains(labels1,'condition')),'Uniformoutput',0);
            if any(~contains(labels1, 'condition'))
                for i_event = 1:numel(EEG(i_block).event)
                    if ~contains(EEG(i_block).event(i_event).type, 'condition') && ~strcmpi(EEG(i_block).event(i_event).type, 'boundary')
                        EEG(i_block).event(i_event).edftype   = str2double(EEG(i_block).event(i_event).type);
                        EEG(i_block).urevent(i_event).edftype = EEG(i_block).event(i_event).edftype;
                        EEG(i_block).event(i_event).type      = ['condition ' EEG(i_block).event(i_event).type];
                        EEG(i_block).urevent(i_event).type    = EEG(i_block).event(i_event).type;
                    end
                end
            end

            % Double-checks
            EEG(i_block) = eeg_checkset(EEG(i_block));

            labels1 = unique({EEG(i_block).event(:).type});
            labels1 = labels1(~contains(labels1, '255'));
            labels1 = labels1(~contains(labels1, 'boundary'));
            assert(all(contains(labels1, 'condition')));

            labels2 = cell2mat(cellfun(@(x) str2double(x(end-1:end)), labels1, 'UniformOutput', false));
            assert(length(labels1) == 3); % 20/30/50, 21/31/51, 22/32/52 == 3 per task/block
            print_events(task, i_block, labels1);

            % Extract events and latencies
            for i_event = 1:numel(EEG(i_block).event)
                if strcmp(EEG(i_block).event(i_event).type,labels1{1}) || strcmp(EEG(i_block).event(i_event).type,labels1{2}) || strcmp(EEG(i_block).event(i_event).type,labels1{3})
                    eventinfo{i_block, 1} = [eventinfo{i_block, 1}, EEG(i_block).event(i_event).edftype];
                    eventinfo{i_block, 2} = [eventinfo{i_block, 2}, EEG(i_block).event(i_event).latency];
                end
            end

            eventinfo{i_block, 3} = sum(eventinfo{i_block,1}==labels2(1) | eventinfo{i_block,1}==labels2(2) | eventinfo{i_block,1}==labels2(3)) / 3;
            eventinfo{i_block, 4} = EEG(i_block).srate;
            fprintf('MT%d has %d trials and length of %0.1f min.\n', i_block, floor(eventinfo{i_block, 3}), block_length(i_block));
        end

    case {'RS', 'EO', 'EC'}
        % Epoch length
        L_seconds = cfg.rs1{1};        % [seconds]
        L_samples = L_seconds * EEG(1).srate; % [s] -> [samples]

        % Overlap
        O_percent  = cfg.rs1{2};         % e.g. 0.5

        % Calculate step size
        step_size = L_samples * (1-O_percent);

        for i_block = 1:num_block
            % Block length
            N = EEG(i_block).pnts;
            block_length(i_block) = round((N / EEG(i_block).srate) / 60, 1);

            % N = length(EEG(i_block).times);
            % T = N ./ EEG(i_block).srate; fprintf('RS%d has %1.0fs of data.\n', i_block, T);

            % Calculate number of trials
            NTRL = floor((N - L_samples) / step_size) + 1;

            % eventinfo{i,3} = floor(N./(L.*EEG(i).srate));
            eventinfo{i_block, 3} = NTRL;
            eventinfo{i_block, 4} = EEG(i_block).srate;
            fprintf('RS%d has %d (%ds long, %1.2f overlap) trials and length of %0.1f min.\n', i_block, eventinfo{i_block, 3}, L_seconds, O_percent, block_length(i_block));
        end
end

% Report total
fprintf('\nTotal number of trials is %d with total length of %0.1f min.\n', floor(sum([eventinfo{:, 3}])), sum(block_length));

% Log
for i_block = 1:num_block
    EEG(i_block).ALSUTRECHT.eventinfo = eventinfo;
end

end


function print_events(task, i_block, labels1)
% Print directly using %s for the string array
fprintf('\n%s%d unique events: \n%s\n', task, i_block, strjoin(labels1, '\n'));
end