function [EEG, eventinfo] = extract_eventinfo(EEG,cfg)

fprintf('\n================================\n');
fprintf('Extracting %s event information\n',EEG(1).ALSUTRECHT.subject.task);
fprintf('================================\n');

% Keep event info
NBLK = length(EEG);
eventinfo = cell(NBLK,4);

switch EEG(1).ALSUTRECHT.subject.task
    case 'SART'
        for i_blk = 1:NBLK
            % Fix response events
            for i_evt = 1:numel(EEG(i_blk).event)
                if strcmp(EEG(i_blk).event(i_evt).type,'condition 11') || strcmp(EEG(i_blk).event(i_evt).type,'condition 12') || strcmp(EEG(i_blk).event(i_evt).type,'condition 2')
                    EEG(i_blk).event(i_evt).edftype   = 1;
                    EEG(i_blk).urevent(i_evt).edftype = 1;
                    EEG(i_blk).event(i_evt).type      = 'condition 1';
                    EEG(i_blk).urevent(i_evt).type    = 'condition 1';
                end

                if ~contains(EEG(i_blk).event(i_evt).type,'condition') && ~strcmpi(EEG(i_blk).event(i_evt).type,'boundary')
                    EEG(i_blk).event(i_evt).edftype   = str2double(EEG(i_blk).event(i_evt).type);
                    EEG(i_blk).urevent(i_evt).edftype = EEG(i_blk).event(i_evt).edftype;
                    EEG(i_blk).event(i_evt).type      = ['condition ' EEG(i_blk).event(i_evt).type];
                    EEG(i_blk).urevent(i_evt).type    = EEG(i_blk).event(i_evt).type;
                end
            end

            % Double-checks
            EEG(i_blk) = eeg_checkset(EEG(i_blk));

            labels1 = unique({EEG(i_blk).event(:).type});
            labels1 = labels1(~contains(labels1,'255'));
            labels1 = labels1(~contains(labels1,'boundary'));
            assert(all(contains(labels1,'condition')));

            labels1 = strjoin(labels1,', ');
            fprintf('SART%d unique events:\n',i_blk); disp(labels1);

            % Extract events and latencies
            for i_evt = 1:numel(EEG(i_blk).event)
                if strcmp(EEG(i_blk).event(i_evt).type,'condition 1') || strcmp(EEG(i_blk).event(i_evt).type,'condition 3') || strcmp(EEG(i_blk).event(i_evt).type,'condition 6')
                    eventinfo{i_blk,1} = [eventinfo{i_blk,1}, EEG(i_blk).event(i_evt).edftype];
                    eventinfo{i_blk,2} = [eventinfo{i_blk,2}, EEG(i_blk).event(i_evt).latency];
                end
            end
            assert(length(eventinfo{i_blk,1}) == length(eventinfo{i_blk,2}));

            eventinfo{i_blk,3} = sum(eventinfo{i_blk,1}==3 | eventinfo{i_blk,1}==6);
            eventinfo{i_blk,4} = EEG(i_blk).srate;
            fprintf('SART%d has %d trials (Go+NoGo).\n',i_blk,eventinfo{i_blk,3});
        end

    case 'MMN'
        for i_blk = 1:NBLK
            % Fix events if needed
            labels1 = {EEG(i_blk).event(:).type};
            if any(~contains(labels1,'condition'))
                for i_evt = 1:numel(EEG(i_blk).event)
                    if ~contains(EEG(i_blk).event(i_evt).type,'condition') && ~strcmpi(EEG(i_blk).event(i_evt).type,'boundary')
                        EEG(i_blk).event(i_evt).edftype   = str2double(EEG(i_blk).event(i_evt).type);
                        EEG(i_blk).urevent(i_evt).edftype = EEG(i_blk).event(i_evt).edftype;
                        EEG(i_blk).event(i_evt).type      = ['condition ' EEG(i_blk).event(i_evt).type];
                        EEG(i_blk).urevent(i_evt).type    = EEG(i_blk).event(i_evt).type;
                    end
                end
            end

            % Double-checks
            EEG(i_blk) = eeg_checkset(EEG(i_blk));

            labels1 = unique({EEG(i_blk).event(:).type});
            labels1 = labels1(~contains(labels1,'255'));
            labels1 = labels1(~contains(labels1,'boundary'));
            assert(all(contains(labels1,'condition')));
            % 255 - start of MMN
            %   4 - not sure ????
            % labels1 = labels1(~(contains(labels1,'255') | contains(labels1,'condition 4')));

            labels1 = strjoin(labels1,', ');
            fprintf('MMN%d unique events:\n',i_blk); disp(labels1);
            % Has to be hard coded, at least one dataset has two tasks (MMN3+SART1)
            % labels2 = cell2mat(cellfun(@(x) str2double(x(end-1:end)),labels1,'UniformOutput',false));
            % assert(length(labels1)==2); % 12/17

            % Extract events and latencies
            for i_evt = 1:numel(EEG(i_blk).event)
                if strcmp(EEG(i_blk).event(i_evt).type,'condition 12') || strcmp(EEG(i_blk).event(i_evt).type,'condition 17')
                    eventinfo{i_blk,1} = [eventinfo{i_blk,1}, EEG(i_blk).event(i_evt).edftype];
                    eventinfo{i_blk,2} = [eventinfo{i_blk,2}, EEG(i_blk).event(i_evt).latency];
                end
            end
            assert(length(eventinfo{i_blk,1})==length(eventinfo{i_blk,2}));

            eventinfo{i_blk,3} = sum(eventinfo{i_blk,1}==12 | eventinfo{i_blk,1}==17);
            eventinfo{i_blk,4} = EEG(i_blk).srate;
            fprintf('MMN%d has %d trials (standard+deviant).\n',i_blk,eventinfo{i_blk,3});
        end

    case 'MT'
        for i_blk = 1:NBLK
            % Fix events if needed
            labels1 = {EEG(i_blk).event(:).type};
            % labels1(~contains(labels1,'condition')) = cellfun(@(x) ['condition ' num2str(x)],labels1(~contains(labels1,'condition')),'Uniformoutput',0);
            if any(~contains(labels1,'condition'))
                for i_evt = 1:numel(EEG(i_blk).event)
                    if ~contains(EEG(i_blk).event(i_evt).type,'condition') && ~strcmpi(EEG(i_blk).event(i_evt).type,'boundary')
                        EEG(i_blk).event(i_evt).edftype   = str2double(EEG(i_blk).event(i_evt).type);
                        EEG(i_blk).urevent(i_evt).edftype = EEG(i_blk).event(i_evt).edftype;
                        EEG(i_blk).event(i_evt).type      = ['condition ' EEG(i_blk).event(i_evt).type];
                        EEG(i_blk).urevent(i_evt).type    = EEG(i_blk).event(i_evt).type;
                    end
                end
            end

            % Double-checks
            EEG(i_blk) = eeg_checkset(EEG(i_blk));

            labels1 = unique({EEG(i_blk).event(:).type});
            labels1 = labels1(~contains(labels1,'255'));
            labels1 = labels1(~contains(labels1,'boundary'));
            assert(all(contains(labels1,'condition')));

            labels2 = cell2mat(cellfun(@(x) str2double(x(end-1:end)),labels1,'UniformOutput',false));
            assert(length(labels1) == 3); % 20/30/50, 21/31/51, 22/32/52 == 3 per task/block

            % Extract events and latencies
            for i_evt = 1:numel(EEG(i_blk).event)
                if strcmp(EEG(i_blk).event(i_evt).type,labels1{1}) || strcmp(EEG(i_blk).event(i_evt).type,labels1{2}) || strcmp(EEG(i_blk).event(i_evt).type,labels1{3})
                    eventinfo{i_blk,1} = [eventinfo{i_blk,1}, EEG(i_blk).event(i_evt).edftype];
                    eventinfo{i_blk,2} = [eventinfo{i_blk,2}, EEG(i_blk).event(i_evt).latency];
                end
            end

            eventinfo{i_blk,3} = sum(eventinfo{i_blk,1}==labels2(1) | eventinfo{i_blk,1}==labels2(2) | eventinfo{i_blk,1}==labels2(3))/3;
            eventinfo{i_blk,4} = EEG(i_blk).srate;
            fprintf('MT%d has %d trials.\n',i_blk,floor(eventinfo{i_blk,3}));
        end

    case {'RS','EO','EC'}
        % Epoch length
        L0 = cfg.rs{1};         % [seconds]
        L  = L0 * EEG(1).srate; % [s] -> [samples]

        % Overlap
        O  = cfg.rs{2};         % e.g. 0.5

        % Calculate step size
        step_size = L * (1-O);

        for i_blk = 1:NBLK
            N = length(EEG(i_blk).times);
            T = N ./ EEG(i_blk).srate;
            fprintf('RS%d has %1.0fs of data.\n',i_blk,T);

            % Calculate number of trials
            NTRL = floor((N - L) / step_size) + 1;

            % eventinfo{i,3} = floor(N./(L.*EEG(i).srate));
            eventinfo{i_blk,3} = NTRL;
            eventinfo{i_blk,4} = EEG(i_blk).srate;
            fprintf('RS%d has %d (%ds long, %1.2f overlap) trials.\n',i_blk,eventinfo{i_blk,3},L0,O);
        end
end

% Log
for i_blk = 1:NBLK
    EEG(i_blk).ALSUTRECHT.eventinfo = eventinfo;
end

end