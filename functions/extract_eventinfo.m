function [EEG, eventinfo] = extract_eventinfo(EEG)

% Keep event info
NBLK = length(EEG);
eventinfo = cell(NBLK,4);

switch EEG(1).ALSUTRECHT.subject.task
    case 'SART'
        disp('Extracting SART event information...');

        for i = 1:NBLK
            % Fix response events
            for j = 1:numel(EEG(i).event)
                if strcmp(EEG(i).event(j).type,'condition 11') || strcmp(EEG(i).event(j).type,'condition 12') || strcmp(EEG(i).event(j).type,'condition 2')
                    EEG(i).event(j).edftype   = 1;
                    EEG(i).urevent(j).edftype = 1;
                    EEG(i).event(j).type      = 'condition 1';
                    EEG(i).urevent(j).type    = 'condition 1';
                end

                if ~contains(EEG(i).event(j).type,'condition')
                    EEG(i).event(j).edftype   = str2double(EEG(i).event(j).type);
                    EEG(i).urevent(j).edftype = EEG(i).event(j).edftype;
                    EEG(i).event(j).type      = ['condition ' EEG(i).event(j).type];
                    EEG(i).urevent(j).type    = EEG(i).event(j).type;
                end
            end

            % Double-checks
            EEG(i) = eeg_checkset(EEG(i));
            labels1 = unique({EEG(i).event(:).type});
            assert(all(contains(labels1,'condition')));

            labels1 = labels1(~contains(labels1,'255'));
            labels1 = strjoin(labels1,', ');
            fprintf('SART%d unique events:\n',i); disp(labels1);

            % Extract events and latencies
            for j = 1:numel(EEG(i).event)
                if strcmp(EEG(i).event(j).type,'condition 1') || strcmp(EEG(i).event(j).type,'condition 3') || strcmp(EEG(i).event(j).type,'condition 6')
                    eventinfo{i,1} = [eventinfo{i,1}, EEG(i).event(j).edftype];
                    eventinfo{i,2} = [eventinfo{i,2}, EEG(i).event(j).latency];
                end
            end
            assert(length(eventinfo{i,1})==length(eventinfo{i,2}));

            eventinfo{i,3} = sum(eventinfo{i,1}==3 | eventinfo{i,1}==6);
            eventinfo{i,4} = EEG(i).srate;
            fprintf('SART%d has %d trials (Go+NoGo).\n',i,eventinfo{i,3});
        end
    case 'MMN'
        disp('Extracting MMN event information...');

        for i = 1:NBLK
            % Fix events if needed
            labels1 = {EEG(i).event(:).type};
            if any(~contains(labels1,'condition'))
                for j = 1:numel(EEG(i).event)
                    if ~contains(EEG(i).event(j).type,'condition')
                        EEG(i).event(j).edftype   = str2double(EEG(i).event(j).type);
                        EEG(i).urevent(j).edftype = EEG(i).event(j).edftype;
                        EEG(i).event(j).type      = ['condition ' EEG(i).event(j).type];
                        EEG(i).urevent(j).type    = EEG(i).event(j).type;
                    end
                end
            end

            % Double-checks
            EEG(i) = eeg_checkset(EEG(i));
            labels1 = unique({EEG(i).event(:).type});
            assert(all(contains(labels1,'condition')));

            % 255 - start of MMN
            %   4 - not sure ????
            % labels1 = labels1(~(contains(labels1,'255') | contains(labels1,'condition 4')));
            labels1 = labels1(~contains(labels1,'255'));
            labels1 = strjoin(labels1,', ');
            fprintf('MMN%d unique events:\n',i); disp(labels1);
            % Has to be hard coded, at least one dataset has two tasks (MMN3+SART1)
            % labels2 = cell2mat(cellfun(@(x) str2double(x(end-1:end)),labels1,'UniformOutput',false));
            % assert(length(labels1)==2); % 12/17

            % Extract events and latencies
            for j = 1:numel(EEG(i).event)
                if strcmp(EEG(i).event(j).type,'condition 12') || strcmp(EEG(i).event(j).type,'condition 17')
                    eventinfo{i,1} = [eventinfo{i,1}, EEG(i).event(j).edftype];
                    eventinfo{i,2} = [eventinfo{i,2}, EEG(i).event(j).latency];
                end
            end
            assert(length(eventinfo{i,1})==length(eventinfo{i,2}));

            eventinfo{i,3} = sum(eventinfo{i,1}==12 | eventinfo{i,1}==17);
            eventinfo{i,4} = EEG(i).srate;
            fprintf('MMN%d has %d trials (standard+deviant).\n',i,eventinfo{i,3});
        end
    case 'MT'
        disp('Extracting MT event information...');

        for i = 1:NBLK
            % Fix events if needed
            labels1 = {EEG(i).event(:).type};
            % labels1(~contains(labels1,'condition')) = cellfun(@(x) ['condition ' num2str(x)],labels1(~contains(labels1,'condition')),'Uniformoutput',0);
            if any(~contains(labels1,'condition'))
                for j = 1:numel(EEG(i).event)
                    if ~contains(EEG(i).event(j).type,'condition')
                        EEG(i).event(j).edftype   = str2double(EEG(i).event(j).type);
                        EEG(i).urevent(j).edftype = EEG(i).event(j).edftype;
                        EEG(i).event(j).type      = ['condition ' EEG(i).event(j).type];
                        EEG(i).urevent(j).type    = EEG(i).event(j).type;
                    end
                end
            end

            % Double-checks
            EEG = eeg_checkset(EEG);
            labels1 = unique({EEG(i).event(:).type});
            assert(all(contains(labels1,'condition')));
            labels2 = cell2mat(cellfun(@(x) str2double(x(end-1:end)),labels1,'UniformOutput',false));
            assert(length(labels1)==3); % 20/30/50, 21/31/51, 22/32/52 == 3 per task/block

            % Extract events and latencies
            for j = 1:numel(EEG(i).event)
                if strcmp(EEG(i).event(j).type,labels1{1}) || strcmp(EEG(i).event(j).type,labels1{2}) || strcmp(EEG(i).event(j).type,labels1{3})
                    eventinfo{i,1} = [eventinfo{i,1}, EEG(i).event(j).edftype];
                    eventinfo{i,2} = [eventinfo{i,2}, EEG(i).event(j).latency];
                end
            end

            eventinfo{i,3} = sum(eventinfo{i,1}==labels2(1) | eventinfo{i,1}==labels2(2) | eventinfo{i,1}==labels2(3))/3;
            eventinfo{i,4} = EEG(i).srate;
            fprintf('MT%d has %d trials.\n',i,floor(eventinfo{i,3}));
        end
    case {'RS','EO','EC'}
        fprintf('Resting-state data does not have events by default. Returning only the number of possible 1s trials in each block.\n');
        L = 1; % length in [s]
        for i = 1:NBLK
            eventinfo{i,3} = floor(length(EEG(i).times)./(L.*EEG(i).srate));
            eventinfo{i,4} = EEG(i).srate;
            fprintf('RS%d has %d (%ds long) trials.\n',i,eventinfo{i,3},L);
        end
end

% Log
for i = 1:NBLK
    EEG(i).ALSUTRECHT.eventinfo = eventinfo;
end

end