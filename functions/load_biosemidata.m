function EEG = load_biosemidata(subject,myPaths)

fprintf('\n================================\n');
fprintf('Loading data\n');
fprintf('================================\n');

% Load data
if strcmpi(myPaths.task,'MT')
    % Motor tasks
    EEG = pop_biosig(subject.datablocks,'channels',1:168);
    assert(EEG(1).nbchan == 168, 'Unexpected number of channels');
    EEG = pop_select(EEG,'rmchannel',[131 132 143:160]);

else
    % % This is usually the case
    % EEG = pop_biosig(subject.datablocks,'channels',1:136);
    %
    % % Check if the number channels is as expected
    % if EEG(1).nbchan == 168
    %     % Rare cases of Utrect RS datasets recorded with Motor task settings
    %     EEG = pop_biosig(subject.datablocks,'channels',1:168);
    %     EEG = pop_select(EEG,'channel',[1:128 161:168]);
    %
    % elseif EEG(1).nbchan == 264
    %     % Some DUB RS datasets
    %     EEG = pop_biosig(subject.datablocks,'channels',1:264);
    %     EEG = pop_select(EEG,'channel',[1:128 257:264]);
    %
    % elseif EEG(1).nbchan == 271
    %     % Some DUB RS datasets
    %     EEG = pop_biosig(subject.datablocks,'channels',1:271);
    %     EEG = pop_select(EEG,'channel',[1:128 257:264]);
    %
    % elseif EEG(1).nbchan == 143
    %     % Some DUB RS datasets
    %     EEG = pop_biosig(subject.datablocks,'channels',1:271);
    %     EEG = pop_select(EEG,'channel',[1:128 257:264]);
    %
    % end

    % This is usually the case
    EEG = pop_biosig(subject.datablocks);
    NBLK = length(EEG);

    for i = 1:NBLK
        % Check if the number channels is as expected
        switch EEG(i).nbchan
            case 136
                % As expected

            case 168
                % Utrect RS datasets recorded with Motor task settings
                % EEG = pop_biosig(subject.datablocks,'channels',1:168);
                EEG(i) = pop_select(EEG(i),'channel',[1:128 161:168]);

            case 264
                % DUB RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels',1:264);
                EEG(i) = pop_select(EEG(i),'channel',[1:128 257:264]);

            case 271
                % DUB RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels',1:271);
                EEG(i) = pop_select(EEG(i),'channel',[1:128 257:264]);

            case 143
                % DUB RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels',1:143);
                EEG(i) = pop_select(EEG(i),'channel',1:136);

            otherwise
                error('Unexpected number of channels...');
        end
    end

    % Assure that data has fs = 512 Hz
    maskfs2k = [EEG(:).srate] > 512;
    if any(maskfs2k)
        pop_editoptions('option_parallel', 0);
        EEG(maskfs2k) = pop_resample(EEG(maskfs2k), 512);
        pop_editoptions('option_parallel', 1);
    end

    % Double-check
    assert(all([EEG(:).srate] == 512));
    assert(all([EEG(:).nbchan] == 136));
end

end