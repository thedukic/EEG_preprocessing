function DATA = load_biosemidata_old(subject)

fprintf('\n================================\n');
fprintf('Loading data\n');
fprintf('================================\n');

% Load data
if strcmpi(subject.task, 'MT')
    % Motor tasks
    DATA = pop_biosig(subject.datablocks, 'channels', 1:168);
    assert(DATA(1).nbchan == 168, 'Unexpected number of channels');

    % fprintf('\nRemoving empty channels:\n');
    % Remove the unused EMG channels
    % 1. The two unused EMG channel for ADM (131-132)
    % 2. The rest of the unused EMG channels (143-160)
    % Loop through the datasets to avoid memory thrashing/throttling
    % for i = 1:length(DATA)
    %     DATA(i) = pop_select(DATA(i), 'rmchannel', [131:132, 143:160]);
    % end

    % % CMC tests
    % [cmc_topo, f_axis, cmc_all, cl95, fh] = compute_cmc_topography(DATA(1));
    % res = plot_cmc_rect_comparison(DATA(1));

    % Remove all EMG data as they will be processed separtely
    for i = 1:length(DATA)
        DATA(i) = pop_select(DATA(i), 'rmchannel', 129:160);
    end

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
    %     % Dublin RS datasets
    %     EEG = pop_biosig(subject.datablocks,'channels',1:264);
    %     EEG = pop_select(EEG,'channel',[1:128 257:264]);
    %
    % elseif EEG(1).nbchan == 271
    %     % Dublin RS datasets
    %     EEG = pop_biosig(subject.datablocks,'channels',1:271);
    %     EEG = pop_select(EEG,'channel',[1:128 257:264]);
    %
    % elseif EEG(1).nbchan == 143
    %     % Dublin RS datasets
    %     EEG = pop_biosig(subject.datablocks,'channels',1:271);
    %     EEG = pop_select(EEG,'channel',[1:128 257:264]);
    %
    % end

    % Load EEG data
    DATA = pop_biosig(subject.datablocks);
    num_blocks = length(DATA);

    % Check if the number channels is as expected
    fprintf('\nChecking for empty electrodes in each block...\n');
    for i_block = 1:num_blocks
        switch DATA(i_block).nbchan
            case 136
                % As expected

            case 168
                % Utrect RS datasets recorded with Motor task settings
                % EEG = pop_biosig(subject.datablocks,'channels', 1:168);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', [1:128 161:168]);

            case 264
                % Dublin RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels', 1:264);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', [1:128 257:264]);

            case 271
                % Dublin RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels', 1:271);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', [1:128 257:264]);

            case 143
                % Dublin RS datasets
                % EEG = pop_biosig(subject.datablocks,'channels', 1:143);
                DATA(i_block) = pop_select(DATA(i_block), 'channel', 1:136);

            otherwise
                error('Unexpected number of channels!');
        end
    end

    assert(all([DATA(:).nbchan] == 136));
    fprintf('Done!\n');

end

% Ensure that data is with fs = 512 Hz
fprintf('Checking the sampling rate in each block...\n');
srate_all = [DATA(:).srate] > 512;
idx_resample = find(srate_all);

if ~isempty(idx_resample)
    % % Temporarily disable EEGLAB parallel processing
    % pop_editoptions('option_parallel', 0);

    % Loop through only the datasets that require downsampling
    for i = 1:length(idx_resample)
        idx = idx_resample(i);
        fprintf('Downsampling dataset %d of %d (from %d Hz to 512 Hz)...\n', i, length(idx_resample), DATA(idx).srate);

        % Process sequentially to keep RAM usage stable
        DATA(idx) = pop_resample(DATA(idx), 512);
    end

    % % Re-enable parallel processing for subsequent steps
    % pop_editoptions('option_parallel', 1);
end

% Assert
assert(all([DATA(:).srate] == 512));
fprintf('Done!\n');

end