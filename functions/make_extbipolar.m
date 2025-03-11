function EEG = make_extbipolar(EEG)
% Put all external electrodes as type EXT?

bipolarEMG = false;

fprintf('\n================================\n');
fprintf('Making EOG/ECG channels bipolar\n');
fprintf('================================\n');

if strcmpi(EEG(1).ALSUTRECHT.subject.task,'MT')
    if bipolarEMG
        fprintf('Making EMG channels bipolar.\n');
    else
        fprintf('Keeping EMG channels monopolar.\n');
    end
end

NBLK = length(EEG);
for i = 1:NBLK
    % Bipolar ECG (or monopolar EL)
    % Early dataset had earlobes instead of ECG
    eleclabels = {EEG(i).chanlocs.labels};
    m1 = ismember(eleclabels,'ECGL');
    m2 = ismember(eleclabels,'ECGR');

    if any(m1) && any(m2)
        % Bipolar ECG
        EEG(i).chanlocs(m1).labels = 'ECG';
        % EEG(i).chanlocs(m1).type   = 'ECG';

        EEG(i).data(m1,:)          = EEG(i).data(m1,:)-EEG(i).data(m2,:);
        EEG(i).data(m2,:)          = [];
        EEG(i).chanlocs(m2)        = [];
    else
        % Monopolar earlobes
        % m1 = ismember(eleclabels,'LEL');
        % m2 = ismember(eleclabels,'REL');
        % EEG(i).chanlocs(m1).type   = 'Earlobe';
        % EEG(i).chanlocs(m2).type   = 'Earlobe';
    end

    % Bipolar VEOG
    eleclabels = {EEG(i).chanlocs.labels};
    m1 = ismember(eleclabels,'VEOGS');
    m2 = ismember(eleclabels,'VEOGI');
    EEG(i).chanlocs(m1).labels = 'VEOG';
    % EEG(i).chanlocs(m1).type   = 'EOG';
    EEG(i).data(m1,:)          = EEG(i).data(m1,:)-EEG(i).data(m2,:);
    EEG(i).data(m2,:)          = [];
    EEG(i).chanlocs(m2)        = [];

    % Bipolar HEOG
    eleclabels = {EEG(i).chanlocs.labels};
    m1 = ismember(eleclabels,'HEOGL');
    m2 = ismember(eleclabels,'HEOGR');
    EEG(i).chanlocs(m1).labels = 'HEOG';
    % EEG(i).chanlocs(m1).type   = 'EOG';
    EEG(i).data(m1,:)          = EEG(i).data(m1,:)-EEG(i).data(m2,:);
    EEG(i).data(m2,:)          = [];
    EEG(i).chanlocs(m2)        = [];

    % % Bipolar mastoid
    % eleclabels = {EEG(i).chanlocs.labels};
    % m1 = ismember(eleclabels,'LM');
    % m2 = ismember(eleclabels,'RM');
    % EEG(i).chanlocs(m1).labels = 'M';
    % EEG(i).chanlocs(m1).type   = 'M';
    % EEG(i).data(m1,:)          = EEG(i).data(m1,:)-EEG(i).data(m2,:);
    % EEG(i).data(m2,:)          = [];
    % EEG(i).chanlocs(m2)        = [];

    % % Left/Right mastoid
    % eleclabels = {EEG(i).chanlocs.labels};
    % m1 = ismember(eleclabels,'LM');
    % m2 = ismember(eleclabels,'RM');
    % if any(m1)
    %     EEG(i).chanlocs(m1).type   = 'Mastoid';
    % end
    % if any(m2)
    %     EEG(i).chanlocs(m2).type   = 'Mastoid';
    % end

    % EMG
    if strcmpi(EEG(1).ALSUTRECHT.subject.task,'MT')
        emglabels = {'APB','FDI','FPB','EPB','EDC','FDS'};
        emgindx   = find(ismember({EEG(i).chanlocs.type},'EMG'));
        if bipolarEMG
            cnt = 0;
            for j = emgindx(1):2:emgindx(end-1)
                cnt = cnt+1;
                EEG(i).data(j,:,:)           = EEG(i).data(j,:,:)-EEG(i).data(j+1,:,:);
                EEG(i).chanlocs(j).labels    = emglabels{cnt};
            end

            emgrmv = emgindx(2):2:emgindx(end);
            EEG(i).data(emgrmv,:)   = [];
            EEG(i).chanlocs(emgrmv) = [];
        else
            cnt = 0;
            for j = emgindx(1):2:emgindx(end-1)
                cnt = cnt+1;
                EEG(i).chanlocs(j).labels   = emglabels{cnt};
                EEG(i).chanlocs(j+1).labels = emglabels{cnt};
            end
        end
    end

    % Finalise
    EEG(i).nbchan = size(EEG(i).data,1);
    EEG(i) = eeg_checkset(EEG(i),'loaddata');

    % assert(size(EEG(i).data,1) == EEG(i).nbchan);
end

fprintf('Done!\n');

end