function EEG = load_biosemidata(subject,myPaths)

fprintf('\n================================\n');
fprintf('Loading data\n');
fprintf('================================\n');

% Load data
if strcmpi(myPaths.task,'MT')
    EEG = pop_biosig(subject.datablocks,'channels',1:168);
    EEG = pop_select(EEG,'rmchannel',[131 132 143:160]);
else
    EEG = pop_biosig(subject.datablocks,'channels',1:136);

    % Rare cases of RS datasets recorded with 168 channels
    if any([EEG(:).srate] == 2048)
        pop_editoptions('option_parallel', 0);

        maskfs2k = [EEG(:).srate] == 2048;
        EEG(maskfs2k) = pop_biosig(subject.datablocks(maskfs2k),'channels',1:168);
        EEG(maskfs2k) = pop_select(EEG(maskfs2k),'channel',[1:128 161:168]);
        EEG(maskfs2k) = pop_resample(EEG(maskfs2k),512);

        pop_editoptions('option_parallel', 1);
        assert(all([EEG(:).srate] == 512));
    end
end

end
