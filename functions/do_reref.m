function EEG = do_reref(EEG,typeRef)

chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
NBLK = length(EEG);

switch typeRef
    case 'aRobust'
        for i = 1:NBLK
            EEG(i).data(chaneeg,:) = EEG(i).data(chaneeg,:) - trimmean(EEG(i).data(chaneeg,:),20,'round',1);
            EEG(i).ref = 'average';
        end

    case 'aRegular'
        for i = 1:NBLK
            EEG(i).data(chaneeg,:) = EEG(i).data(chaneeg,:) - mean(EEG(i).data(chaneeg,:),1);
            EEG(i).ref = 'average';
        end
    case 'aRegular2'
        error('Not supported!');
        % EEG.data = EEG.data - (sum(EEG.data,1)/(EEG.nbchan+1));
        % EEG.ref  = 'average';
end

end