function EEG = do_reref(EEG,typeMethod)

fprintf('\n================================\n');
fprintf('Rereferencing data (%s)\n',typeMethod);
fprintf('================================\n');

chaneeg = strcmp({EEG(1).chanlocs.type},'EEG');
NBLK = length(EEG);

switch typeMethod
    case 'aRobust'
        for i = 1:NBLK
            NDIM = ndims(EEG(i).data);
            if NDIM == 2
                EEG(i).data(chaneeg,:) = EEG(i).data(chaneeg,:) - trimmean(EEG(i).data(chaneeg,:),20,'round',1);
            elseif NDIM == 3
                EEG(i).data(chaneeg,:,:) = EEG(i).data(chaneeg,:,:) - trimmean(EEG(i).data(chaneeg,:,:),20,'round',1);
            end
            EEG(i).ref = 'average';
        end

    case 'aRegular'
        for i = 1:NBLK
            NDIM = ndims(EEG(i).data);
            if NDIM == 2
                EEG(i).data(chaneeg,:) = EEG(i).data(chaneeg,:) - mean(EEG(i).data(chaneeg,:),1);
            elseif NDIM == 3
                EEG(i).data(chaneeg,:,:) = EEG(i).data(chaneeg,:,:) - mean(EEG(i).data(chaneeg,:,:),1);
            end
            EEG(i).ref = 'average';
        end

    case 'aRegular2'
        error('Method not supported!');
        % EEG.data = EEG.data - (sum(EEG.data,1)/(EEG.nbchan+1));
        % EEG.ref  = 'average';
end

fprintf('Done!\n');

end