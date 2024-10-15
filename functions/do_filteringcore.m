function data = do_filteringcore(b,a,data,eventStruct,srate)

assert(isstable(b,a));
FNYQ = srate/2;

eventLabels = {eventStruct(:).type};
eventTimes  = [eventStruct(:).latency];
mask = contains(eventLabels,'boundary');

jumpsBad = floor([0 eventTimes(mask) size(data,2)]);
jumpsBad = unique(jumpsBad);

NChunk = length(jumpsBad)-1;
if NChunk>1, fprintf('Boundaries detected! Each chunk (N = %d) will be filtered separately.\n',NChunk); end

for j = 1:NChunk
    % Indices
    dataInd = [jumpsBad(j)+1 jumpsBad(j+1)];
    dataInd = dataInd(1):dataInd(2);

    if length(dataInd)>srate
        dataTmp = double(data(:,dataInd));

        % Remove DC (big offsets can cause artifacts)
        dataTmp = remove_dcsignal(dataTmp, FNYQ);

        % Filter
        dataTmp = filtfilt(b,a, dataTmp')';

        % Place back
        data(:,dataInd) = dataTmp;
    else
        warning('Chunk %d: Data too small (<=%d) for filtering (L = %d samples)!',j,srate,length(dataInd));
    end
end

% fprintf('Done!\n');

end