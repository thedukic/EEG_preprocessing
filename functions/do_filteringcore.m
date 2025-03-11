function data = do_filteringcore(b,a,data,eventStruct,srate)

% Initialise
assert(isstable(b,a));
data = double(data);
FNYQ = srate/2;

if ~isempty(eventStruct)
    eventLabels = {eventStruct(:).type};
    eventTimes  = [eventStruct(:).latency];
    isBoundary  = contains(eventLabels,'boundary'); % disp(isBoundary);

    % Avoid the previous errors in renaming the event names
    % -> this has happened rarly if there were CMS drop-outs
    assert(length(eventLabels) == length(eventTimes));
    assert(all(strcmpi(eventLabels(isBoundary),'boundary')));

    jumpsBad = floor([0 eventTimes(isBoundary) size(data,2)]);
    jumpsBad = unique(jumpsBad);
else

    jumpsBad = [0 size(data,2)];
end

NChunk = length(jumpsBad)-1;
if NChunk>1
    fprintf('Boundaries detected! Each chunk (N = %d) will be filtered separately.\n',NChunk);
else
    fprintf('Great! Boundaries are not detected. The whole recording will be filtered at once.\n');
end

for i = 1:NChunk
    % Indices
    dataInd = [jumpsBad(i)+1 jumpsBad(i+1)];
    dataInd = dataInd(1):dataInd(2);

    if length(dataInd) > srate
        % Remove DC (big offsets can cause artifacts)
        dataTmp = remove_dcsignal(data(:,dataInd), FNYQ);

        % Filter and place back
        data(:,dataInd) = filtfilt(b,a,dataTmp')';
    else
        warning('Chunk %d: Data too small (<=%d) for filtering (L = %d samples)!',i,srate,length(dataInd));
    end
end

% fprintf('Done!\n');

end