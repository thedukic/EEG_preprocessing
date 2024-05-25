function EEG = make_rsmasks(EEG)
% Make a mask for each trial
% eg
% 111 000 000
% 000 111 000
% 000 000 111
% N = cellfun(@(x,y) size(x,2),{EEG(:).data});

% Mark where each block starts/stops
N = cellfun(@(x,y) size(x,2),{EEG(:).data});
NBLK = length(N);

rs_mask = false(NBLK,sum(N));
for j = 1:NBLK
    if j == 1
        rs_mask(j,1:N(1)) = true;
    else
        rs_mask(j,sum(N(1:j-1))+1:sum(N(1:j))) = true;
    end
end

% Check
assert(sum(rs_mask,"all")==sum(N));

% Which are eyes-open blocks
eo_mask = contains(EEG(1).ALSUTRECHT.subject.datablocks,'EO');

% Log
for i = 1:NBLK
    EEG(i).ALSUTRECHT.blockinfo.eo_mask = eo_mask;
    EEG(i).ALSUTRECHT.blockinfo.rs_mask = rs_mask;
end

end