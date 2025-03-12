function EEG = make_rsmasks(EEG)
% Make a mask for each block
% e.g.,
% 111 000 000
% 000 111 000
% 000 000 111

fprintf('\n================================\n');
fprintf('Making resting-state block masks\n');
fprintf('================================\n');

% Mark start/stop of each block
N = cellfun(@(x,y) size(x,2),{EEG(:).data});
NBLK = length(N);

rs_mask = false(NBLK,sum(N));
for i = 1:NBLK
    if i == 1
        rs_mask(i,1:N(1)) = true;
    else
        rs_mask(i,sum(N(1:i-1))+1:sum(N(1:i))) = true;
    end
end

% Double-check
assert(sum(rs_mask,"all") == sum(N));

% Mark eyes-open blocks
eo_mask = contains(EEG(1).ALSUTRECHT.subject.datablocks,'EO');

% Remove those that were removed completely due to very high noise
maskRemoveblock = EEG(1).ALSUTRECHT.extremeNoise.maskRemoveblock;
assert(length(eo_mask) == length(maskRemoveblock));
eo_mask(maskRemoveblock) = [];

% Log
for i = 1:NBLK
    EEG(i).ALSUTRECHT.blockinfo.eo_mask = eo_mask;
    EEG(i).ALSUTRECHT.blockinfo.rs_mask = rs_mask;
end

fprintf('Done!\n');

end