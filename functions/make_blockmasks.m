function EEG = make_blockmasks(EEG)
% Make a mask for each block
% e.g.,
% 111 000 000
% 000 111 000
% 000 000 111

fprintf('\n================================\n');
fprintf('Making recording block masks\n');
fprintf('================================\n');

% Mark start/stop of each block
N = cellfun(@(x) size(x,2), {EEG.data});
NBLK = length(N);

% Otherwise it does not make sense to run this function
assert(NBLK > 1);

mask_rs = false(NBLK, sum(N));
for i_block = 1:NBLK
    if i_block == 1
        mask_rs(i_block, 1:N(1)) = true;
    else
        mask_rs(i_block, sum(N(1:i_block-1)) + 1:sum(N(1:i_block))) = true;
    end
end

% Double-check
assert(sum(mask_rs, "all") == sum(N));

% Mark eyes-open blocks
eo_mask = contains(EEG(1).ALSUTRECHT.subject.datablocks, 'EO');

% Remove those that were removed completely due to very high noise
if isfield(EEG(1).ALSUTRECHT, "extremeNoise")
    maskRemoveblock = EEG(1).ALSUTRECHT.extremeNoise.maskRemoveblock;
    assert(length(eo_mask) == length(maskRemoveblock));
    eo_mask(maskRemoveblock) = [];
end

% Log
for i_block = 1:NBLK
    EEG(i_block).ALSUTRECHT.blockinfo.eo_mask = eo_mask;
    EEG(i_block).ALSUTRECHT.blockinfo.rs_mask = mask_rs;
end

% Print summary info
fprintf('Total blocks: %d\n', NBLK);
fprintf('Mean block length: %.0f samples (%.1f min, fs = %.0f Hz)\n', mean(N), (mean(N) / EEG(1).srate) / 60, EEG(1).srate);
fprintf('Done!\n');

end