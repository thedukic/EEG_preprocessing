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
num_block = length(N);

% Otherwise it does not make sense to run this function
assert(num_block > 1);

mask_block = false(num_block, sum(N));
for i_block = 1:num_block
    if i_block == 1
        mask_block(i_block, 1:N(1)) = true;
    else
        mask_block(i_block, sum(N(1:i_block-1)) + 1:sum(N(1:i_block))) = true;
    end
end

% Double-check
assert(sum(mask_block, "all") == sum(N));

% Mark eyes-open blocks
mask_eo = contains(EEG(1).ALSUTRECHT.subject.datablocks, 'EO');

% % Remove those that were removed completely due to very high noise
% if isfield(EEG(1).ALSUTRECHT, "bad_segments")
%     mask_extreme = EEG(1).ALSUTRECHT.bad_segments.maskRemoveblock;
%     assert(length(mask_eo) == length(mask_extreme));
%     mask_eo(mask_extreme) = [];
% end

% Log
for i_block = 1:num_block
    EEG(i_block).ALSUTRECHT.blocks.mask_eo    = mask_eo;
    EEG(i_block).ALSUTRECHT.blocks.mask_block = mask_block;
end

% Print summary info
fprintf('Total blocks: %d\n', num_block);
fprintf('Mean block length: %.0f samples (%.1f min, fs = %.0f Hz)\n', mean(N), (mean(N) / EEG(1).srate) / 60, EEG(1).srate);
fprintf('Done!\n');

end