function EEG = do_reref(EEG, typeMethod)

fprintf('\n================================\n');
fprintf('Rereferencing EEG electrodes (%s)\n', typeMethod);
fprintf('================================\n');

chaneeg = strcmp({EEG(1).chanlocs.type}, 'EEG');
num_channel = sum(chaneeg);
assert(num_channel <= 128);
num_block = length(EEG);
fprintf('Referencing %d EEG channels across %d blocks.\n', num_channel, num_block);

switch typeMethod
    case 'aRobust'
        for i_block = 1:num_block
            % Isolate the target EEG channels
            data_curr = EEG(i_block).data(chaneeg, :, :);
            [n_chan, n_samples, n_trials] = size(data_curr); % n_trials automatically returns 1 for 2D data

            % Flatten to a 2D matrix: [channels x (samples * trials)]
            data_2d = reshape(data_curr, n_chan, []);

            % Estimate vectorised Huber mean across channels (dimension 1)
            mu_2d = estimate_huber_mean(data_2d, 1.345);

            % Subtract the robust average reference
            data_ref_2d = data_2d - mu_2d;

            % Reshape back to original dimensions (handles both 2D and 3D perfectly)
            EEG(i_block).data(chaneeg, :, :) = reshape(data_ref_2d, n_chan, n_samples, n_trials);
            EEG(i_block).ref = 'average';
        end

    case 'aRobust_trim'
        for i_block = 1:num_block
            % Isolate the target EEG channels
            data_curr = EEG(i_block).data(chaneeg, :, :);
            [n_chan, n_samples, n_trials] = size(data_curr);

            % Flatten to a 2D matrix: [channels x (samples * trials)]
            data_2d = reshape(data_curr, n_chan, []);

            % Estimate trimmed mean across channels (dimension 1)
            mu_2d = trimmean(data_2d, 20, 'round', 1);

            % Subtract the robust trimmed average reference
            data_ref_2d = data_2d - mu_2d;

            % Reshape back to original dimensions (handles both 2D and 3D)
            EEG(i_block).data(chaneeg, :, :) = reshape(data_ref_2d, n_chan, n_samples, n_trials);
            EEG(i_block).ref = 'average';
        end

    case 'aRegular'
        for i_block = 1:num_block
            % Isolate the target EEG channels
            data_curr = EEG(i_block).data(chaneeg, :, :);
            [n_chan, n_samples, n_trials] = size(data_curr);

            % Flatten to a 2D matrix: [channels x (samples * trials)]
            data_2d = reshape(data_curr, n_chan, []);

            % Estimate standard mean across channels (dimension 1)
            mu_2d = mean(data_2d, 1);
            % mu_2d = sum(data_2d, 1) / (EEG(i_block).nbchan + 1);

            % Subtract the regular average reference
            data_ref_2d = data_2d - mu_2d;

            % Reshape back to original dimensions (handles both 2D and 3D)
            EEG(i_block).data(chaneeg, :, :) = reshape(data_ref_2d, n_chan, n_samples, n_trials);
            EEG(i_block).ref = 'average';
        end

    case 'rest'
        % Load leadfield
        load('leadfield_template.mat', 'L');

        for i_block = 1:num_block
            % Isolate the target EEG channels
            data_curr = EEG(i_block).data(chaneeg, :, :);
            [n_chan, n_samples, n_trials] = size(data_curr);

            % Flatten to a 2D matrix: [channels x (samples * trials)]
            data_2d = reshape(data_curr, n_chan, []);

            % Apply REST reference on the flattened 2D matrix
            data_rest_2d = rest_refer(data_2d, L);

            % Reshape back to original dimensions (handles both 2D and 3D)
            EEG(i_block).data(chaneeg, :, :) = reshape(data_rest_2d, n_chan, n_samples, n_trials);
        end

    otherwise
        error('!');

end

fprintf('Done!\n');

end

function [data_z] = rest_refer(data,G)
%   Main function of Reference Electrode Standardization Technique
%   Corrdinate system: Three-concentric-sphere head volume conductor
%       model. The triangle shows the nose. The centre of the spheres is
%       defined as the coordinate origin. The axis directed away from the
%       origin toward the left ear is defined as the -x axis, and that
%       from the origin to the nasion is the +y axis. The +z axis is defined
%       as the axis that is perpendicular to both these axes and directed
%       from the origin to the vertex.
%   Input:
%         data:  The EEG potentials with average reference,channels X time points,
%                e.g. 62 channels X 10000 time points.The original reference must be
%                average reference.
%         G: Lead Field matrix, channels X sources, e.g. 62 channels X 3000 sources.
%   Output:
%         data_z: The EEG potentials with zero reference,
%                channels X time points.
%
%  Author: Shiang Hu
%      Date: Sep. 22, 2016
%  Edit by Li Dong (Sep.24,2016)
%  Edit by Li Dong (Aug. 28, 2016)
%   For more see http://www.neuro.uestc.edu.cn/rest/
%   Reference: Yao D (2001) A method to standardize a reference of scalp EEG recordings to a point at infinity.
%                       Physiol Meas 22:693?11. doi: 10.1088/0967-3334/22/4/305

if nargin < 2
    error('Please input the Lead Field matrix!');
end
if size(data,1) ~= size(G,1)
    error('No. of Channels of lead field matrix and data are NOT equal!');
end
if size(data,1) > size(data,2)
    warning('No. of channels > No. of time points in data???');
end

% UNIT CONVERSION: Scale leadfield to match the 0.05 regularisation gate
% Compute the singular values of the demeaned leadfield
% Divide G by its maximum singular value to normalise it to a max scale of 1.0
leadfield_singular_values = svd(G - mean(G, 1));
G = G / max(leadfield_singular_values);

% Enforce a strict average reference across all channels
data = data - mean(data, 1);
Gar = G - mean(G, 1);

% The value 0.05 is for real data; for simulated data, it may be set as zero
data_z = G * pinv(Gar, 0.05) * data;

% V = V_avg + AVG(V_0)
data_z = data + repmat(mean(data_z), size(G,1), 1);

end