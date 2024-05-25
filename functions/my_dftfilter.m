function filt = my_dftfilter(dat,Fs,Fl,dftreplace,dftbandwidth,dftneighbourwidth)

% FT_PREPROC_DFTFILTER reduces power line noise (50 or 60Hz) using a
% discrete Fourier transform (DFT) filter, or spectrum interpolation.
%
% Use as
%   [filt] = my_dftfilter(dat, Fsample)
% or
%   [filt] = my_dftfilter(dat, Fsample, Fline)
% or
%   [filt] = my_dftfilter(dat, Fsample, Fline, 'dftreplace', 'zero')
% or
%   [filt] = my_dftfilter(dat, Fsample, Fline, 'dftreplace', 'neighbour')
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fline      frequency of the power line interference (if omitted from the input
%              the default European value of 50 Hz is assumed)
%
% Additional optional arguments are to be provided as key-value pairs:
%   dftreplace = 'zero' (default), 'neighbour' or 'neighbour_fft'.
%
% If dftreplace = 'zero', the powerline component's amplitude is estimated by
% fitting a sine and cosine at the specified frequency, and subsequently
% this fitted signal is subtracted from the data. The longer the sharper
% the spectral notch will be that is removed from the data.
% Preferably the data should have a length that is an integer multiple of the
% oscillation period of the line noise (i.e. 20ms for 50Hz noise). If the
% data is of different length, then only the first N complete periods are
% used to estimate the line noise. The estimate is subtracted from the
% complete data.
%
% If dftreplace = 'neighbour' the powerline component is reduced via
% spectrum interpolation (Leske & Dalal, 2019, NeuroImage 189,
% doi: 10.1016/j.neuroimage.2019.01.026), estimating the required signal
% components by fitting sines and cosines. The algorithmic steps are
% described in more detail below. % Preferably the data should have a length
% that is an integer multiple of the oscillation period of the line noise
% (i.e. 20ms for 50Hz noise). If the data is of different length, then only
% the first N complete periods are used to estimate the line noise. The
% estimate is subtracted from the complete data. Due to the possibility of
% using slightly truncated data for the estimation of the necessary signal
% components, this method is more forgiving with respect to numerical
% issues with respect to the sampling frequency, and suboptimal epoch
% lengths, in comparison to the method below.
%
% If dftreplace = 'neighbour_fft' the powerline component is reduced via spectrum
% interpolation, in its original implementation, based on an algorithm that
% uses an FFT and iFFT for the estimation of the spectral components. The signal is:
% I)   transformed into the frequency domain via a fast Fourier
%       transform (FFT),
% II)  the line noise component (e.g. 50Hz, dftbandwidth = 1 (±1Hz): 49-51Hz) is
%       interpolated in the amplitude spectrum by replacing the amplitude
%       of this frequency bin by the mean of the adjacent frequency bins
%       ('neighbours', e.g. 49Hz and 51Hz).
%       dftneighbourwidth defines frequencies considered for the mean (e.g.
%       dftneighbourwidth = 2 (±2Hz) implies 47-49 Hz and 51-53 Hz).
%       The original phase information of the noise frequency bin is
%       retained.
% III) the signal is transformed back into the time domain via inverse FFT
%       (iFFT).
% If Fline is a vector (e.g. [50 100 150]), harmonics are also considered.
% Preferably the data should be continuous or consist of long data segments
% (several seconds) to avoid edge effects. If the sampling rate and the
% data length are such, that a full cycle of the line noise and the harmonics
% fit in the data and if the line noise is stationary (e.g. no variations
% in amplitude or frequency), then spectrum interpolation can also be
% applied to short trials. But it should be used with caution and checked
% for edge effects.
%
% When using spectral interpolation, additional arguments are:
%
%   dftbandwidth      half bandwidth of line noise frequency bands, applies to spectrum interpolation, in Hz
%   dftneighbourwidth width of frequencies neighbouring line noise frequencies, applies to spectrum interpolation (dftreplace = 'neighbour'), in Hz
%
% If the data contains NaNs, the output of the affected channel(s) will be
% all(NaN).
%
% See also PREPROC

% Undocumented option:
%   Fline can be a vector, in which case the regression is done for all
%   frequencies in a single shot. Prerequisite is that the requested
%   frequencies all fit with an integer number of cycles in the data.
%
% Copyright (C) 2003, Pascal Fries
% Copyright (C) 2003-2015, Robert Oostenveld
% Copyright (C) 2016, Sabine Leske
% Copyright (C) 2021, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% defaults
% dftreplace        = 'zero';
% dftbandwidth      = [1 2 3];
% dftneighbourwidth = [2 2 2];

if isequal(dftreplace, 'zero')
    dftbandwidth      = NaN(size(Fl));
    dftneighbourwidth = NaN(size(Fl));
end

% determine the size of the data
[nchans, nsamples] = size(dat);

% set the default filter frequency
if nargin<3 || isempty(Fl)
    Fl = 50;
end

% ensure to be a column  vector
Fl = Fl(:);

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
    ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% determine the largest integer number of line-noise cycles that fits in
% the data, allowing for some numerical noise
n    = round(floor(nsamples .* (Fl./Fs + 100*eps)) .* Fs./Fl);

% check whether the filtering can be done in a single step, this can be
% done for the zero method if all n==n(1)
if (~strcmp(dftreplace, 'zero') && numel(n)>1) || ~all(n==n(1))
    % the different frequencies require different numbers of samples, apply the filters sequentially
    filt = dat;
    for i = 1:numel(Fl)
        filt = my_dftfilter(filt, Fs, Fl(i), 'dftreplace', dftreplace, 'dftbandwidth', dftbandwidth(i), 'dftneighbourwidth', dftneighbourwidth(i)); % enumerate all options
    end
    return
end

if strcmp(dftreplace,'zero')
    % Method A): DFT filter
    sel = 1:n(1);

    % temporarily remove mean to avoid leakage
    meandat = mean(dat(:,sel),2);
    dat = bsxfun(@minus, dat, meandat);

    % fit a sine and cosine to each channel in the data and subtract them
    time = (0:nsamples-1)/Fs;
    tmp  = exp(1i*2*pi*Fl*time);       % complex sin and cos
    ampl = 2*dat(:, sel)/tmp(:, sel);  % estimated amplitude of complex sin and cos on integer number of cycles
    est  = ampl*tmp;                   % estimated signal at this frequency
    filt = dat - est;                  % subtract estimated signal
    filt = real(filt);

    % add the mean back to the filtered data
    filt = bsxfun(@plus, filt, meandat);

elseif strcmp(dftreplace,'neighbour')
    dftbandwidth      = dftbandwidth(1:numel(Fl)); % this is a check due to the clunky defaults
    dftneighbourwidth = dftneighbourwidth(1:numel(Fl));

    % Method B1): DFT filter based estimation of stopband phase and DFT based
    % estimation of outside band power
    sel = 1:n;

    % temporarily remove the mean to avoid leakage
    meandat = mean(dat(:, sel),2);
    dat = bsxfun(@minus, dat, meandat);

    % fit a sine and cosine to the requested set of frequencies
    R    = Fs/n; % Rayleigh frequency
    time = (0:nsamples-1)/Fs;

    nbin = round((dftneighbourwidth+dftbandwidth)/R); % number of bins to be estimated at each side of the centre frequency
    freqs = (-nbin:nbin)*R + Fl;

    tmp  = exp(2*1i*pi*freqs(:)*time); % complex sin and cos
    beta = 2*dat(:, sel)/tmp(:, sel);  % estimated amplitude of complex sin and cos on integer number of cycles

    % boolean variable that indicates which frequency bins are to be replaced
    stopband = my_nearest(freqs - Fl, dftbandwidth.*[-1 1]);
    stopbool = false(1,numel(freqs));
    stopbool(stopband(1):stopband(2)) = true;

    % bandstop signal
    stopsignal = real(beta(:, stopbool)*tmp(stopbool, :));

    % retain the phase information
    beta(:, stopbool) = beta(:, stopbool)./abs(beta(:, stopbool));

    % estimate the amplitude from the flanking bands
    amp = mean(abs(beta(:, ~stopbool)),2);
    beta(:, stopbool) = beta(:, stopbool).*amp(:, ones(1, sum(stopbool)));

    % replacement signal
    replacesignal = real(beta(:, stopbool)*tmp(stopbool, :));

    filt = dat - stopsignal + replacesignal;

    % add the mean back to the filtered data
    filt = bsxfun(@plus, filt, meandat);

elseif strcmp(dftreplace,'neighbour_fft')
    % Method B): Spectrum Interpolation
    dftbandwidth      = dftbandwidth(:);
    dftneighbourwidth = dftneighbourwidth(:);
    if numel(Fl)<numel(dftbandwidth)
        dftbandwidth = dftbandwidth(1:numel(Fl));
    end
    if numel(Fl)<numel(dftneighbourwidth)
        dftneighbourwidth = dftneighbourwidth(1:numel(Fl));
    end

    % error message if periodicity of the interference frequency doesn't match the DFT length
    if n ~= nsamples
        error('Spectrum interpolation requires that the data length fits complete cycles of the powerline frequency, e.g., exact multiples of 20 ms for a 50 Hz line frequency (sampling rate of 1000 Hz).');
        % nfft = round(ceil(nsamples .* (Fl./Fs + 100*eps)) .* Fs./Fl);
        nfft = nsamples;
    else
        nfft = nsamples;
    end

    if (length(Fl) ~= length(dftbandwidth)) || (length(Fl) ~= length(dftneighbourwidth))
        ft_error('The number of frequencies to interpolate (cfg.dftfreq) should be the same as the number of bandwidths (cfg.dftbandwidth) and bandwidths of neighbours (cfg.neighbourwidth)');
    end

    % frequencies to interpolate
    f2int = [Fl(:)-dftbandwidth(:) Fl(:)+dftbandwidth(:)];

    % frequencies used for interpolation
    f4int = [f2int(:,1)-dftneighbourwidth(:) f2int f2int(:,2)+dftneighbourwidth(:)];

    data_fft = fft(dat, nfft, 2); % calculate fft to obtain spectrum that will be interpolated
    frq      = Fs*linspace(0, 1, nfft+1);

    % interpolate 50Hz (and harmonics) amplitude in spectrum
    for i = 1:length(Fl)
        smpl2int = nearest(frq,f2int(i,1)):nearest(frq,f2int(i,2));                                                   % samples of frequencies that will be interpolated
        smpl4int = [(nearest(frq,f4int(i,1)):nearest(frq,f4int(i,2))-1), (nearest(frq,f4int(i,3))+1:nearest(frq,f4int(i,4)))]; % samples of neighbouring frequencies used to calculate the mean

        % new amplitude is calculated as the mean of the neighbouring frequencies
        mns4int= bsxfun(@times, ones(size(data_fft(:,smpl2int))), mean(abs(data_fft(:,smpl4int)),2));

        % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
        data_fft(:, smpl2int) = bsxfun(@times, exp(bsxfun(@times,angle(data_fft(:, smpl2int)),1i)), mns4int);
    end

    % complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
    % to ensure a real valued signal after iFFT
    filt = ifft(data_fft,[],2,'symmetric');
    filt = filt(:, 1:nsamples);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indx] = my_nearest(array, val, insideflag, toleranceflag)

% NEAREST return the index of an array nearest to a scalar
%
% Use as
%   [indx] = nearest(array, val, insideflag, toleranceflag)
%
% The second input val can be a scalar, or a [minval maxval] vector for
% limits selection.
%
% If not specified or if left empty, the insideflag and the toleranceflag
% will default to false.
%
% The boolean insideflag can be used to specify whether the value should be
% within the array or not. For example nearest(1:10, -inf) will return 1,
% but nearest(1:10, -inf, true) will return an error because -inf is not
% within the array.
%
% The boolean toleranceflag is used when insideflag is true. It can be used
% to specify whether some tolerance should be allowed for values that are
% just outside the array. For example nearest(1:10, 0.99, true, false) will
% return an error, but nearest(1:10, 0.99, true, true) will return 1. The
% tolerance that is allowed is half the distance between the subsequent
% values in the array.
%
% See also FIND

% Copyright (C) 2002-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


if numel(val)==2
    % interpret this as a range specification like [minval maxval]
    % see also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1431
    intervaltol = eps;
    sel = find(array>=val(1) & array<=val(2));
    if isempty(sel)
        ft_error('The limits you selected are outside the range available in the data');
    end
    indx = sel([1 end]);
    if indx(1)>1 && abs(array(indx(1)-1)-val(1))<=intervaltol
        indx(1) = indx(1)-1;
    end
    if indx(2)<length(array) && abs(array(indx(2)+1)-val(2))<=intervaltol
        indx(2) = indx(2)+1;
    end
    return
end

mbscalar(val);

if nargin<3 || isempty(insideflag)
    insideflag = false;
end

if nargin<4 || isempty(toleranceflag)
    toleranceflag = false;
end

% ensure that it is a column vector
array = array(:);

% determine the most extreme values in the array
minarray = min(array);
maxarray = max(array);

% do some strict checks whether the value lies within the min-max range
if insideflag
    if ~toleranceflag
        if val<minarray || val>maxarray
            if numel(array)==1
                ft_warning('the selected value %g should be within the range of the array from %g to %g', val, minarray, maxarray);
            else
                ft_error('the selected value %g should be within the range of the array from %g to %g', val, minarray, maxarray);
            end
        end
    else
        if ~isequal(array, sort(array))
            ft_error('the input array should be sorted from small to large');
        end
        if numel(array)<2
            ft_error('the input array must have multiple elements to compute the tolerance');
        end
        mintolerance = (array(2)-array(1))/2;
        maxtolerance = (array(end)-array(end-1))/2;
        if val<(minarray-mintolerance) || val>(maxarray+maxtolerance)
            ft_error('the value %g should be within the range of the array from %g to %g with a tolerance of %g and %g on both sides', val, minarray, maxarray, mintolerance, maxtolerance);
        end
    end % toleragceflag
end % insideflag

% FIXME it would be possible to do some soft checks and potentially give a
% warning in case the user did not explicitly specify the inside and
% tolerance flags

% note that [dum, indx] = min([1 1 2]) will return indx=1
% and that  [dum, indx] = max([1 2 2]) will return indx=2
% whereas it is desired to have consistently the match that is most towards the side of the array

if val>maxarray
    % return the last occurrence of the largest number
    [dum, indx] = max(flipud(array));
    indx = numel(array) + 1 - indx;

elseif val<minarray
    % return the first occurrence of the smallest number
    [dum, indx] = min(array);

else
    % implements a threshold to correct for errors due to numerical precision
    % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=498 and http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1943
    %   if maxarray==minarray
    %     precision = 1;
    %   else
    %     precision = (maxarray-minarray) / 10^6;
    %   end
    %
    %   % return the first occurrence of the nearest number
    %   [dum, indx] = min(round((abs(array(:) - val)./precision)).*precision);

    % use find instead, see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1943
    wassorted = true;
    if ~issorted(array)
        wassorted = false;
        [array, xidx] = sort(array);
    end

    indx2 = find(array<=val, 1, 'last');
    indx3 = find(array>=val, 1, 'first');
    if abs(array(indx2)-val) <= abs(array(indx3)-val)
        indx = indx2;
    else
        indx = indx3;
    end
    if ~wassorted
        indx = xidx(indx);
    end
end