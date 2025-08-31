function [result, units] = lfp_bandpassPower(filenum, freqlim, ...
    smoothing, varargin)
%result = lfp_bandpassPower(filenum, freqlim, smoothing)
%   Function intended for use with lfp_createWave.
%   Applies lfp_bandpass2 to the wave in <filenum>, then squares each
%   sample, then smooths as in lfp_wavesmooth with parameter <smoothing>
%   and gain compensation (lfp_wavesmooth's default mode).
%OPTIONS
% All options except those listed below are passed through verbatim to
% lfp_bandpass2, which will raise an error if an option is not recognized.
% 'window', windowshape - allows you to select different shapes of
%   smoothing windows.  <windowshape> is a string, which can be 'hann'
%   (the default) for a hanning window, or 'box' for a rectangular window.

%$Rev: 359 $
%$Date: 2015-07-20 18:43:50 -0400 (Mon, 20 Jul 2015) $
%$Author: dgibson $

global lfp_SamplesUnits

opts2delete = [];
windowshape = 'hann';

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'window'
            opts2delete(end+1) = argnum; %#ok<AGROW>
            argnum = argnum + 1;
            opts2delete(end+1) = argnum; %#ok<AGROW>
            windowshape = varargin{argnum};
    end
end
varargin(opts2delete) = [];

result = lfp_bandpass2(filenum, freqlim, varargin{:});
result = result.^2;
switch windowshape
    case 'hann'
        winfunc = hanning(2*smoothing+1);
    case 'box'
        winfunc = ones(2*smoothing+1,1);
    otherwise
        error('lfp_bandpassPower:winshape', ...
            'Unrecognized window shape: "%s"', windowshape);
end
s1 = conv(result, winfunc);
% s1 has an "extra" <smoothing> points on each end, so the first
% "valid" point is (smoothing + 1) and the last is (end - smoothing).
result = s1(smoothing + 1 : end - smoothing);
if smoothing > 0
    result = result / sum(winfunc);
end
units = sprintf('%s^2', lfp_SamplesUnits{filenum});

