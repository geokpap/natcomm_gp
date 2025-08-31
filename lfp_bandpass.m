function [result, units] = lfp_bandpass(filenum, freqlim)
%result = lfp_bandpass(filenum, freqlim)
% WARNING: YOU PROBABLY WANT TO USE lfp_bandpass2 INSTEAD OF THIS!
%
% Creates a 4th-order Butterworth filter with cutoff freqs <freqlim>,
% applies it via filtfilt in forwards and backwards time to the waveform in
% <filenum>, and returns result.  The phase difference is thus zero at all
% frequencies and filtering becomes 8th order.  A lowpass filter can be
% specified by giving just the upper cutoff frequency.  Do NOT specify 0
% as the lower cutoff frequency for this purpose!
%
%   WARNING: for certain values of B that are a bit too close to zero,
% Matlab's filtfilt function can still return all zeros even when B is
% not zero, if its input is greater than a certain length.  Therefore,
% if this function bombs, tweak the frequency limits by inconsequential
% amounts.

%   NOTE:  the Matlab BUTTER function can return a filter design with
% all zeros in the numerator, resulting in an "all-stop filter", for certain
% magical combinations of values.  One such combination is 
%   [5 10]*lfp_SamplePeriod*2
% when lfp_SamplePeriod is approximately 0.000501666015624735 (there seems 
% to be some roundoff error in the decimal conversion, so if you try it 
% with this decimal value, you will not see the problem).  Tweaking it ever
% so slightly fixes it: e.g. multiplying by 1+eps or dividing by 1+2*eps 
% (dividing by 1+eps makes no difference).

%$Rev: 251 $
%$Date: 2011-12-20 19:55:13 -0500 (Tue, 20 Dec 2011) $
%$Author: dgibson $

if ~isequal(size(filenum), [1 1])
    error('lfp_bandpass:badfilenum', ...
        '<filenum> must contain exactly one value');
elseif ~strcmp(class(filenum), 'double')
    error('lfp_bandpass:badfilenum2', ...
        '<filenum> must be a number');
end
% if ~isequal(size(freqlim), [1 2])
%     error('lfp_bandpass:badfreqlim', ...
%         '<freqlim> must be row vector of two values');
% elseif ~strcmp(class(freqlim), 'double')
%     error('lfp_bandpass:badfreqlim2', ...
%         '<freqlim> must be numeric');
% end

lfp_declareGlobals;

% 'butter' requires freqs spec'd with 1.0 corresponding to half the sample rate. 
freqlim = freqlim*lfp_SamplePeriod*2;
[B,A] = butter(4, freqlim);
if ~any(B)
    OK = false;
    for k = 1:10 
        [B,A] = butter(4, freqlim*(1+k*eps));
        if any(B)
            OK = true;
            break;
        end
        if ~OK
            error('lfp_bandpass:arghh', ...
            'Incomprehensible problems with BUTTER' );
        end
    end
end
result = filtfilt(B, A, reshape(lfp_Samples{filenum}, 1, []));
units = lfp_SamplesUnits{filenum};
