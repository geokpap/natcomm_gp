function [result, units] = lfp_bandpass2(filenum, freqlim, varargin)
%result = lfp_bandpass2(filenum, freqlim)
%
% Creates a 4th-order Butterworth filter with cutoff freqs <freqlim>,
% applies it in forwards and backwards time to the waveform in
% <filenum>, and returns result.  The phase difference is thus zero at all
% frequencies and filtering becomes 8th order.  A lowpass filter can be
% specified by giving just the upper cutoff frequency.  Do NOT specify 0
% as the lower cutoff frequency for this purpose!
%
%   NOTE:  identical to lfp_bandpass, except instead of passing a TF design
% to filtfilt (which seems to require one), I ignore the startup transients
% and just run it backwards and forwards through "filter".  Seems to behave
% itself better than lfp_bandpass.
%
%OPTIONS
% 'fwdonly' - do not do the reverse-time pass.
% {'high' 'low' 'stop'} - any of these can be used to specify a type of
%   Butterworth filter other than bandpass.  'high' and 'low' both require
%   a scalar value for <freqlim>.  'stop' produces an 8th-order bandstop
%   filter; see Matlab help for details.  Behavior is undefined if more
%   than one of these options is given, but as of this writing would be the
%   last one given. 
% 'order', N - Instead of creating a 4th-order filter, creates a filter of
%   order <N> (so after the double filtering, the order is 2*<N>).
% 'revonly' - do not do the forwards-time pass.

%$Rev: 392 $
%$Date: 2018-01-08 12:32:41 -0500 (Mon, 08 Jan 2018) $
%$Author: dgibson $

global lfp_SamplePeriod lfp_Samples lfp_SamplesUnits

if ~isequal(size(filenum), [1 1])
    error('lfp_bandpass:badfilenum', ...
        '<filenum> must contain exactly one value');
elseif ~strcmp(class(filenum), 'double')
    error('lfp_bandpass:badfilenum2', ...
        '<filenum> must be a number');
end

argnum = 1;
order = 4;
buttertype = '';
fwdflag = true;
revflag = true;
while argnum <= length(varargin)
    switch varargin{argnum}
        case {'fwdonly'}
            revflag = false;
        case {'high' 'low' 'stop'}
            buttertype = varargin{argnum};
        case 'order'
            argnum = argnum + 1;
            order = varargin{argnum};
        case {'revonly'}
            fwdflag = false;
        otherwise
            error('lfp_bandpass2:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% 'butter' requires freqs spec'd with 1.0 corresponding to half the sample rate. 
freqlim = freqlim*lfp_SamplePeriod*2;
if isempty(buttertype)
    [z, p, k] = butter(order, freqlim);
else
    [z, p, k] = butter(order, freqlim, buttertype);
end
[sos,g]=zp2sos(z,p,k);
h2=dfilt.df2sos(sos,g);
if revflag
    result = filter(h2, lfp_Samples{filenum}(end:-1:1));
    if ~fwdflag
        result = result(end:-1:1);
    end
end
if fwdflag
    if revflag
        result = filter(h2, result(end:-1:1));
    else
        result = filter(h2, lfp_Samples{filenum}(:));
    end
end
% Someday I need to figure out how to normalize the gain
%units = lfp_SamplesUnits{filenum};
% 4/15/09: it turns out that the passband gain is unity for 'low', 'high',
% 'stop', and bandpass.
units = lfp_SamplesUnits{filenum};
