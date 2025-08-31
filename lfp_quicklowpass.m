function [result, units] = lfp_quicklowpass(filenum, freqlim, varargin)
%[result, units] = lfp_quicklowpass(filenum, freqlim)
% This is a simpler, faster replacement for lfp_bandpass2 when lowpass
% filtering is desired.  It is essentially just an lfp_createWave-friendly
% wrapper for dg_quicklowpass.

%$Rev: 384 $
%$Date: 2016-08-25 17:52:27 -0400 (Thu, 25 Aug 2016) $
%$Author: dgibson $

global lfp_SamplePeriod lfp_Samples lfp_SamplesUnits

if ~isequal(size(filenum), [1 1])
    error('lfp_bandpass:badfilenum', ...
        '<filenum> must contain exactly one value');
elseif ~isnumeric(filenum)
    error('lfp_bandpass:badfilenum2', ...
        '<filenum> must be a number');
end
if ~isequal(size(freqlim), [1 1])
    error('lfp_bandpass:badfilenum', ...
        '<freqlim> must contain exactly one value');
elseif ~isnumeric(freqlim)
    error('lfp_bandpass:badfilenum2', ...
        '<freqlim> must be a number');
end

result = dg_quicklowpass( ...
    freqlim, 1/lfp_SamplePeriod, lfp_Samples{filenum});
units = lfp_SamplesUnits{filenum};
