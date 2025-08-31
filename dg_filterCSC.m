function [filteredSamples, T] = dg_filterCSC(TS, Samples, freqlim, varargin)
%filteredSamples = dg_filterCSC(TS, Samples)
% By default, performs the same filtering operation as lfp_bandpass2. Gain
% in the passband is unity (but note that this is only an average value for
% filter designs other than 'butter').  The filter is bandpass unless
% specified otherwise using OPTIONS.
%INPUTS
% TS - frame timestamps in microseconds, as returned by dg_readCSC.
% Samples - in samples X frames format, as returned by dg_readCSC.
% freqlim - in Hz.
%OUTPUTS
% filteredSamples - result of filtering <Samples>, same size as <Samples>.
% T - sample period.
%OPTIONS
% {'high' 'low' 'stop'} - any of these can be used to specify a type of
%   filter other than bandpass.  'high' and 'low' both require
%   a scalar value for <freqlim>.  'stop' produces an 8th-order bandstop
%   filter; see Matlab help for details.  Behavior is undefined if more
%   than one of these options is given, but as of this writing would be the
%   last one given. 
% {'cheby1' 'cheby2' 'ellip'} - any of these can be used to specify a
%   filter design other than Butterworth.
%NOTES
% For LFP processing, I recommend the default ('butter').  -DG

%$Rev: 169 $
%$Date: 2013-03-01 17:48:02 -0500 (Fri, 01 Mar 2013) $
%$Author: dgibson $

argnum = 1;
filtertype = 'bandpass';
filterdesign = 'butter';
R = 3;  % peak-to-peak ripple in dB for chebys and ellip
Rs = 10; % stopband attenuation in dB re: peak in passband for ellip
while argnum <= length(varargin)
    switch varargin{argnum}
        case {'high' 'low' 'stop'}
            filtertype = varargin{argnum};
        case {'cheby1' 'cheby2' 'ellip'}
            filterdesign = varargin{argnum};
        otherwise
            error('dg_filterCSC:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

T = 1e-6 * round(median(diff(TS))) / size(Samples,1);

% 'butter' requires freqs spec'd with 1.0 corresponding to half the sample rate. 
freqlim = freqlim*T*2;
designargs = {4};
switch filterdesign
    case {'cheby1' 'cheby2'}
        designargs{end+1} = R;
    case 'ellip'
        designargs{end+1} = R;
        designargs{end+1} = Rs;
end
designargs{end+1} = freqlim;
designargs{end+1} = filtertype;
[z, p, k] = feval(filterdesign, designargs{:});
[sos,g]=zp2sos(z,p,k);
h2=dfilt.df2sos(sos,g);
filteredSamples = filter(h2, Samples(end:-1:1));
filteredSamples = reshape( filter(h2, filteredSamples(end:-1:1)), ...
    size(Samples) );
