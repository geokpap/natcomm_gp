function [S,f,Serr]=chronux_mtspectrumc(data,tapers,pad,Fs,fpass,err,trialave)
% Multi-taper spectrum - continuous process
%
% Usage:
%
% [S,f,Serr]=chronux_mtspectrumc(data,tapers,pad,Fs,fpass,err,trialave)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (in form samples x channels/trials) -- required
%       tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	    pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%       Fs   (sampling frequency) - optional. Default 1.
%       fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%       err  (error calculation [1 p] - Theoretical error bars; [2 p] Jackknife error bars, 
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%       trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Output:
%       S       (spectrum in form frequency x channels/trials)
%       f       (frequencies)
%       Serr    (error bars)

if nargin < 1; error('Need data'); end;
if nargin < 2; tapers=[3 5]; end;
if nargin < 3;pad=0;end;
if nargin < 4; Fs=1; end;
if nargin < 5; fpass=[0 Fs/2]; end;
if nargin < 6; err=0; end;
if nargin < 7; trialave=0; end;
if isempty(tapers); tapers=[3 5]; end;
if isempty(pad);pad=0;end;
if isempty(Fs); Fs=1; end;
if isempty(fpass); fpass=[0 Fs/2]; end;
if isempty(err); err=0; end;
if isempty(trialave); trialave=0;end;

[N,C]=size(data);
nfft=2^(nextpow2(N)+pad);
[f,findx]=chronux_getfgrid(Fs,nfft,fpass); 
tapers=chronux_dpsschk(tapers,N)/sqrt(Fs); % check tapers
J=chronux_mtfftc(data,tapers,nfft);
J=J(findx,:,:);
S=squeeze(mean(conj(J).*J,2));
if trialave; S=squeeze(mean(S,2));end;
if nargout==3; 
   Serr=specerr(S,J,err,trialave);
end;