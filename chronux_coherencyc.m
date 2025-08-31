function [C,phi,f,confC,phierr,Cerr]=chronux_coherencyc(data1,data2,tapers,pad,Fs,fpass,err,trialave)
% Multi-taper coherency - continuous process
%
% Usage:
% [C,phi,f,confC,phierr,Cerr]=coherencyc(data1,data2,tapers,nfft,Fs,fpass,err,trialave)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data1 (in form samples x channels/trials) -- required
%       data2 (in form samples x channels/trials) -- required
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
%       err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%       trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Output:
%       C (abs of coherency frequency index x channels/trials)
%       phi (phase of coherency frequency x channels/trials)
%       f (frequencies)
%       confC (confidence level for c at 1-p %)
%       phierr (error bars for phi)
%       Cerr  (Jackknife error bars for C - use only for Jackknife)
if nargin < 2; error('Need data1 and data2'); end;warning off MATLAB:divideByZero
[N,dum]=size(data1);
if nargin < 3; tapers=[3 5]; end;
if nargin < 4;pad=0;end;
if nargin < 5; Fs=1; end;
if nargin < 6; fpass=[0 Fs/2]; end;
if nargin<7; err=0; end;
if nargin<8; trialave=0;end;
if nargout > 5 & err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;

if isempty(tapers); tapers=[3 5]; end;
if isempty(pad);pad=0;end;
if isempty(Fs); Fs=1; end;
if isempty(fpass); fpass=[0 Fs/2]; end;
if isempty(err); err=0; end;
if isempty(trialave); trialave=0;end;

nfft=2^(nextpow2(N)+pad);
[f,findx]=chronux_getfgrid(Fs,nfft,fpass); 
tapers=chronux_dpsschk(tapers,N); % check tapers
J1=chronux_mtfftc(data1,tapers,nfft);
J2=chronux_mtfftc(data2,tapers,nfft);
J1=J1(findx,:,:); J2=J2(findx,:,:);
S12=squeeze(mean(conj(J1).*J2,2));
S1=squeeze(mean(conj(J1).*J1,2));
S2=squeeze(mean(conj(J2).*J2,2));
if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
C12=S12./sqrt(S1.*S2);
%if trialave; C12=squeeze(mean(C12,2)); end;
C=abs(C12); 
phi=angle(C12);
if nargout==6; 
     [confC,phierr,Cerr]=coherr(C,J1,J2,err,trialave);
elseif nargout==5;
     [confC,phierr]=coherr(C,J1,J2,err,trialave);
end;
