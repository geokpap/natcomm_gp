function [C,phi,t,f,confC,phierr,Cerr]=dg_cohgramc(data1,data2,movingwin,tapers,pad,Fs,fpass,err,trialave,rmdc)
% Multi-taper time-frequency coherence - continuous processes
%
% Usage:
%
% [C,phi,t,f,confC,phierr,Cerr]=cohgramc(data1,data2,movingwin,tapers,nfft,Fs,fpass,err,trialave,rmdc)
% Input: 
% Note units have to be consistent. Thus, if movingwin is in seconds, Fs
% has to be in Hz. see chronux.m for more information.
%
%       data1 (in form samples x channels/trials) -- required
%       data2 (in form samples x channels/trials) -- required
%       movingwin (in the form [window winstep] -- required
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
%       C (abs of coherency, windows x frequency index x channels/trials)
%       phi (phase of coherency, windows x  frequency x channels/trials)
%       t (time in inverse units of Fs, row vector; relative to the sample
%           before the start of data1) 
%       f (frequencies, row vector)
%       confC (significance thresh for c at (1-p) level, scalar)
%       phierr (error bars for phi)
%       Cerr  (Jackknife confidence limits for C - use only for Jackknife:
%           {two rows, one for upper and one for lower} x windows x
%           frequency index x channels/trials)

% Copied from
% chronux\chronux_code\exploratory\correlations\continuous\cohgramc.m
% DG added <rmdc> param 2/28/2005
% Negative phase means that the second channel is lagging the first.

%$Rev: 327 $
%$Date: 2014-05-29 17:22:03 -0400 (Thu, 29 May 2014) $
%$Author: dgibson $

if nargin < 3; error('Need data1 and data2 and window parameters'); end;
[N1,C1]=size(data1);[N2,C2]=size(data2);
if N1~=N2 | C1~=C2; error('data incompatible'); end;
if nargin < 4; tapers=[3 5]; end;
if nargin < 5;pad=0;end;
if nargin < 6; Fs=1; end;
if nargin < 7; fpass=[0 Fs/2]; end;
if nargin < 8; err=0; end;
if nargin < 9; trialave=0;end;
if nargout > 6 & err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;

if isempty(tapers); tapers=[3 5]; end;
if isempty(pad);pad=0;end;
if isempty(Fs); Fs=1; end;
if isempty(fpass); fpass=[0 Fs/2]; end;
if isempty(err); err=0; end;
if isempty(trialave); trialave=0;end;


Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=2^(nextpow2(Nwin)+pad);
[f,findx]=chronux_getfgrid(Fs,nfft,fpass); 
tapers=chronux_dpsschk(tapers,Nwin)/sqrt(Fs); % check tapers

winstart=[1:Nstep:N1-Nwin+1];
nw=length(winstart);
for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin1=data1(indx,:);datawin2=data2(indx,:);
   if rmdc
       datawin1 = datawin1 - repmat(mean(datawin1, 1), size(datawin1, 1), 1);
       datawin2 = datawin2 - repmat(mean(datawin2, 1), size(datawin2, 1), 1);
   end
   if nargout==7;
     [c,ph,f,confc,phie,cerr]=chronux_coherencyc(datawin1,datawin2,tapers,pad,Fs,fpass,err,trialave);
     confC=confc;
     phierr(1,n,:,:)=squeeze(phie(1,:,:));
     phierr(2,n,:,:)=squeeze(phie(2,:,:));
     Cerr(1,n,:,:)=squeeze(cerr(1,:,:));
     Cerr(2,n,:,:)=squeeze(cerr(2,:,:));
   elseif nargout==6;
     [c,ph,f,confc,phie]=chronux_coherencyc(datawin1,datawin2,tapers,pad,Fs,fpass,err,trialave);
     confC=confc;
     phierr(1,n,:,:)=squeeze(phie(1,:,:));
     phierr(2,n,:,:)=squeeze(phie(2,:,:));
   else
     [c,ph,f]=chronux_coherencyc(datawin1,datawin2,tapers,pad,Fs,fpass,err,trialave);
   end;
   C(n,:,:)=c;
   phi(n,:,:)=ph;
end;
C=squeeze(C); phi=squeeze(phi);if nargout==7;Cerr=squeeze(Cerr);end;
if nargout==6; phierr=squeeze(phierr);end
winmid=winstart+round(Nwin/2);
t=winmid/Fs;