function [S,winmid,f,Serr]=lfp_mtspecgram2(data,movingwin,tapers,pad,Fs,...
    fpass,err,trialave,rmdc,norm)
% Multi-taper time-frequency spectrum - continuous process
% Modified March 2005 by Dan Gibson from version downloaded from
% chronux.org dated 8/20/2004.  Added <rmdc> param, which is required.
% June 2005 Dan Gibson added <norm> param, which is required.
% 5/17/2007: DG removed code assigning default values to params that are no
% longer optional, and replaced return value <t> with <winmid>, where
% t=winmid/Fs.  <winmid> is a list of sample numbers corresponding to the
% middle of each window.
% 1/15/2008: DG removed normalization of tapers by sqrt(Fs) because this
% version of chronux disagrees with itself as to where and how many times
% that should be done, and chronux 1.50 moved the entire issue into chronux_dpsschk
% itself, which appears to do the opposite normalization (but that code is
% not compatible here anyway).
%
% Usage:
% [S,winmid,f,Serr]=mtspectgramc(data,movingwin,tapers,pad,Fs,fpass,err,trialave)
% Input: 
% Note units have to be consistent. Thus, if movingwin is in seconds, Fs
% has to be in Hz. see chronux.m for more information.
%       data        (in form samples x channels/trials) -- required
%       movingwin         (in the form [window winstep] i.e length of moving
%                                                 window and step size)
%                                                 Note that units here have
%                                                 to be consistent with
%                                                 units of Fs
%       tapers      (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	    pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	   e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	   to 512 points; if PAD = 2, we pad the FFT
%			      	   to 2048 points, etc.
%       Fs          (sampling frequency) - optional. Default 1.
%       fpass       (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%       err         (error calculation [1 p] - Theoretical error bars; [2 p] Jackknife error bars,
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%       trialave    (average over trials when 1, don't average when 0) - optional. Default 0
%       norm        norm > 0: normalize each time window to the average
%                   power over norm <= f <= Fs/2.  norm = 0: no
%                   normalization.  norm = -1: normalize each time window
%                   to sum of samples in that window, computed before
%                   applying 'rmdc'.
% Output:
%       S       (spectrum in form time x frequency x channels/trials)
%       winmid  sample indices of window midpoints [-DG]
%       f       (frequencies)
%       Serr    (error bars)
%

%$Rev: 327 $
%$Date: 2014-05-29 17:22:03 -0400 (Thu, 29 May 2014) $
%$Author: dgibson $

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; tapers=[3 5]; end;
if nargin < 4;pad=0;end;
if nargin < 5; Fs=1; end;
if nargin < 6; fpass=[0 Fs/2]; end;
if nargin < 7; err=0; end;
if nargin < 8; trialave=0; end;
if isempty(tapers); tapers=[3 5]; end;
if isempty(pad);pad=0;end;
if isempty(Fs); Fs=1; end;
if isempty(fpass); fpass=[0 Fs/2]; end;
if isempty(err); err=0; end;
if isempty(trialave); trialave=0;end;


[N,C]=size(data);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=2^(nextpow2(Nwin)+pad);
[f,findx]=chronux_getfgrid(Fs,nfft,fpass); 
tapers=chronux_dpsschk(tapers,Nwin); % check tapers

winstart=[1:Nstep:N-Nwin+1];
nw=length(winstart);
% Just to get the correct size of the padded s:
[s,f]=chronux_mtspectrumc(data(winstart(1):winstart(1)+Nwin-1,:),tapers,pad,Fs,fpass,0,trialave);
if norm > 0
    savefpass = fpass;
    fpass(2) = Fs/2;
end
% Pre-allocate S for speed:
S = zeros(nw, size(s, 1), size(s, 2));
for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   if norm == -1
       winsum = sum(datawin, 1);
       if trialave
           winsum = sum(winsum,2);
       end
   end
   if rmdc
       datawin = datawin - repmat(mean(datawin, 1), size(datawin, 1), 1);
   end
   if nargout==4;
     [s,f,serr]=chronux_mtspectrumc(datawin,tapers,pad,Fs,fpass,err,trialave);
     Serr(1,n,:,:)=squeeze(serr(1,:,:));
     Serr(2,n,:,:)=squeeze(serr(2,:,:));
   else;
     [s,f]=chronux_mtspectrumc(datawin,tapers,pad,Fs,fpass,err,trialave);
   end;
   if norm == -1
       s = s./repmat(winsum, size(s,1), 1);
   elseif norm > 0
       hfmean = mean(s(f>=norm & f<=Fs/2,:), 1);
       s = s./repmat(hfmean, size(s,1), 1);
       s = s(f >= savefpass(1) & f <= savefpass(2), :);
   end
   S(n,:,:)=s;
end;
S=squeeze(S); if nargout==4;Serr=squeeze(Serr);end;
winmid=winstart+round(Nwin/2);
if norm > 0
    f = f(f >= savefpass(1) & f <= savefpass(2));
end
