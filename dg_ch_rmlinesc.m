function data=dg_ch_rmlinesc(data,tapers, Fs, passband, pad, p)
% DEPRECATED.  See lfp_rmlnspec.
% Removes significant sine waves from data (continuous data).
% 3-Feb-2007: Hacked from 8/24/2004 version of rmlinesc in chronux.
%NOTES
% As far as I can tell, I abandoned this in favor of lfp_rmlnspec which
% both does more (removes harmonic series) and less (does not attempt to
% automatedly identify the frequencies to remove).  DG 11-Mar-2016.
% Also see lfp_rmlinesc (another deprecated function).
%
% Usage: data = dg_ch_rmlinesc(data,tapers,Fs,passband,pad,p)
%
%  Inputs:  
% Note that units of Fs, passband have to be consistent.
%       data        (data in [N,C] i.e. time x channels/trials) - required.
%	    tapers 	    (parameters for calculating tapers [NW,K]) - optional. Defaults to [3 5]
%	    Fs 	        (sampling frequency) -- optional. Defaults to 1.
%	    passband 	(band of frequencies to be kept [fmin fmax]) - optional. Defaults to [0 Fs/2]
%	    pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%	    p		    (P-value to calculate error bars for) - optional. Defaults to 0.05 (95% confidence).
%
%
%  Outputs: 
%       data        (data with significant lines removed)
%

%$Rev: 258 $
%$Date: 2017-06-21 17:36:42 -0400 (Wed, 21 Jun 2017) $
%$Author: dgibson $

[N,C]=size(data);
if nargin<2; tapers=[3 5]; end;
if nargin<3; Fs=1;end;
if nargin<4;passband=[0 Fs/2];end;
if nargin<5;pad=0;end;
if nargin<6;p=0.05;end;
tapers=dpsschk(tapers,N); % calculate the tapers
[Fval,A,f,sig,sd] = ftestc(data,tapers,Fs,passband,pad,p);
fmax=findpeaks(Fval,sig);
for ch=1:C;
    fsig=f(fmax(ch).loc);
    Nf=length(fsig);
    msg = sprintf('The significant lines for channel %d and the amplitudes are \n',ch);
    for nf=1:Nf;
        msg = [msg sprintf('f=%12.8f\n',fsig(nf))];
        msg = [msg sprintf('Re(A)=%12.8f\n',real(A(fmax(ch).loc,ch)))];
        msg = [msg sprintf('Im(A)=%12.8f\n',imag(A(fmax(ch).loc,ch)))];
    end;
    lfp_log(msg);
    datasine(:,ch)=exp(i*2*pi*[1:N]'*fsig/Fs)*A(fmax(ch).loc,ch)+exp(-i*2*pi*[1:N]'*fsig/Fs)*conj(A(fmax(ch).loc,ch));
end;
datan=data-datasine;
if nargout==0; 
   figure;subplot(211); plot(f,Fval); line(get(gca,'xlim'),[sig sig]);
   px1=pmtm(data(:,1));
   px2=pmtm(datan(:,1));
   subplot(212);plot(1:length(px1),10*log10(px1),1:length(px2),10*log10(px2));
end;
data=datan;   
