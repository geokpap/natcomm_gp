function [mag, units, phase] = lfp_bandpasscoh(filenums, moving_win, ...
    freqlim, nwk)
%result = lfp_bandpasscoh(filenums, moving_win, freqlim)

%result = lfp_bandpasscoh(filenums, moving_win, freqlim)
%   Function intended for use with lfp_createWave. Applies lfp_bandpass2
%   with filter frequencies <freqlim> (if <freqlim> contains only one freq,
%   then it's a lowpass filter) to the two waves in <filenums>, computes
%   multitaper coherence spectra for the moving window <moving_win>
%   spectra, and returns the magnitude averaged over <freqlim>.  <nwk> is a
%   two element row vector containing the parameters that are set by the
%   'nw' and 'k' options of lfp_spec('mt',...).  Each window has DC removed
%   and is padded to next power of 2.
%OPTIONS
%   'complex' - instead of averaging the coherence magnitude over
%       <freqlim>, computes the complex average of the coherency over
%       <freqlim> and returns its magnitude as <mag> and its phase as
%       <phase>.  Note that this may radically alter the <mag> values
%       returned if the phase is not approximately constant over <freqlim>.
%   'pad', padfactor - works as for lfp_spec('coh', ...)

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~(isequal(size(filenums), [1 2]) && isa(filenums, 'numeric'))
    error('lfp_bandpasscoh:badfilenums', ...
        '<filenums> must be a 1 x 2 numeric array.');
end
if ~(isequal(size(moving_win), [1 2]) && isa(moving_win, 'numeric'))
    error('lfp_bandpasscoh:badmoving_win', ...
        '<moving_win> must be a 1 x 2 numeric array.');
end
if ~(isequal(size(nwk), [1 2]) && isa(nwk, 'numeric'))
    error('lfp_bandpasscoh:badnwk', ...
        '<nwk> must be a 1 x 2 numeric array.');
end

units = 'unitless';
complexflag = false;
padfactor = 0;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'complex'
            complexflag = true;
        case 'pad'
            padflag = true;
            argnum = argnum + 1;
            if ~strcmp(class(varargin{argnum}), 'double')
                error('lfp_bandpasscoh:badpad', ...
                    'The padding factor must be a number');
            end
            padfactor = varargin{argnum};
        otherwise
            error('lfp_bandpasscoh:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

wave1 = lfp_bandpass2(filenums(1), freqlim);
wave2 = lfp_bandpass2(filenums(2), freqlim);

nfft = 2^(nextpow2(winpts) + pad);
[f,findx] = getfgrid(Fs,nfft,fpass); 
selfreqs = find(f >= freqlim(1) & f <= freqlim(2));
if isempty(selfreqs)
    error('lfp_bandpasscoh:nofreqs', ...
        'No frequencies selected, check parameters and lfp_FreqLim.' );
end

approx_num_windows = length(lfp_TimeStamps) * lfp_SamplesPerFrame ...
    * lfp_SamplePeriod / moving_win(2);
minperpt = 2*4e-5;
if approx_num_windows > 0.1/minperpt
    warning('lfp_bandpasscoh:manywindows', ...
        'You are computing approximately %.1d points\n(expected to run ~%.1f min @ 2.4 GHz)', ...
        approx_num_windows, approx_num_windows*minperpt);
end

[S,wsl,f]=lfp_mtspecgram2( reshape(lfp_Samples{filenums},[],1),...
    moving_win,...
    nwk,0,1/lfp_SamplePeriod,...
    lfp_FreqLim,0,0,1,0 );
[C,phi,t,f,confC,phierr,Cerr] = dg_cohgramc(data1, data2, movingwin, ...
    tapers, pad, Fs, fpass, err, trialave, rmdc);

tapers = dpsschk(tapers,N); % check tapers
for n = 1:numwin
    J1 = mtfftc(data1,tapers,nfft);
    J2 = mtfftc(data2,tapers,nfft);
    J1 = J1(findx,:,:); 
    J2 = J2(findx,:,:);
    S12 = squeeze(mean(conj(J1).*J2,2));
    S1 = squeeze(mean(conj(J1).*J1,2));
    S2 = squeeze(mean(conj(J2).*J2,2));
    C12 = S12./sqrt(S1.*S2);
    C = abs(C12);
    phi = angle(C12);
end

selfreqs = find(f >= freqlim(1) & f <= freqlim(2));
disp('The frequencies to be averaged are:');
disp(reshape(f(selfreqs),[],1));
wavedata = sum(S(:, selfreqs), 2);
t = zeros(size(wsl));
% Find timestamps of start of each window:
hWaitBar = waitbar(0, '', 'Name', 'Finding timestamps of windows');
for idx = 1:length(wsl)
    waitbar(idx/length(wsl), hWaitBar );
    t(idx) = lfp_index2time(wsl(idx));
end
close(hWaitBar);

result = lfp_interpolateWave(t, wavedata);

