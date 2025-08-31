function result = lfp_bandpower(filenum, moving_win, freqlim, varargin)
%result = lfp_bandpower(filenum, moving_win, freqlim)

%result = lfp_bandpower(filenum, moving_win, freqlim)
%   Function intended for use with lfp_createWave.
%   Computes a series of FFTs of the wave in filenum using time intervals
%   specified by moving_win (same format as lfp_spec, but no default
%   values) and a Hamming window.  Each window has DC removed and is padded
%   to next power of 2.  The energies in all frequency components that are
%   in the range of freqlim (specified as [min max]) are summed to yield a
%   time series with time points spaced at intervals of <moving_win>(2).
%   That result is linearly interpolated between the time points to yield a
%   waveform with the same sampling rate as the CSC data. WARNING: the wave
%   in <filenum> is treated as if it were continuous in time, without
%   regard to recording stops and restarts. NOTE: Cannot disable action of
%   lfp_FreqLim due to call to lfp_stspecgram, which respects lfp_FreqLim.
%OPTIONS
%  'mt', [nw k] - computes multitaper spectra instead of Hamming window
%       spectra.  [nw k] is a two element row vector containing the
%       parameters that are set by the 'nw' and 'k' options of
%       lfp_spec('mt',...).  Each window has DC removed and is padded to
%       next power of 2.

% 5/16/07: DG changed rmdcflag from 0 to 1 in call to lfp_stspecgram.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~isequal(size(filenum), [1 1])
    error('lfp_bandpower:badfilenum', ...
        '<filenum> must be a single value.');
end

mtflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'mt'
            mtflag = true;
            argnum = argnum + 1;
            if length(varargin) < argnum
                error('lfp_bandpower:badoption2', ...
                    'You must specify [nw k] for ''mt''.' );
            end
            nwk = varargin{argnum};
        otherwise
            error('lfp_bandpower:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

approx_f = 0 : 1/moving_win(1) : 1/(2*lfp_SamplePeriod);
selfreqs = find(approx_f >= freqlim(1) & approx_f <= freqlim(2));
if ~isempty(lfp_FreqLim)
    selfreqs = find(approx_f(selfreqs) >= lfp_FreqLim(1) ...
        & approx_f(selfreqs) <= lfp_FreqLim(2));
end
if isempty(selfreqs)
    error('lfp_bandpower:nofreqs', ...
        'No frequencies selected, check parameters and lfp_FreqLim.' );
end

approx_num_windows = length(lfp_TimeStamps) * lfp_SamplesPerFrame ...
    * lfp_SamplePeriod / moving_win(2);
if mtflag
    minperpt = 4e-5;
else
    minperpt = 1e-4;
end
if approx_num_windows > 0.1/minperpt
    warning('lfp_bandpower:manywindows', ...
        'You are computing approximately %.1d points\n(expected to run ~%.1f min @ 2.4 GHz)', ...
        approx_num_windows, approx_num_windows*minperpt);
end

if mtflag
    [S,wsl,f]=lfp_mtspecgram2( reshape(lfp_Samples{filenum},[],1),...
        moving_win,...
        nwk,0,1/lfp_SamplePeriod,...
        lfp_FreqLim,0,0,1,0 );
else
    winpoints = round(moving_win / lfp_SamplePeriod);
    exp = 0;
    while 2^exp < winpoints(1)
        exp = exp + 1;
    end
    padlength = 2^exp - winpoints(1);
    [S f wsl] = lfp_stspecgram(@hamming, ...
        reshape(lfp_Samples{filenum},[],1), moving_win, ...
        1/lfp_SamplePeriod, [], padlength, 1);
end
selfreqs = find(f >= freqlim(1) & f <= freqlim(2));
disp('The frequencies to be summed are:');
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
if ~mtflag
    % Convert to middle of each window:
    t = t + moving_win(1)/2;
end

result = lfp_interpolateWave(t, wavedata);

