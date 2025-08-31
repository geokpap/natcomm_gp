function result = lfp_movxcov(filenums, moving_win, lag)
%result = lfp_movxcov(filenums, moving_win, lag)

%result = lfp_movxcov(filenums, moving_win, lag)
%   Function intended for use with lfp_createWave. Computes a series of
%   cross-covariances of the two waves in <filenums> using time intervals
%   specified by moving_win (same format as lfp_spec, but no default
%   values).  <lag> may be positive or negative, and is an offset measured
%   in seconds that is applied to the filenums(1) before retrieving the
%   sample data for each window.  It thus plays the same role as the
%   "lags" output of Matlab's xcorr and xcov functions.  Note that this
%   implies that the nominal time at which the xcov is computed is the
%   center of the window in filenums(2), which for non-zero lag is NOT the
%   same as the center of the window in filenums(1). For each filenum,
%   after a windowful of sample data is retrieved, it is multiplied by a
%   Hanning window.  That result is then divided by the square root of the
%   sum of squares of its elements to normalize, so that the final
%   cross-covariance value would be 1 for two (possibly scaled) copies of
%   the same waveform as input. The cross-covariances for the series of
%   windows are then linearly interpolated to yield a waveform with the
%   same sampling rate as the CSC data. WARNING: the wave in <filenums> is
%   treated as if it were continuous in time, without regard to recording
%   stops and restarts.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~isequal(size(filenums), [1 2])
    error('lfp_movxcov:badfilenums', ...
        '<filenums> must be a two-element row vector.');
end

approx_num_windows = length(lfp_TimeStamps) * lfp_SamplesPerFrame ...
    * lfp_SamplePeriod / moving_win(2);
if approx_num_windows > 10000
    warning('lfp_movxcov:manywindows', ...
        'You are computing approximately %d points, at about 20 sec per 10000', ...
        approx_num_windows);
end

% compute wsl (window shift list)
nlag = round(lag/lfp_SamplePeriod);
nsamples = numel(lfp_Samples{filenums(1)});
winwidth = round(moving_win(1) / lfp_SamplePeriod);
if nlag > 0
    wsl = 0 : round(moving_win(2)/lfp_SamplePeriod) : ...
        (nsamples - 1) - winwidth - nlag;
else
    wsl = -nlag : round(moving_win(2)/lfp_SamplePeriod) : ...
        (nsamples - 1) - winwidth;
end

% compute xcov at nlag for each window
hWaitBar = waitbar(0, '', 'Name', 'Computing xcovs');
wavedata = zeros(length(wsl), 1);
hw = hanning(winwidth);
for winnum = 1:length(wsl)
    samplerange = (1 : winwidth) + wsl(winnum);
    s1 = lfp_Samples{filenums(1)}(samplerange + nlag) ...
        - mean(lfp_Samples{filenums(1)}(samplerange + nlag)) ;
    s2 = lfp_Samples{filenums(2)}(samplerange) ...
        - mean(lfp_Samples{filenums(2)}(samplerange));
    s1 = hw .* s1';
    s2 = hw .* s2';
    wavedata(winnum) = sum(s1 .* s2) ...
        / (sqrt(sum(s1.^2)) * sqrt(sum(s2.^2)));
    waitbar(winnum/length(wsl), hWaitBar );
end
close(hWaitBar);

% Find timestamps of start of each window:
t = zeros(size(wsl));
hWaitBar = waitbar(0, '', 'Name', 'Finding timestamps of windows');
for idx = 1:length(wsl)
    waitbar(idx/length(wsl), hWaitBar );
    t(idx) = lfp_index2time(wsl(idx));
end
close(hWaitBar);
% Convert to middle of each window:
t = t + moving_win(1)/2;

result = lfp_interpolateWave(t, wavedata);
