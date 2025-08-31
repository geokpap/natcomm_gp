function [S, f, winshiftlist] = lfp_stspecgram(hWinfunc, samples, ...
    moving_win, samplefreq, labelling, padlength, rmdcflag)
%[S, f, winshiftlist] = lfp_stspecgram(hWinfunc, samples, ...
%     moving_win, samplefreq, labelling, padlength)

%  Constructs a spectrogram from the signal in the column vector <samples>
%  using a sliding window <moving_win(1)> seconds wide that moves
%  by increments <moving_win(2)> seconds.  The window function itself is
%  created by calling the function handle <hWinfunc>; the function must
%  return a column vector n elements long when called as hWinfunc(n).
%  <samplefreq>, specified in Hz, sets the time scale.  <title> is passed
%  as the title to plotspecgram.  <labelling> is a structure containing
%  fields 'title' and 'start'; <labelling.start> is the nominal time value
%  for the first datapoint (this affects the labelling of the x axis).  If
%  <labelling> is [], then no figure is displayed.  Returns the exact matrix
%  of values that are normally displayed in the spectrogram regardless of
%  the value of <labelling>, and the exact set of frequencies, and the
%  number of points by which the start of each window is shifted relative
%  to the start of <samples>.  Note that to display correctly, the matrix S
%  must be transposed before passing to imagesc and then invoke 'axis xy'
%  after calling imagesc; i.e. freq increases with column number, and time
%  increases with row number.
%
%  <winshiftlist> is the number of samples by which each window is shifted
%  relative to <samples> (so 0 indicates no shift).
%
%  If <samples> has more than one column, then each column is FFT'd
%  separately, and then the resulting power spectra are averaged to compute
%  the average spectrum.
%
%  <padlength> is optional, and defaults to 0 if not given.  It specifies
%  the number of samples by which to pad each window with zeroes before
%  computing the FFT.
%
%  Note that the window movement is done in terms of sample number in
%  integer arithmetic; consequently, if <moving_win(2)> is small, and there
%  is a large amount of time between the left end of the x-axis and the
%  reference event time, then there can be substantial truncation error by
%  the time the window gets to the reference event, and so the window that
%  includes the reference event may not be well-centered on it at all.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 6
    padlength = 0;
end

% set up
displayflag = true;
if isempty(labelling)
    displayflag = false;
    labelling.start = 0;
end
winpoints = round(moving_win * samplefreq);
wf = feval(hWinfunc, winpoints(1));
t=[ 1/samplefreq+moving_win(1)/2 : moving_win(2) : ...
        size(samples,1)/samplefreq-moving_win(1)/2 ];
t=t+labelling.start;
paddedwinpoints = winpoints(1) + padlength;
f = samplefreq * (0:floor(paddedwinpoints/2)) / paddedwinpoints;

% compute
hWaitBar = waitbar(0, '', 'Name', 'Computing spectrogram');
t_winshiftlist = 0 : moving_win(2) : ...
    lfp_SamplePeriod * (size(samples,1) - 1) - moving_win(1);
% calculate winshiftlist in samples:
winshiftlist = round(t_winshiftlist * samplefreq);
S = zeros(length(winshiftlist), length(f));
windownum = 1;
for winshift = winshiftlist
    waitbar(windownum/length(winshiftlist), hWaitBar );
    samplerange = 1 + winshift : winpoints(1) + winshift;
    winwave = samples(samplerange,:);
    if any(any(isnan(winwave)))
        warning('lfp_stspecgram:gotnan', ...
            'There is at least one NaN in the input data' );
    end
    if rmdcflag
        winwave = winwave - repmat(mean(winwave,1), size(winwave, 1), 1);
    end
    winwave = [ winwave .* repmat(wf,1,size(samples,2))
        zeros(padlength, size(samples, 2)) ];
    spectrum = fft(winwave);
    spectrum = spectrum(1:length(f),:);
    power = spectrum .* conj(spectrum) / paddedwinpoints;
    if size(samples,2) > 1
        power = mean(power, 2);
    end
    S(windownum, :) = power';
    windownum = windownum + 1;
end
close(hWaitBar);

% plot
fstart = 1;
fend = length(f);
if ~isempty(lfp_FreqLim)
    fstart = max(dg_binsearch(f, lfp_FreqLim(1)) - 1, 1);
    fend = min(dg_binsearch(f, lfp_FreqLim(2)) - 1, length(f));
    if f(fend - 1) == lfp_FreqLim(2)
        fend = fend - 1;
    end
    f = f(fstart:fend);
end
S = S(:,fstart:fend);
if displayflag
    plotspecgram(S,t,f,[],labelling.title);
end
