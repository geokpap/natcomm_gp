function ts = lfp_signifRateChange(trial, clustnum, window, varargin)
% Computes the average number of spikes per bin during a baseline period,
% and searches for "a significant increase in firing rate" as follows.  A
% threshold is set either at a specified number of standard deviations
% from the baseline mean or at a specified p-level calculated using Poisson
% statistics.  Starting at each bin, the average firing rate over a series
% of N consecutive bins is computed, and the first bin in the series is
% marked as to whether the series starting at that bin was above threshold
% or not.  When M consecutive bins are found to be above threshold, the
% start time of the first of the M consecutive bins is reported as the time
% of the "significant increase in firing rate".  <ts> is NaN if no such bin
% is found.

%INPUTS
% trial, clustnum, window - as usual; <window> denotes the time period
%   relative to the alignment event which is to be searched for a
%   significant increase.
%OUTPUT
% ts - the timestamp in seconds of the earliest significant increase in
%   firing rate
%OPTIONS
% 'baseline', basewindow - <basewindow> is a 2-element array like <window>
%   specifying the baseline period relative to the laignment event.
%   Default [-300 0].
% 'binwidth', binwidth - <binwidth> is in milliseconds; default 10.
% 'm', M - <M> is the number of consecutive moving average firing rates
%   that must be above threshold.  Default 5.
% 'n', N - <N> is the number of bins used for calculating the moving
%   average firing rate.  Default 5.
% 'p', p - <p> is the p-level for the Poissonian threshold calculation. 
%   Default is 0, which invokes the mean + K*SD formula.
% 'plothis' - creates a figure window showing the histogram of counts
%   in each bin, and the threshold.
% 'SDs', K - <K> is the value used for threshold = mean + K*SD, where <SD>
%   is the standard deviation.  Note that this has no effect if <p> is
%   nonzero.  Default 3.

%$Rev: 120 $
%$Date: 2010-04-23 20:59:08 -0400 (Fri, 23 Apr 2010) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 3
    window = [];
end

if isempty(window) && ~isempty(lfp_XLimAll)
    window = lfp_XLimAll;
end
if ~isempty(window) && ...
        ~(length(window)==2 && numel(window)==2 && isnumeric(window))
    error('lfp_signifRateChange:window', ...
        '<window> must be a two-element numeric array');
end

ts = NaN;
basewindow = [-300 0];
M = 5;
N = 5;
p = 0;
plothisflag = false;
K = 3;
binwidth = 10;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'baseline'
            argnum = argnum + 1;
            basewindow = varargin{argnum};
        case 'binwidth'
            argnum = argnum + 1;
            binwidth = varargin{argnum};
        case 'm'
            argnum = argnum + 1;
            M = varargin{argnum};
        case 'n'
            argnum = argnum + 1;
            N = varargin{argnum};
        case 'p'
            argnum = argnum + 1;
            p = varargin{argnum};
        case 'plothis'
            plothisflag = true;
        case 'SDs'
            argnum = argnum + 1;
            K = varargin{argnum};
        otherwise
            error('lfp_signifRateChange:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end
bw = binwidth * 1e-3;

trialevents = ...
    lfp_Events(lfp_TrialIndex(trial,1):lfp_TrialIndex(trial,2), :);
alignevtidx = find(ismember(trialevents(:,2), lfp_AlignmentRef));
if isempty(alignevtidx)
    error('lfp_signifRateChange:noalign', ...
        'Could not find alignment event %s in trial %d', ...
        dg_thing2str(lfp_AlignmentRef), trial);
end
alignevtidx = alignevtidx(1) + lfp_TrialIndex(trial,1) - 1;
alignTS = lfp_Events(alignevtidx,1);
if isempty(window)
    window = [lfp_Events(lfp_TrialIndex(trial,1),1) ...
        lfp_Events(lfp_TrialIndex(trial,2),1)] - alignTS;
end

binsafter = floor(window(2)/bw);
binsbefore = floor(-window(1)/bw);
binedges = (-binsbefore:binsafter) * bw;
basebinsafter = floor(basewindow(2)/bw);
basebinsbefore = floor(-basewindow(1)/bw);
basebinedges = (-basebinsbefore:basebinsafter) * bw;

binedgeTS = binedges + alignTS;
basebinedgeTS = basebinedges + alignTS;
starttime = min(binedgeTS(1), basebinedgeTS(1));
endtime = max(binedgeTS(end), basebinedgeTS(end));
spikes2count = lfp_Spikes{clustnum}( ...
    lfp_Spikes{clustnum} >= starttime & lfp_Spikes{clustnum} <= endtime );

basecounts = histc(spikes2count, basebinedgeTS);
basecounts(end) = [];
counts = histc(spikes2count, binedgeTS);
counts(end) = [];
overthresh = false(1, length(counts) - N + 1);
if p == 0
    meanbaserate = sum(basecounts)/(basebinedgeTS(end) - basebinedgeTS(1));
    basestd = std(basecounts/bw);
    thresh_rate = meanbaserate + K * basestd;
    threshline = thresh_rate * bw;
    for k = 1:length(overthresh)
        overthresh(k) = sum(counts(k:(k+N-1)))/(N*bw) > thresh_rate;
        if k >= M && all(overthresh((k-M+1):k))
            break;
        end
    end
else
    meanbasecount = mean(basecounts) * N;   % expected count in N bins
    maxspikes = 2000 * bw; % I assume nothing can fire above 2 kHz
    mycdf = poisscdf(1:maxspikes, meanbasecount);
    thresh_count = find(mycdf > (1 - p));
    if isempty(thresh_count)
        error('lfp_signifRateChange:thresh_count', ...
            'You will never see the rate change you are looking for.');
    end
    thresh_count = thresh_count(1);
    threshline = thresh_count/N;
    for k = 1:length(overthresh)
        overthresh(k) = sum(counts(k:(k+N-1))) > thresh_count;
        if k >= M && all(overthresh((k-M+1):k))
            break;
        end
    end
end
if k >= M && all(overthresh((k-M+1):k))
    ts = binedgeTS(k-M+1);
end
if plothisflag
    hF = figure;
    hA = axes('Parent', hF);
    bar(hA, binedges(1:end-1) + median(diff(binedges))/2, counts);
    set(hA, 'NextPlot', 'add');
    plot(get(hA, 'XLim'), [threshline threshline], 'r');
    title(hA, sprintf('lfp_signifRateChange(%d, %d, [%.3g %.3g], %s)', ...
        trial, clustnum, window(1), window(2), dg_thing2str(varargin)), ...
        'Interpreter', 'none');
end


