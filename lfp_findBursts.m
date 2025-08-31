function timestamps = lfp_findBursts(filenum, testband, refband, ...
    testthresh, refthresh, clipvals, moving_win)
%timestamps = lfp_findBursts(filenum, testband, refband, thresh, clipvals)
% Function intended for use with lfp_createEvents.  Finds band-limited
% bursts in wave data.  <timestamps> is a cell vector whose first element
% contains onsets and second element contains offsets.
% 1.	Replace clipped values with NaN (note this will not work on "weird"
%       clipping), defined as <= clipvals(1) or >= clipvals(2)
% 2.	compute power & stats in test band and reference band
% 3.	create discriminator signals at thresh*std for both bands
% 4.	mark onset and offset of (test>thresh*std) && (ref<thresh*std); by
%       definition, onset and offset refer to the comparison between the
%       current window and the previous, so the first window cannot be an
%       onset or offset, whereas the last can be either.
%DEFAULTS
% testthresh, refthresh = 7
% clipvals = [-2048 2047]
% moving_win = [1 0.25];

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 7
    moving_win = [1 0.25];
end
if nargin < 6
    clipvals = [-2048 2047];
end
if nargin < 4
    testthresh = 7;
    refthresh = 7;
end

if ~isequal(size(filenum), [1 1])
    error('lfp_findBursts:badfilenums', ...
        '<filenum> must have exactly 1 element' );
end

% 1.	Replace clipped values with NaN
samples = reshape(lfp_Samples{filenum}, [], 1);
samples(samples <= clipvals(1) | samples >= clipvals(2)) = NaN;

% 2.	compute power & stats in test band and reference band
[bandpower, wints] = lfp_bandpower2(samples, moving_win, ...
    [testband; refband], 0);
if all(isnan(bandpower(:,1)))
    error('lfp_findBursts:baddata1', ...
        'There were no time windows free of clipped values' )
end
goodwindows = ~isnan(bandpower(:,1));
testmean = mean(bandpower(goodwindows,1));
teststd = std(bandpower(goodwindows,1));
refmean = mean(bandpower(goodwindows,2));
refstd = std(bandpower(goodwindows,2));

% 3.	create discriminator signals at thresh*std for both bands
testdisc = bandpower(:,1) > testthresh * teststd + testmean;
if ~any(testdisc)
    error('lfp_findBursts:badthresh1', ...
        'No time window had suprathreshold power in test band' );
end
refdisc = bandpower(:,2) < refthresh * refstd + refmean;
if ~any(testdisc)
    error('lfp_findBursts:badthresh2', ...
        'No time window had subthreshold power in reference band' );
end
gotburst = testdisc & refdisc;
if ~any(gotburst)
    error('lfp_findBursts:badthresh3', ...
        'No time window had a burst' );
end

% 4.	mark onset and offset of gotburst
onsets = [false; gotburst(2:end) & ~gotburst(1:end-1)];
offsets = [false; ~gotburst(2:end) & gotburst(1:end-1)];
timestamps = { wints(onsets) wints(offsets) };
