function ts = lfp_findSpikes(filenum, thresh)
%ts = lfp_findSpikes(filenum, thresh)
%   Searches CSC channel <filenum> for waveforms that resemble action
% potentials.  It is assumed that the channel is suitably high-pass
% filtered to remove high amplitude low frequency components, but not so
% much as to distort a typical spike waveform (e.g. with a cutoff frequency
% in the 50 - 100 Hz range).
%INPUTS
% filenum: filenum of the CSC channel to analyze.
% thresh: the spike waveform must exceed this threshold in the
%   positive-going direction.
%OUTPUTS
% ts: vector of timestamps denoting the peak time of each spike.
%NOTES
% For each sample at which the signal crosses <thresh>, the interval
% <win> around that sample is analyzed to determine whether it contains a
% plausible spike waveform by applying the following criteria:
%   0. The entire window around the trigger must exist.
%   1. It must not be a local minimum that just touches <thresh>; this
%   makes threshold crossing a three-point criterion, so it is defined to
%   be unsatisfied at the first sample.
%   2. It must not contain any samples of magnitude greater than
%   10*<thresh>.
%   3. The maximum must be in the middle 50% of <win>.
%   4. The last sample must be within <thresh>/2 of the first sample.
% A good default value for <thresh> is 4*SD, where lfp_XLimAll and
% lfp_AlignmentRef have been set appropriately to select the time window to
% be analyzed, and then lfp_findbadtrials with an appropriate 'magfactor'
% is used to eliminate trials with huge signals in the analysis window:
%     lfp_BadTrials = lfp_findbadtrials(filenum, ...
%         'windows', {lfp_AlignmentRef}, lfp_XLimAll, ...
%         'magfactor', 10);
%     sampledata = lfp_getSamples([], filenum, []);
%     SD = std(sampledata(:));
% Some highpass filtering is necessary to achieve reliable triggering on
% potential spikes, but the cutoff frequency should be no higher than 100
% Hz to avoid grossly distorting the spike waveform (50 Hz is better if
% that doesn't disturb the triggering too much).  Some lowpass filtering is
% also desirable so that <peakidx> doesn't get thrown off by high frequency
% noise.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_Samples lfp_SamplePeriod

win = [-0.0005 0.0015]; % window around trigger in seconds
winidx = round(win(1)/lfp_SamplePeriod) : round(win(end)/lfp_SamplePeriod);
q1idx = round(length(winidx)/4);
q3idx = round(3*length(winidx)/4);

% <xingidx> points at the last sample that was below <thresh>:
xingidx = find( [ false
    (lfp_Samples{filenum}(3 : end) > thresh) ...
    & (lfp_Samples{filenum}(2 : end - 1) <= thresh) ...
    & (lfp_Samples{filenum}(1 : end - 2) <= thresh) ] );
peakidx = zeros(size(xingidx));
peakidx2 = 1;
for xingidx2 = 1:length(xingidx)
    if xingidx(xingidx2) + winidx(1) < 1 || ...
            xingidx(xingidx2) + winidx(end) > numel(lfp_Samples{filenum})
        continue
    end
    if any( abs(lfp_Samples{filenum}(xingidx(xingidx2) + winidx)) ...
            > 10 * thresh ) ...
            || lfp_Samples{filenum}(xingidx(xingidx2) + winidx(end)) ...
            - lfp_Samples{filenum}(xingidx(xingidx2) + winidx(1)) ...
            > thresh/2
        continue
    else
        [~, maxidx] = max(lfp_Samples{filenum}(xingidx(xingidx2) + winidx));
        if maxidx < q1idx || maxidx > q3idx
            continue
        end
        peakidx(peakidx2) = maxidx + xingidx(xingidx2) + winidx(1) - 1;
        peakidx2 = peakidx2 + 1;
    end
end
peakidx(find(peakidx==0, 1) : end) = [];
ts = lfp_index2time(peakidx);
        

