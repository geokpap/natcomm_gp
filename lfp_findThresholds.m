function [pos_xings, neg_xings] = lfp_findThresholds(posthresh, negthresh, ...
    filenum, samplerange)
%LFP_FINDTHRESHOLDS finds threshold crossings in a waveform
%[pos_xings neg_xings] = lfp_findThresholds(posthresh, negthresh, ...
%     filenum, samplerange)

% Finds the positive-going crossings of <posthreshold> and the
% negative-going crossings of <negthreshold> in the waveform in <filenum>,
% limited to the range of samples starting at index number <samplerange>(1)
% and ending at index number <samplerange>(2).  If <samplerange> is missing
% or [], then the entire waveform is processed.  For each crossing, returns
% the sample index of the sample before the first sample that is beyond the
% threshold, with positive-going crossings in pos_xings and negative-going
% crossings in neg_xings.  pos_xings and neg_xings are both column vectors.
%
% Note that since a crossing requires that the second of the two points be
% strictly greater than threshold, whereas the first can be equal to
% threshold, it is possible to get e.g. two positive crossings in a row
% when <posthreshold> = <negthreshold> if the wave just touches threshold
% and then goes up again.
%
% pos_xings is a list of indices into
% lfp_Samples(samplerange(1):samplerange(2)), so the indices into
% lfp_Samples of those same samples are pos_xings + samplerange(1) - 1.
% The indices into lfp_Samples of the first samples beyond threshold is
% pos_xings + samplerange(1).  Similarly for neg_xings.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if (nargin < 4) || (isempty(samplerange))
    samplerange(2) = numel(lfp_Samples{filenum});
    samplerange(1) = 1;
end

pos_xings = find( ...
    (lfp_Samples{filenum}(samplerange(1) + 1 : samplerange(2)) > posthresh) ...
    & (lfp_Samples{filenum}(samplerange(1) : samplerange(2) - 1) <= posthresh) );
neg_xings = find( ...
    (lfp_Samples{filenum}(samplerange(1) + 1 : samplerange(2)) < negthresh) ...
    & (lfp_Samples{filenum}(samplerange(1) : samplerange(2) - 1) >= negthresh) );
