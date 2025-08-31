function samples = lfp_findmax(filenum, starts, stops)
%LFP_FINDMAX finds the maximum value of the waveform in filenum in each of
% the intervals between corresponding elements of <starts> and <stops>.
% Returns a row vector of the sample numbers of the maximum values.
%NOTE: This differs from lfp_findlocalmax in that it simply finds the
%absolute maximum in each interval.  DG 29-Oct-2008 

%samples = lfp_findmax(filenum, starts, stops)
% <starts> and <stops> are arrays of sample indices containing the same
% numbers of elements.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_Samples;

if ~(numel(starts) == numel(stops))
    error('lfp_findmax:mismatch', ...
        '<starts> and <stops> have different numbers of elements.');
end

samples = zeros(1, numel(starts));
for idx = 1 : numel(starts)
    [c samples(idx)] = max(lfp_Samples{filenum}(starts(idx):stops(idx)));
    samples(idx) = samples(idx) + starts(idx) - 1;
    if samples(idx) == starts(idx) ...
            || samples(idx) == stops(idx)
        warning('lfp_findmax:endpoint', ...
            'The max value was at %d in %d:%d', ...
            samples(idx), starts(idx), stops(idx) );
    end
end