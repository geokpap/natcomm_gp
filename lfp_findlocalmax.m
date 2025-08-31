function samples = lfp_findlocalmax(filenum, starts, stops)
%LFP_FINDLOCALMAX finds a local maximum value of the waveform in
% filenum in each of the intervals between corresponding elements of
% <starts> and <stops>. Returns a row vector of the sample numbers of the
% maximum values.

%samples = lfp_findlocalmax(filenum, starts, stops)
% <starts> and <stops> are arrays of sample indices containing the same
% numbers of elements.  The local maximum is found by starting at the
% center of the interval and doing a linear search in the direction of
% the greatest neighboring value.  Any values in <starts> or <stops> that
% exceed the bounds of lfp_Samples{filenum} are replaced by the
% corresponding endpoint, i.e. values less than 1 become 1, values greater
% than numel(lfp_Samples{filenum}) become numel(lfp_Samples{filenum}).  In
% any case where there are 3 or fewer points between starts(k) and stops(k)
% inclusive, the index of the maximum value is returned.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_Samples;

if ~(numel(starts) == numel(stops))
    error('lfp_findlocalmax:mismatch', ...
        '<starts> and <stops> have different numbers of elements.');
end

% Truncate all start and stop values to bounds of lfp_Samples{filenum}:
badindices = find(starts < 1);
if ~isempty(badindices)
    starts(badindices) = 1;
end
badindices = find(stops < 1);
if ~isempty(badindices)
    stops(badindices) = 1;
end
badindices = find(starts > numel(lfp_Samples{filenum}));
if ~isempty(badindices)
    starts(badindices) = numel(lfp_Samples{filenum});
end
badindices = find(stops > numel(lfp_Samples{filenum}));
if ~isempty(badindices)
    stops(badindices) = numel(lfp_Samples{filenum});
end

samples = zeros(1, numel(starts));
for idx = 1 : numel(starts)
    if (stops(idx) - starts(idx)) < 3
        [y, samples(idx)] = max(lfp_Samples{filenum}(starts(idx):stops(idx)));
    else
        center = round((starts(idx) + stops(idx))/2);
        if lfp_Samples{filenum}(center+1) > lfp_Samples{filenum}(center-1)
            incr = 1;
            boundary = stops(idx);
        else
            incr = -1;
            boundary = starts(idx);
        end
        for snum = center + incr : incr : boundary - incr
            if lfp_Samples{filenum}(snum + incr) < lfp_Samples{filenum}(snum)
                break
            end
        end
        if lfp_Samples{filenum}(snum + incr) >= lfp_Samples{filenum}(snum)
            snum = boundary;
        end
        samples(idx) = snum;
    end
    if samples(idx) == starts(idx) ...
            || samples(idx) == stops(idx)
        warning('lfp_findlocalmax:endpoint', ...
            'The max value was at %d in %d:%d', ...
            samples(idx), starts(idx), stops(idx) );
    end
end