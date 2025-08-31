function ts = lfp_minDur(filenum, thresh, mindur)
%LFP_MINDUR finds points where a waveform exceeds a threshold for a minimum
% amount of time.
%ts = lfp_minDur(filenum, thresh, mindur)
% For every positive-going threshold crossing (<thresh>) in <filenum>,
% determines whether the next negative-going crossing of the same threshold
% is at least <mindur> seconds later than the positive-going crossing.
% Returns a list of timestamps for such points.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~( strcmp(class(filenum), 'double') ...
        && strcmp(class(thresh), 'double') ...
        && strcmp(class(mindur), 'double') )
    error('lfp_minDur:badarg', ...
        'All arguments must be numeric');
end

if ~(isequal(size(filenum), [1 1]) && (fix(filenum) == filenum))
    error('lfp_minDur:badfilenum', ...
        '<filenum> must be a single integer');
end

if ~ismember(filenum, lfp_ActiveFilenums)
    error('lfp_minDur:badfilenum2', ...
        'filenum %d does not exist', filenum);
end

[pxings nxings] = lfp_findThresholds(thresh, thresh, filenum, ...
    [1 numel(lfp_Samples{filenum})] );
goodxing = repmat(false, size(pxings));
for pxindex = 1:length(pxings)
    nxindex = find(nxings > pxings(pxindex));
    if isempty(nxindex)
        goodxing(pxindex) = lfp_index2time(numel(lfp_Samples{filenum})) ...
            > (lfp_index2time(pxings(pxindex)) + mindur);
    else
        goodxing(pxindex) = lfp_index2time(nxings(nxindex(1))) ...
            > (lfp_index2time(pxings(pxindex)) + mindur);
    end
end
goodpxindex = find(goodxing);
ts = lfp_index2time(pxings(goodpxindex));
