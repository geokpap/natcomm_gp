function ts = lfp_minDurNeg(filenum, thresh, mindur, varargin)
%LFP_MINDURNEG finds points where a waveform remains below a threshold
% for a minimum amount of time.
%ts = lfp_minDurNeg(filenum, thresh, mindur)
% For every negative-going threshold crossing (<thresh>) in <filenum>,
% determines whether the next positive-going crossing of the same threshold
% is at least <mindur> seconds later than the negative-going crossing.
% Returns a list of timestamps for such points.  The sample timestamps
% returned are the ones immediately preceding the threshold crossings.
%OPTIONS:
% 'excludeNaN' - causes points containing NaN to be treated as if they were
%   positive-going threshold crossings.
% 'bounds'- <ts> is a two-element cell vector where ts{1} contains
%   the same values that are returned by default, and ts{2} contains
%   timestamps of the positive thresh crossings that go with the negative
%   thresh crossings in ts{1}, marking the ends of the minimum duration
%   below-threshold periods.

%$Rev: 123 $
%$Date: 2010-05-06 19:56:19 -0400 (Thu, 06 May 2010) $
%$Author: dgibson $

lfp_declareGlobals;

argnum = 1;
excludeNaNflag = false;
boundsflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'excludeNaN'
            excludeNaNflag = true;
        case 'bounds'
            boundsflag = true;
        otherwise
            error('lfp_minDurNeg:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~( strcmp(class(filenum), 'double') ...
        && strcmp(class(thresh), 'double') ...
        && strcmp(class(mindur), 'double') )
    error('lfp_minDurNeg:badarg', ...
        'All arguments must be numeric');
end

if ~(isequal(size(filenum), [1 1]) && (fix(filenum) == filenum))
    error('lfp_minDurNeg:badfilenum', ...
        '<filenum> must be a single integer');
end

if ~ismember(filenum, lfp_ActiveFilenums)
    error('lfp_minDurNeg:badfilenum2', ...
        'filenum %d does not exist', filenum);
end

[pxings nxings] = lfp_findThresholds(thresh, thresh, filenum, ...
    [1 numel(lfp_Samples{filenum})] );
if excludeNaNflag
    nanstarts = dg_findNaNruns(lfp_Samples{filenum});
end
goodxing = repmat(false, size(nxings));
if boundsflag
    goodendTS = NaN(size(nxings));
    goodendptr = 1;
end
for nxindex = 1:length(nxings)
    pxindex = find(pxings > nxings(nxindex));
    if excludeNaNflag
        nanindex = find(nanstarts > nxings(nxindex));
    else
        nanindex = [];
    end
    if isempty(pxindex)
        if isempty(nanindex)
            endsample = numel(lfp_Samples{filenum});
        else
            endsample = nanstarts(nanindex(1));
        end
    else
        if isempty(nanindex)
            endsample = pxings(pxindex(1));
        else
            endsample = min(nanstarts(nanindex(1)), pxings(pxindex(1)));
        end
    end
    goodxing(nxindex) = lfp_index2time(endsample) ...
        > (lfp_index2time(nxings(nxindex)) + mindur);
    if boundsflag
        goodendTS(nxindex) = lfp_index2time(endsample);
        goodendptr = goodendptr + 1;
    end
end
goodnxindex = find(goodxing);
if boundsflag
    ts{1} = lfp_index2time(nxings(goodnxindex));
    ts{2} = goodendTS(goodnxindex);
else
    ts = lfp_index2time(nxings(goodnxindex));
end