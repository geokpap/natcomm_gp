function [thresh, cnts, binctrs] = dg_S2M_findvalley(x)
% Computes 500-point histogram of values in <x>, smooths it, and finds the
% last valley.  Returns empty if valley search fails.
%INPUTS
% x
%OUTPUTS
% thresh: upper limit of the x-value range covered by the bin at the bottom
%   of the valley.
% cnts: vector of counts in 500 bins.
% binctrs: x value at center of each of the 500 bins.

%$Rev: 259 $
%$Date: 2017-08-02 14:37:57 -0400 (Wed, 02 Aug 2017) $
%$Author: dgibson $

thresh = [];
[cnts, binctrs] = hist(x, 500);
binedges = binctrs + median(diff(binctrs))/2; % upper edges
hw = hanning(101);
smcnts = conv(cnts, hw, 'same') / sum(hw);
diffsmcnts = diff(smcnts);
% If there is one and only one run of zeros in <smcnts>, then we
% use it as our valley:
zeroidx = find(smcnts == 0);
isbadthresh = true;
if ~isempty(zeroidx) && all(diff(zeroidx) == 1)
    threshbin = round(mean(zeroidx([1 end])));
    thresh = binedges(threshbin);
    if sum(x > thresh) >= 3
        isbadthresh = false;
    end
end
if isbadthresh
    % resort to desperate measures: keep finding the last valley
    % before the current <threshbin> that could be used as a threshold
    % until there are at least three <x> values and at least one bin above
    % it.
    threshbin = length(diffsmcnts);
    while sum(x > binedges(threshbin)) < 3 ...
            || length(diffsmcnts) == threshbin
        pkidx = find(diffsmcnts(1:threshbin) > 0, 1, 'last');
        if isempty(pkidx)
            return;
        end
        % search backwards for valley start and valley end (which should be
        % the same point if the value of smcnts is non-zero at the bottom
        % of the valley).
        srchidx = pkidx-1;
        valend = srchidx;
        valstart = 0;
        while (valstart == 0)
            if srchidx < 1
                break
            end
            if diffsmcnts(srchidx) > 0
                valend = srchidx;
                srchidx = srchidx - 1;
            elseif diffsmcnts(srchidx) == 0
                srchidx = srchidx - 1;
            else
                valstart = srchidx;
            end
        end
        if valstart == 0
            % There is no valley
            break
        end
        threshbin = round((valstart + valend) / 2);
    end
    thresh = binedges(threshbin);
end
