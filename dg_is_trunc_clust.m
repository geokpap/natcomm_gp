function [clipping, fracclust] = dg_is_trunc_clust(pk_height)
% Applies one of two heuristics to determine whether the distribution of
% values in <pk_height> is clipped on the low-valued side, and if so, by
% how much.
%INPUTS
% pk_height: column vector with one row for each observation, containing
%   the spike's peak height.
%OUTPUTS
% clipping: 0 = not clipped; 1 = clipped symmetrical distribution; 2 =
%   clipped according to linear extrapolation test.
% fracclust: fraction of cluster present (see NOTES).
%NOTES
%   Modus operandi: estimate the fraction of the cluster that's there based
% on the assumption that a whole cluster is symmetrically distributed and
% only the bottom will get clipped.  Note that it is perfectly possible for
% <fracclust> to be greater than 1, indicating that the distribution is
% bottom-heavy rather than top-heavy.
%   Sets smoothing width to produce a smoothed histogram that does not have
% too many 3-point maxima (see "while sum(islocmax)" in code), and uses
% the tallest smoothed peak as an initial estimate of the mode of the
% distribution.  Then adjusts the position of the estimated mode so that
% the total of the raw counts within 3 bins of the estimated mode is as
% close to equal below the mode and above the mode as possible.  Finally,
% computes the ratio of the number of points that are less than the mode to
% the number that are not less than the mode; that ratio would be 1 for a
% perfectly symmetrical distribution, less than 1 for a clipped symmetrical
% distribution, and greater than 1 for a bottom-heavy distribution.
%   For top-heavy distributions, the value below which the distribution is
% clipped is estimated based on the assumption that the distribution is
% otherwise symmetrical, and if there are in fact no values below the
% estimated clipping value, then the distribution is marked as clipped
% (i.e. <clipping> is returned <true>). 
%   An alternative test (see 'lin_extrap_test' in this file) is used for
% bottom heavy distributions and cases where there are values below the
% estimated clipping value (i.e. cases where the symmetry assumption
% fails).  We determine whether it's truncated by making a linear
% extrapolation through the half-maximum and quarter-maximum points on the
% left side of the smoothed histogram and determining whether the minimum
% value in the data is to the right or left of the x-intercept.

%$Rev: 293 $
%$Date: 2022-06-01 18:51:22 -0400 (Wed, 01 Jun 2022) $
%$Author: dgibson $

% Make histogram of <pk_height>, excluding 0.1% upper outliers:
prctiles = prctile(pk_height, [50 99.9]);
binwidth = prctiles(2) / 50;
binedges = 0:binwidth:prctiles(2);
binctrs = binedges + binwidth/2;
counts = histc( pk_height( ...
    pk_height < prctiles(2) ), binedges );
counts(end) = [];
binctrs(end) = [];
% Titrate the smoothing to find the smallest odd-length window that
% produces an acceptably small number of 3-point local maxima.
smoothlen = 1;
smoo = counts;
difcnt = diff(smoo);
islocmax = difcnt(1:end-1) > 0 & difcnt(2:end) <=0;
while sum(islocmax) >= ceil(length(counts)/10)
    smoothlen = smoothlen + 2;
    hw = hanning(smoothlen);
    smoo = conv(counts, hw, 'same') / sum(hw);
    difcnt = diff(smoo);
    islocmax = difcnt(1:end-1) > 0 & difcnt(2:end) <=0;
end
locmaxidx = find(islocmax) + 1;
% Treat the tallest local maximum as the mode:
[smoomax, locmaxidx2] = max(smoo(locmaxidx));
mode_idx = locmaxidx(locmaxidx2);
% Now adjust it to find the best balance possible within a +/- 3 bin range:
idxrange = max(1, mode_idx - 3) : min(length(counts), mode_idx + 3);
newbal = NaN(size(idxrange));
for idx2 = 1:length(idxrange)
    newidx = idxrange(idx2);
    numbins = min([5, newidx - 1, length(counts) - newidx]);
    newbal(idx2) = sum(smoo(newidx - numbins : newidx - 1)) ...
        / sum(smoo(newidx + 1 : newidx + numbins));
end
[~, idx2] = min(abs(newbal - 1));
mode_idx = idxrange(idx2);

% Find peak of histo, and measure the estimated fraction of the cluster
% that's present in the raw data based on the assumption that the
% distribution is symmetrical:
numbottom = sum(pk_height < binctrs(mode_idx));
numtop = numel(pk_height) - numbottom;
fracbottom = numbottom/numtop;
fracclust = (fracbottom + 1) / 2;
% The fraction can't be greater than 1 under our assumption, but we will
% tolerate up to 1.02 because we have a finite sample:
if fracclust > 1.02
    % top-heavy distribution, linearly extrapolate:
    [clipping, fracclust] = lin_extrap_test(pk_height, smoo, ....
        smoomax, binctrs);
    return
end

% If the distribution is symmetrical other than being clipped, then there
% will be an abrupt total loss of counts in the tail below the mode, which
% accounts completely for the difference between <numtop> and <numbottom>.
% To verify our assumptions, we need to find the threshold value
% <loclipthresh> below which there will be no counts if the assumptions are
% correct.  We estimate where the lower clipping threshold <loclipthresh>
% should be by finding a higher clipping threshold (i.e. above the mode)
% <hiclipthresh> such that the number of points between the mode and
% <hiclipthresh> is the same as the number of points in the lower tail, and
% reflect it around the mode to find <loclipthresh>.  This can only be done
% if <fracbottom> is strictly less than one.  We allow for one binwidth of
% error in the estimate because it's all based on analysis of the binned
% histogram.
clipping = fracclust < 0.98;
if fracbottom < 1
    hiclipthresh = prctile( pk_height(pk_height > binctrs(mode_idx)), ...
        100 * fracbottom );
    loclipthresh = 2 * binctrs(mode_idx) - hiclipthresh - binwidth;
    if sum(pk_height < loclipthresh) > 0
        % "clipped but otherwise symmetrical" model is violated, linearly
        % extrapolate:
        [clipping, fracclust] = lin_extrap_test(pk_height, smoo, ...
            smoomax, binctrs);
    else
        % clipped but otherwise symmetrical
        clipping = 1;
    end
end

end


function [clipping, fracclust] = lin_extrap_test(pk_height, smoo, ...
    smoomax, binctrs)
% Linearly extrapolate the half-maximum and quarter-maximum points on the
% left side of the smoothed histogram, and if the minimum value in
% <pk_height> is greater than the x-intercept, then it's clipped.
halfpt = find(smoo > smoomax/2, 1);
quartpt = find(smoo > smoomax/4, 1);
slope = (smoo(halfpt) - smoo(quartpt)) ...
    / (binctrs(halfpt) - binctrs(quartpt));
intercept = binctrs(quartpt) - smoo(quartpt)/slope;
if min(pk_height) > intercept
    clipping = 2;
    % We must at least be missing the triangle from
    % <min(pk_height)> down to <intercept>, so this is an upper
    % bound on <fracclust>:
    num_missing = (min(pk_height) - intercept)^2 * slope / 2;
    fracclust = (length(pk_height) - num_missing) ...
        / length(pk_height);
else
    clipping = 0;
    fracclust = 1;
end

end

