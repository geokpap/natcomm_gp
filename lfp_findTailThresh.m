function [thresh, hF] = lfp_findTailThresh(filenum, varargin)

% Algorithm: start at the end of the vhisto, search downwards for a local
% max in the smoothed histo, then search the raw histo downwards from that
% value to find a string of zeros that occupies a range of values amounting
% to at least 25% of the value previously found for the local max.  The
% middle of that range will be our blink v threshold.  If no such range
% exists, then we search for the minimum in the smoothed histo that lies
% between the previously found local max and 50% of its value, and use the
% value where the minimum lies.  If there are multiple bins with the same
% minimal value, the highest bin is used.  Returns NaN if no such threshold
% exists. If lfp_vhisto figure display is enabled, a red threshold marker
% and a black curve representing the raw (unsmoothed) histogram are added
% to the plot.
%OPTIONS
% Accepts all lfp_vhisto options.  The normal defaults for 'nbins' and
% 'smooth' are overridden with 'nbins', 400, 'smooth', .05 (which work with
% both unsmoothed v and 60-point smoothed v).  Also accepts the following
% without passing it on to lfp_vhisto:
% 'gapfrac', gapfrac - <gapfrac> specifies the fraction of the value
%   found for the local max that must be occupied by zeros when identifying
%   the gap where thresh should be placed; default 0.25.  Also sets the
%   range searched (2*gapfrac) in the smoothed histo when a gap is not
%   found. 
% 'gapbins', gapbins - same as 'gapfrac' except that the number of zero
%   bins required is specified explicitly rather than as a fraction of the
%   local max position.  Specify 'gapbins', 2 to guarantee that a min will
%   be found if one exists.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

gapfrac = 0.25;
gapbins = 0;
nbins = 400;
nofigflag = false;
smooth = .05;
argnum = 1;
args2delete = [];
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'gapbins'
            argnum = argnum + 1;
            gapbins = varargin{argnum};
            args2delete(end+1:end+2) = [argnum-1 argnum];
        case 'gapfrac'
            argnum = argnum + 1;
            gapfrac = varargin{argnum};
            args2delete(end+1:end+2) = [argnum-1 argnum];
        case 'nbins'
            argnum = argnum + 1;
            nbins = varargin{argnum};
        case 'nofig'
            nofigflag = true;
        case 'smooth'
            argnum = argnum + 1;
            smooth = varargin{argnum};
    end
    argnum = argnum + 1;
end
varargin(args2delete) = [];

thresh = NaN;
[counts, smoothedcounts, binedges, values, hF] = ...
    lfp_vhisto(filenum, 'nbins', nbins, 'smooth', smooth, varargin{:});
localmaxidx = find([ false
    smoothedcounts(2:end-1) > smoothedcounts(1:end-2) ...
    & smoothedcounts(2:end-1) > smoothedcounts(3:end)
    false ]);
if isempty(localmaxidx)
    error('The smoothed histo is monotonic.');
end
blinkv = binedges(localmaxidx(end));
vperpt = (binedges(end) - binedges(1)) / length(binedges);
if gapbins
    numpts = gapbins;
else
    numpts = ceil(gapfrac*blinkv/vperpt);
end
startofzeros = 0;
for k = (localmaxidx(end)-numpts) : -1 : 1
    if ~any(counts(k : k+numpts-1))
        startofzeros = k;
        break
    end
end
if startofzeros
    thresh = binedges(round(startofzeros + numpts/2));
else
    minstart = ceil(localmaxidx(end)*(1-2*gapfrac));
    minend = localmaxidx(end);
    [C,I] = min(smoothedcounts(minend:-1:minstart));
    localminidx = minend - I + 1;
    if localminidx ~= minstart && localminidx ~= minend
        thresh = binedges(localminidx);
    end
end
if ~nofigflag
    hA = get(hF, 'CurrentAxes');
    hold(hA, 'on');
    plot(hA, binedges, counts, 'k');
    plot(hA, [thresh thresh], get(hA, 'YLim'), 'r');
    if startofzeros
        set(hA, 'YLim', ...
            [0 max(max(counts(startofzeros:startofzeros+2*numpts)), ...
            2*smoothedcounts(localmaxidx(end)) )]);
    end
end

