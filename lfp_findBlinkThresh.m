function [thresh, hF] = lfp_findBlinkThresh(filenum, varargin)

% Algorithm: start at the end of the vhisto, search downwards for a local
% max in the smoothed histo, then search the raw histo downwards from that
% value to find a string of zeros that occupies a range of values amounting
% to at least 25% of the value previously found for the local max.  The
% middle of that range will be our blink v threshold.  Returns NaN if no
% such threshold exists.
%OPTIONS
% Accepts all lfp_vhisto options.  However, specifying 'nbins' or 'smooth'
% will result in a confusing figure title and possibly also poor
% performance (the values specified in the call to lfp_vhisto work with
% both unsmoothed v and 60-point smoothed v).

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

nofigflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'nofig'
            nofigflag = true;
    end
    argnum = argnum + 1;
end

thresh = NaN;
[counts, smoothedcounts, binedges, values, hF] = ...
    lfp_vhisto(filenum, 'nbins', 400, 'smooth', .05, varargin{:});
localmaxidx = find([ false
    smoothedcounts(2:end-1) > smoothedcounts(1:end-2) ...
    & smoothedcounts(2:end-1) > smoothedcounts(3:end)
    false ]);
blinkv = binedges(localmaxidx(end));
vperpt = (binedges(end) - binedges(1)) / length(binedges);
numpts = ceil(blinkv/(4*vperpt));
startofzeros = 0;
for k = (localmaxidx(end)-numpts) : -1 : 1
    if ~any(counts(k : k+numpts-1))
        startofzeros = k;
        break
    end
end
if startofzeros
    thresh = binedges(round(startofzeros + numpts/2));
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
        set(hA, 'XLim', [0 500]);
    end
end

