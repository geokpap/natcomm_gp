function [bands, extrastarts, extraends] = dg_findfreqbands(corr, f, varargin)
%[bands, extrastarts, extraends] = findfreqbands(corr, f)
% Plots a fig showing band starts with green lines and band ends with red
% lines.  The contour line used to define bands is shown in black, the band
% starts are shown in green, and the band ends in red, and the threshold is
% shown in grey: a band is a region where the black line is further from
% the diagonal than the grey line. The threshold is defined as 2 standard
% deviations greater than the median distance within <freqlim>, and the 1%
% outliers (smallest 0.5% and largest 0.5%) are removed from the distance
% data before calculating SD.
%INPUTS
% corr: cross-covariance matrix
% f: frequencies for elements of <corr>
%OUTPUTS
% bands: 2 column list of lower and upper limits of auto-covarying
%   frequency bands.
% extrastarts, extraends: parts of bands that were incompletely identified.
%OPTIONS
% 'freqlim', freqlim: the band within which to count <bands> and within
%   which to compute the median distance of the contour line from the
%   diagonal in the <corr> matrix.
% 'maxbands', maxbands: maximum number of <bands> to tolerate within
%   <freqlim>.
% 'minbw', bw - removes any frequency bands that are less then <bw> Hz wide
%   before counting the number of bands in <freqlim>.  Default <bw> = 0.
% 'numSDs', N - set the number of SDs from the median used for finding
%   bands.  Default = 2.
% 'plotfig': what it says.
% 'xcovthresh', thresh: the xcov value at which to plot the contour line
%   that defines the bands, relative to the maximum value in the matrix.
% 'titlestr', titlestr: extra stuff to tack onto the end of the figure
%   title.

%$Rev: 227 $
%$Date: 2015-09-22 17:34:59 -0400 (Tue, 22 Sep 2015) $
%$Author: dgibson $

bw = 0;
freqlim = [5 100];
maxbands = 5;
numSDs = 2;
plotfig = false;
xcovthresh = 0.7;
titlestr = '';

if ~isempty(varargin)
    argnum = 1;
    while argnum <= length(varargin)
        switch varargin{argnum}
            case 'freqlim'
                argnum = argnum + 1;
                freqlim = varargin{argnum};
            case 'numSDs'
                argnum = argnum + 1;
                numSDs = varargin{argnum};
            case 'maxbands'
                argnum = argnum + 1;
                maxbands = varargin{argnum};
            case 'minbw'
                argnum = argnum + 1;
                bw = varargin{argnum};
            case 'plotfig'
                plotfig = true;
            case 'xcovthresh'
                argnum = argnum + 1;
                xcovthresh = varargin{argnum};
            case 'titlestr'
                argnum = argnum + 1;
                titlestr = varargin{argnum};
        end
        argnum = argnum + 1;
    end
end

[bands, extrastarts, extraends, distthresh, C] = ...
    findfreqbands2(corr, f, xcovthresh, freqlim, numSDs);
if bw > 0
    istoosmall = bands(:,2) - bands(:,1) < bw;
    bands(istoosmall, :) = [];
end
band_in_freqlim = all(bands >= freqlim(1) & bands <= freqlim(2), 2);
bands = bands(band_in_freqlim, :);
if plotfig
    f_in_freqlim = f >= freqlim(1) & f <= freqlim(2);
    figure;
    imagesc(f(f_in_freqlim), f(f_in_freqlim), ...
        corr(f_in_freqlim, f_in_freqlim));
    axis xy;
    colorbar;
    caxis([0 1]);
    hold on;
    for k = 1:size(bands,1);
        plot([0 2*bands(k,1)], [2*bands(k,1) 0], ...
            'Color', [0 .7 0]);
        plot([0 2*bands(k,2)], [2*bands(k,2) 0], 'r');
    end
    dg_plotContours(C, gca, 1, 'k');
    if ~isempty(distthresh)
        plot(f, f + distthresh * sqrt(2), 'Color', [1 1 1]*0.5);
    end
    if isempty(titlestr)
        titlestr = sprintf( ...
            'xcovthresh=%.2g, freqlim=%s, maxbands=%d', ...
            xcovthresh, dg_thing2str(freqlim), maxbands);
    else
        titlestr = sprintf( ...
            'xcovthresh=%.2g, freqlim=%s, maxbands=%d\n%s', ...
            xcovthresh, dg_thing2str(freqlim), maxbands, titlestr);
    end
    title(titlestr);
end
end

function [bands, extrastarts, extraends, distthresh, C] = ...
    findfreqbands2(corr, f, xcovthresh, freqlim, numSDs)
bands = [];
extrastarts = [];
extraends = [];
distthresh = [];
normthresh = max(corr(:)) * xcovthresh;
C = contourc(f, f, corr, normthresh * [1 1]);
% Find each contour line in <C> together with its length
contouridx = [];
contourlen = [];
idx = 1;
while idx <= size(C,2)
    contouridx(end+1,1) = idx; %#ok<AGROW>
    contourlen(end+1, 1) = C(2,idx); %#ok<AGROW>
    idx = idx + C(2,idx) + 1;
end
[longestlen, longestidx] = max(contourlen);
x = C(1, contouridx(longestidx) + 1 : contouridx(longestidx) + longestlen);
y = C(2, contouridx(longestidx) + 1 : contouridx(longestidx) + longestlen);
% <dist2diag>, <freq> are of same length as <x>, <y>, i.e. the number of
% points in the contour.  This can be much larger than <f>.
dist2diag = reshape((x - y) / sqrt(2), [], 1); % column vector
freq = reshape((x + y) / 2, [], 1);
lastpt = longestlen;
freqdiff = diff(freq);
% Matlab 'ctourc' contours can start anywhere, and include any portion of a
% closed loop up to a whole one.  That is a pain in the ass to deal with,
% to wit:
if any(freqdiff < 0)
    % Check to see if freqs are a symmetric (but not necessarily closed)
    % loop that should be cut in half. There are two possibilities: the max
    % and min freqs are represented only once, which means that one of them
    % is the exact center of symmetry, or they are represented twice, which
    % means that one of them should be a consecutive pair of points which
    % together are the center of symmetry.
    [c,minidx] = min(freq); %#ok<ASGLU>
    [c,maxidx] = max(freq); %#ok<ASGLU>
    if minidx(1) > maxidx(1)
        % The main sequence is descending
        editidx = minidx(1):-1:maxidx(1);
    else
        editidx = minidx(1):maxidx(1);
    end
    dist2diag = dist2diag(editidx);
    freq = freq(editidx);
    lastpt = length(freq);
    freqdiff = diff(freq);
    if sum(freqdiff < 0) > sum (freqdiff > 0)
        % freqs are reversed
        freq = freq(end:-1:1);
        freqdiff = diff(freq);
        dist2diag = dist2diag(end:-1:1);
    end
    if any(freqdiff < 0)
        warning('findfreqbands:nonmono', ...
            'Contour is non-monotonic with f, truncating to monotonic part');
        lastpt = find(freqdiff < 0) - 1;
        lastpt = lastpt(1);
    end
end
if sum(dist2diag(1:lastpt) < 0) > sum(dist2diag(1:lastpt) > 0)
    % Got the upper side of the diagonal instead of lower; invert.
    dist2diag = -dist2diag;
end
is_in_freqlim = freq >= freqlim(1) & freq <= freqlim(2);
if ~any(is_in_freqlim)
    return
end
dist2diag = dist2diag(is_in_freqlim);
freq = freq(is_in_freqlim);
% out-of-band freq points have now been trimmed away from <freq>
% and <dist2diag>
prct = prctile(dist2diag, [0.5 99.5]);
SD = std(dist2diag(dist2diag > prct(1) & dist2diag < prct(2)));
distthresh = nanmedian(dist2diag) + numSDs * SD;
if isnan(distthresh)
    error('findfreqbands:ohdear', ...
        'Oh dear, oh dear, oh dear...');
end
reldist = dist2diag - distthresh;
bandendidx = reshape(find( ...
    reldist(2:end) <= 0 & reldist(1:end-1) > 0), [], 1);
% <bandendidx> is an index into <reldist>, <dist2diag>, <freq>,
% <x>, or <y>.
bandends = NaN(size(bandendidx));
for k = 1:length(bandends)
    rd = reldist(bandendidx(k) + [0 1]);
    F = freq(bandendidx(k) + [0 1]);
    [rd, sortidx] = sort(rd);
    bandends(k) = interp1(rd, F(sortidx), 0);
end
bandstartidx = reshape(find( ...
    reldist(1:end-1) <= 0 & reldist(2:end) > 0), [], 1);
bandstarts = NaN(size(bandstartidx));
for k = 1:length(bandstarts)
    rd = reldist(bandstartidx(k) + [0 1]);
    F = freq(bandstartidx(k) + [0 1]);
    [rd, sortidx] = sort(rd);
    bandstarts(k) = interp1(rd, F(sortidx), 0);
end
if ~isempty(bandstarts) && ~isempty(bandends)
    [pairs, extraL, extraR] = dg_zip(bandstarts, bandends);
    bands = [bandstarts(pairs(:,1)) bandends(pairs(:,2))];
    extrastarts = bandstarts(extraL);
    extraends = bandends(extraR);
else
    bands = [];
    extrastarts = bandstarts;
    extraends = bandends;
end
end
