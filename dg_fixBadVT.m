function badtimes = dg_fixBadVT(src, dest, badtimesdest, plotflag)
%badtimes = dg_fixBadVT(src)
%dg_fixBadVT(src, dest)
%dg_fixBadVT(src, dest, badtimesdest)
%dg_fixBadVT(src, dest, badtimesdest, plotflag)
% This function is meant specifically to fix Neuralynx video tracker files
% where the animal's location is frequently given as (0,0).  All segments
% of both channels (x and y) where the value is 0 are replaced by a simple
% linear interpolation from the surrounding nonzero frames.  Leading and
% trailing runs of zeros are left untouched. If <dest> is not given, then
% the destination is the same as <src> but with the string '_repaired'
% appended to the filename before the extension.  No checking is done to
% see if <dest> already exists. <badtimesdest> is the pathname to a .mat
% file containing <badtimes>.  If not given, then it is the same as <dest>
% but with '_badtimes' appended to the filename before the extension.  
%INPUTS
% src: path to Neuralynx video tracker file
% dest: where to put the repaired video tracker file; specifying ''
%   suppresses output
% badtimesdest: where to put the badtimes file; specifying ''
%   suppresses output
% plotflag: if true, plots histograms of the lengths of runs of consecutive
%   zeros together with the threshold for including a run in <badtimes>.
%OUTPUT
% badtimes: a 2-column array containing the Neuralynx (microsecond)
%   timestamps of the starting frame and ending frame of each run of zeros
%   that is in the longest 5% of runs.

%$Rev: 74 $
%$Date: 2010-08-20 17:17:42 -0400 (Fri, 20 Aug 2010) $
%$Author: dgibson $

if nargin < 2
    [p n e] = fileparts(src);
    dest = fullfile(p, [n '_repaired' e]);
end
if nargin < 3
    [p n e] = fileparts(dest);
    badtimesdest = fullfile(p, [n '_badtimes.mat']);
end
if nargin < 4
    plotflag = false;
end
   
[ts x y phi t pts h] = Nlx2MatVT_411(src, [1 1 1 1 1 1], 1, 1);
[x, zerostartidx, zeroendidx, zerolengths, hAx] = ...
    interpzeropts(x, 'X', plotflag);
[y, s, e, l, hAy] = interpzeropts(y, 'Y', plotflag);
discordant = (x==0) & (y~=0) | (x~=0) & (y==0);
if any(discordant)
    fprintf(1, 'There were %d frames where only one of X and Y was zero', ...
        sum(discordant));
end
thresh = prctile(zerolengths, 95);
islongrun = zerolengths > thresh;
badtimes = [reshape(ts(zerostartidx(islongrun)), [], 1) ...
    reshape(ts(zeroendidx(islongrun)), [], 1)];
if plotflag
    set(hAx, 'NextPlot', 'add');
    plot(hAx, [thresh thresh], get(hAx, 'YLim'), 'r');
    set(hAy, 'NextPlot', 'add');
    plot(hAy, [thresh thresh], get(hAy, 'YLim'), 'r');
    dg_plotBadVT(ts, x, badtimes);
    title('X');
    dg_plotBadVT(ts, y, badtimes);
    title('Y');
end
if ~isempty(dest)
    Mat2NlxVT_411(dest, 0, 1, 1, length(ts), [1 1 1 1 1 1 1], ...
        ts, x, y, phi, t, pts, h);
end
if ~isempty(badtimesdest)
    save(badtimesdest, 'badtimes', '-mat');
end
end

function [v, zerostartidx, zeroendidx, zerolengths, hA] = ...
    interpzeropts(v, label, plotflag)
iszero = reshape(v == 0, 1, []);
numzeros = sum(iszero);
fprintf(1, '%2.0f%% (%d/%d) of %s positions were zero\n', ...
    100*numzeros/length(v), numzeros, length(v), label);
zerostartidx = find([false iszero(2:end) & ~iszero(1:end-1)]);
zeroendidx = find([iszero(1:end-1) & ~iszero(2:end) false]);
if zeroendidx(1) < zerostartidx(1)
    zeroendidx(1) = [];
end
if zerostartidx(end) > zeroendidx(end)
    zerostartidx(end) = [];
end
if length(zerostartidx) ~= length(zeroendidx)
    error('dg_fixBadVT:oops', ...
        'zeroidx length mismatch');
end
zerolengths = zeroendidx - zerostartidx + 1;
for runnum = 1:length(zerostartidx)
    v(zerostartidx(runnum) - 1 : zeroendidx(runnum) + 1) = ...
        linspace(v(zerostartidx(runnum) - 1), ...
        v(zeroendidx(runnum) + 1), ...
        zeroendidx(runnum) - zerostartidx(runnum) + 3);
end
if plotflag
    hF = figure;
    hA = axes('Parent', hF);
    hist(hA, zerolengths, min(zerolengths):max(zerolengths));
    title(label);
    ylabel('Number of runs of zeros');
    xlabel('Length of run of zeros');
else
    hA = [];
end
end

