function dg_homeTimeAnalysis(winlen, radius, varargin)
%OPTIONS
% 'label', labelflag - puts title on each plot and labels on x axis if
%   <labelflag> is 1, name of file's directory as ylabel if 0 (default).
% 'maxrows', maxrows - creates figures with no more than <maxrows> rows of
%   plots on each page.
% 'daysfirst' - plots subdirectories whose names begin "day" in a separate
%   figure first, then plots others.
% 'daysonly' - plots only subdirectories whose names begin "day".

%$Rev: 37 $
%$Date: 2009-07-30 14:49:55 -0400 (Thu, 30 Jul 2009) $
%$Author: dgibson $

daysfirstflag = false;
daysonlyflag = false;
labelflag = false;
maxrows = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'daysfirst'
            daysfirstflag = true;
        case 'daysonly'
            daysonlyflag = true;
        case 'label'
            argnum = argnum + 1;
            labelflag = varargin{argnum};
        case 'maxrows'
            argnum = argnum + 1;
            maxrows = varargin{argnum};
            if ~isnumeric(maxrows)
                error('dg_homeTimeAnalysis:maxrows', ...
                    'maxrows must be a number' );
            end
            maxrows = maxrows(1);
        otherwise
            error('dg_homeTimeAnalysis:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
dataroot = uigetdir('', 'Choose Data Root Directory');

% Enter an animal ID (typical filename: TT14534--040031.R01, ID = 4534)
animalID = inputdlg('Enter animal ID:');
animalID = animalID{1};

% Program finds all files that belong to that animal ID.
[files, dates, dirnames] = findActivityFiles(dataroot, animalID);
if daysfirstflag || daysonlyflag
    isday = cellfun(@strncmpi, ...
        dirnames, ...
        repmat({'day'}, size(dirnames)), ...
        repmat({3}, size(dirnames)) );
    makeFigures(files(isday), dates(isday), dirnames(isday), ...
        maxrows, winlen, radius, ...
        labelflag, animalID);
    if ~daysonlyflag
        makeFigures(files(~isday), dates(~isday), dirnames(~isday), ...
            maxrows, winlen, radius, ...
            labelflag, animalID);
    end
else
    makeFigures(files, dates, dirnames, maxrows, winlen, radius, ...
        labelflag, animalID);
end

end

function makeFigures(files, dates, dirnames, maxrows, winlen, radius, ...
    labelflag, animalID)
xmax = 16.5;
ymax = 16.5;
[dates, IX] = sort(dates);
files = files(IX);
dirnames = dirnames(IX);

% Run the 'splitT', 20 analysis on all files and leave results in open fig
% windows, with plots grouped into before-split in col. 1 and after-split
% in col. 2, with sessions chronologically by row.  The sessions go:
% saline, day1 - day7, challenge. Occupancy histos on one fig window,
% fractional home time on other fig window.  Must be flexible in case this
% turns out to be too many subplots per fig window.
hFocc = figure;
hFhometime = figure;
hFpolar = figure;
if isempty(maxrows)
    maxrows = length(files);
end
rownum = 1;
for k = 1:length(files)
    if rownum > maxrows
        hFocc = figure;
        hFhometime = figure;
        hFpolar = figure;
        rownum = 1;
    end
    hA(1,1) = subplot(maxrows, 2, rownum*2-1, 'Parent', hFocc);
    hA(2,1) = subplot(maxrows, 2, rownum*2, 'Parent', hFocc);
    hA(1,2) = subplot(maxrows, 2, rownum*2-1, 'Parent', hFhometime);
    hA(2,2) = subplot(maxrows, 2, rownum*2, 'Parent', hFhometime);
    if labelflag
        labelval = 1;
    else
        labelval = dirnames{k};
    end
    [homebase, isathome, xy, splitpoint] = dg_homeTime(files{k}, ...
        0.5, winlen, radius,'axes', hA, 'splitT', 20, ...
        'xmax', xmax, 'ymax', ymax, 'label', labelval);
    hApolar = plotpolarposition(0.5, xy, splitpoint, xmax, ymax, ...
        maxrows, rownum, hFpolar, labelval);
    if ~labelflag && rownum == 1
        hT = title(hA(1,1), sprintf('Animal ID: %s', animalID));
        set(hT,'Units', 'normalized', 'Position', [1 1.065]);
        hT = title(hA(1,2), sprintf('Animal ID: %s', animalID));
        set(hT,'Units', 'normalized', 'Position', [1 1.065]);
        hT = title(hApolar(1,1), sprintf('Animal ID: %s', animalID));
        set(hT,'Units', 'normalized', 'Position', [1 1.065]);
    end
    rownum = rownum + 1;
end

end


function [files, dates, dirnames] = findActivityFiles(datadir, animalID)
[p, dirname] = fileparts(datadir);
files = cell(0,1);
dates = zeros(0,1);
dirnames = cell(0,1);
filelist = dir(datadir);
for k=1:length(filelist)
    if filelist(k).isdir
        if ~isequal(filelist(k).name, '.') ...
                && ~isequal(filelist(k).name, '..')
            [f, dt, dn] = findActivityFiles( ...
                fullfile(datadir, filelist(k).name), animalID );
            files = [files; f];
            dates = [dates; dt];
            dirnames = [dirnames; dn];
        end
    elseif ~isempty(regexpi(filelist(k).name, ...
            [animalID '[a-zA-Z]?--\d*.R\d*$'] ))
        files{end+1} = fullfile(datadir, filelist(k).name);
        dates(end+1) = filelist(k).datenum;
        dirnames{end+1} = dirname;
    end
end
end

function hA = plotpolarposition(Ts, xy, splitpoint, xmax, ymax, ...
    numrows, rownum, hFpolar, labelval)
% <rownum> refers to the rows in the other figures; each "row" in this
% figure is actually two rows in the paired grid.  <numrows> likewise.
% <hA> contains axes handles for before splitpoint on row 1, after on row
% 2; radius in col. 1, angle in col. 2.  <labelval> as for dg_homeTime.
center = [xmax ymax]/2;
xy = xy - repmat(center, size(xy, 1), 1);
r = sqrt(sum(xy.^2, 2));
rmax = sqrt(sum(center.^2));
theta = -360 * atan2(xy(:,2), xy(:,1)) / (2*pi); % "-" for image coords
% Find discontinuities caused by wrapping around the 180-degree mark: for
% this purpose, any pair of successive values of opposite sign where at
% least one member of the pair has magnitude greater than 90 degrees.  For
% computational convenience, if either member of the pair is exactly zero,
% then it isn't treated as a discontinuity.
disconts = find( (sign(theta(1:end-1)) == -sign(theta(2:end))) & ...
    (abs(theta(1:end-1)) > 90 | abs(theta(2:end)) > 90) );
if isequal(labelval, 1)
    xmargin = 0.1;
    ymargin = 0.1;
else
    xmargin = 0.1;
    ymargin = 0.02;
end
if isequal(labelval, 0)
    topmargin = ymargin;
    bottommargin = ymargin;
elseif isequal(labelval, 1)
    topmargin = 0.1;
    bottommargin = 0.1;
else
    topmargin = 0.05;
    bottommargin = 0.05;
end
timepts = (0:length(r)-1) * Ts;
hA(1,1) = dg_subplot(2*numrows, 2, [2*rownum-1 1], 'paired', ...
    'xmargin', xmargin, 'ymargin', ymargin, ...
    'topmargin', topmargin, 'bottommargin', bottommargin);
plot(hA(1,1), timepts(1:splitpoint), r(1:splitpoint));
ylim(hA(1,1), [0 rmax]);
set(hA(1,1), 'XTickLabel', '');
hA(2,1) = dg_subplot(2*numrows, 2, [2*rownum-1 2], 'paired', ...
    'xmargin', xmargin, 'ymargin', ymargin, ...
    'topmargin', topmargin, 'bottommargin', bottommargin);
plot(hA(2,1), timepts(splitpoint+1:end), r(splitpoint+1:end));
ylim(hA(2,1), [0 rmax]);
set(hA(2,1), 'XTickLabel', '');
% Plot theta in segments separated by discontinuities
hA(1,2) = dg_subplot(2*numrows, 2, [2*rownum 1], 'paired', ...
    'xmargin', xmargin, 'ymargin', ymargin, ...
    'topmargin', topmargin, 'bottommargin', bottommargin);
hold(hA(1,2), 'on');
numseg = sum(disconts < splitpoint);
segstart = 1;
for segnum = 1:numseg
    plot(hA(1,2), timepts(segstart:disconts(segnum)), ...
        theta(segstart:disconts(segnum)));
    segstart = disconts(segnum) + 1;
end
plot(hA(1,2), timepts(segstart:splitpoint), theta(segstart:splitpoint));
ylim(hA(1,2), [-180 180]);
hA(2,2) = dg_subplot(2*numrows, 2, [2*rownum 2], 'paired', ...
    'xmargin', xmargin, 'ymargin', ymargin, ...
    'topmargin', topmargin, 'bottommargin', bottommargin);
hold(hA(2,2), 'on');
numseg = sum(disconts > splitpoint);
if numseg
    firstseg = find(disconts > splitpoint);
    firstseg = firstseg(1);
else
    firstseg = 0;
end
segstart = splitpoint + 1;
for segnum = firstseg:(firstseg + numseg - 1)
    plot(hA(2,2), timepts(segstart:disconts(segnum)), ...
        theta(segstart:disconts(segnum)));
    segstart = disconts(segnum) + 1;
end
plot(hA(2,2), timepts(segstart:end), theta(segstart:end));
ylim(hA(2,2), [-180 180]);
switch labelval
    case 1
        title(hA(1,1), 'Position in Polar Coordinates');
        ylabel(hA(1,1), 'r');
        ylabel(hA(1,2), 'theta');
        xlabel(hA(1,2), 'Elapsed Time, s');
        ylabel(hA(2,1), 'r');
        ylabel(hA(2,2), 'theta');
        xlabel(hA(2,2), 'Elapsed Time, s');
    case 0
        set(hA(1,2), 'XTickLabel', '');
        set(hA(2,2), 'XTickLabel', '');
    otherwise
        ylabel(hA(1,2), labelval);
        ylabel(hA(2,2), labelval);
        if rownum ~= numrows
            set(hA(1,2), 'XTickLabel', '');
            set(hA(2,2), 'XTickLabel', '');
        end
end
end


