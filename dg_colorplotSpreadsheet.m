function dg_colorplotSpreadsheet(filename, varargin)
%dg_colorplotSpreadsheet creates a plot like Yasuo population z-score
%dg_colorplotSpreadsheet(filename)
%dg_colorplotSpreadsheet(..., 'bintuple', n)
%dg_colorplotSpreadsheet(..., 'bipolarcolor')
%dg_colorplotSpreadsheet(..., 'clim', clim)
%dg_colorplotSpreadsheet(..., 'cmap', cmap)
%dg_colorplotSpreadsheet(..., 'labels', labeltable)
%dg_colorplotSpreadsheet(..., 'midbins', n)
%dg_colorplotSpreadsheet(..., 'nrows', n)
%dg_colorplotSpreadsheet(..., 'rowheight', r)
%INPUTS
% filename - pathname to a tab-delimited text or *.xls binary Excel
%   spreadsheet file.  If there is a header row, then it must contain a
%   non-numerical string (i.e. a string that cannot be converted to a
%   number by Matlab's 'str2num' function) for every column of data. Column
%   1 contains stage number, which is used to group rows as follows.  A
%   series of consecutive rows having the same stage number are assumed to
%   represent different event-centered time windows for that stage.  Column
%   2 contains the event ID on which the time window is centered (a warning
%   is issued if this changes between stages).  Column 3 contains the
%   number of units for that stage (a warning is issued if this changes
%   within the stage).  Column 4 is empty.  The remaining columns contain
%   time-binned data (e.g. z-score averaged over all units) or empty cells
%   that print as a white space.
%OPTIONS
% 'bintuple', n - instead of displaying each bin, displays the average of
%   each group of <n> bins.  Bins are grouped subject to the constraint
%   that there be an even number of groups (and therefore that the center
%   of the series of grouped bins is at a group boundary).  If grouping
%   results in bins being left over, then the same number of bins are
%   dropped from the beginning and the end of the series, unless there is
%   an odd number to be dropped, in which case one more is dropped from the
%   end.
% 'bipolarcolor', nyellow - uses dg_bipolarcolorscale as colormap, with
%   the same number of colors as the specified (or default) colormap.
%   <nyellow> becomes the <numyellows> parameter to dg_bipolarcolorscale
%   if it is non-negative; otherwise the default <numyellows> is used.
% 'labels', labeltable - overrides the default value for the table that
%   converts from event IDs to event names.  <labeltable> is simply a cell
%   string array whose kth element contains the name for event ID k.
% 'clim', clim - sets limits on color scale to range from clim(1) to
%   clim(2).  (Note that this is difficult to do after the figure has
%   already been drawn because of the special colorbar that doesn't show
%   the NaN color and the fact that if clim(1) is set to a value higher
%   than the minimum of the data, clipped values will be shown in the same
%   color as NaN.)
% 'cmap', cmap - uses the colormap <cmap> instead of the default.
% 'midbins', n - uses only the <n> bins from the middle range of bins on
%   each row; if this results in an odd number of bins being dropped, then
%   one more is dropped from the end than from the beginning.
% 'nrows', n - add extra empty rows to bring the total number of rows to at
%   least abs(n); if n>0 then the empties are added after the end, and if 
%   n<0 then they are added before the beginning.
% 'rowheight', r - set rows to fixed height, expressed as a fraction of
%   window height (i.e. a number between 0 and 1; 0.2 approximates the
%   behavior of Yasuo's Delphi program).
% 'ttext', option that triggers the labeling of each line with its session
%   number (as they are not all continuous) and changes side text from
%   "Stage" to "sess".

%$Rev: 153 $
%$Date: 2012-07-17 18:40:53 -0400 (Tue, 17 Jul 2012) $
%$Author: dgibson $

% Default values:
bintuple = 1;
bipolarcolor = false;
clim = [];
cmap = [];
labeltable{1} = 'Start Rec';
labeltable{2} = 'End Rec';
labeltable{3} = 'BL On';
labeltable{4} = 'BL Off';
labeltable{10} = 'Click';
labeltable{11} = 'Gate';
labeltable{12} = 'Locomotion';
labeltable{13} = 'Start';
labeltable{14} = 'Turn';
labeltable{15} = 'RT Off';
labeltable{16} = 'LT Off';
labeltable{17} = 'R Goal';
labeltable{18} = 'L Goal';
labeltable{21} = 'Rough';
labeltable{22} = 'Smooth';
labeltable{30} = 'Noise';
labeltable{31} = '1 kHz';
labeltable{38} = '8 kHz';
labeltable{40} = 'Noise Off';
labeltable{41} = '1 kHz Off';
labeltable{48} = '8 kHz Off';
labeltable{39} = 'Tone On'; % substitute for [31 38]
labeltable{49} = 'Tone Off'; % substitute for [41 48]
labeltable{54} = 'Turn'; % substitute for 14
labeltable{56} = 'Turn Off'; % substitute for [15 16]
labeltable{78} = 'Goal'; % substitute for [17 18]
midbins = 0;
nrows = 0;
rowheight = 0;
ttext = false;

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'bintuple'
            argnum = argnum + 1;
            bintuple = varargin{argnum};
        case 'bipolarcolor'
            argnum = argnum + 1;
            nyellow = varargin{argnum};
            bipolarcolor = true;
        case 'clim'
            argnum = argnum + 1;
            clim = varargin{argnum};
        case 'cmap'
            argnum = argnum + 1;
            cmap = varargin{argnum};
        case 'labels'
            argnum = argnum + 1;
            labeltable = varargin{argnum};
        case 'midbins'
            argnum = argnum + 1;
            midbins = varargin{argnum};
        case 'nrows'
            argnum = argnum + 1;
            nrows = varargin{argnum};
        case 'rowheight'
            argnum = argnum + 1;
            rowheight = varargin{argnum};
        case 'ttext'
            ttext = true;
        otherwise
            error('dg_colorplotSpreadsheet:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

[pathstr, name, ext] = fileparts(filename);
switch lower(ext)
    case '.xls'
        spreadsheet = xlsread(filename);
    otherwise
        spreadsheet = dg_tabread(filename);
end

% <spreadsheet> now contains an exact image of the spreadsheet minus the
% first row, with NaN in place of empty cells.

% <startnewstage> is a list of row numbers that start a new stage
startnewstage = ...
    [ 1; find(spreadsheet(2:end, 1) ~= spreadsheet(1:end-1, 1)) + 1 ];
rowsPerStage = startnewstage(2:end) - startnewstage(1:end-1);
rowsPerStage(end+1) = size(spreadsheet,1) + 1 - startnewstage(end);
if any(rowsPerStage ~= rowsPerStage(1))
    badstageidx = find(rowsPerStage ~= rowsPerStage(1));
    error('dg_colorplotSpreadsheet:badstage', ...
        'Stage %d contains a different number of rows', ...
        spreadsheet(startnewstage(badstageidx(1)), 1) );
end
rowsPerStage = rowsPerStage(1);

% Calculate bin parameters
nbinsgross = size(spreadsheet,2) - 4;
if midbins == 0
    midbins = nbinsgross;
end
ngroups = floor(midbins/bintuple);
if bintuple > 1 && mod(ngroups, 2) == 1
    ngroups = ngroups - 1;
end
nbins2use = ngroups * bintuple;
firstbin = floor((nbinsgross - nbins2use)/2) + 1;
lastbin = firstbin + nbins2use - 1;

% Convert spreadsheet bin contents to display format:
events = spreadsheet(1:rowsPerStage, 2);
displaydata = [];  % matrix for display
unitcountstr = '';  % string to display unit counts
for firstrow = startnewstage(1:end)'
    if ~isequal(...
            spreadsheet(firstrow : (firstrow + rowsPerStage - 1), 2), ...
            events )
        warning('dg_colorplotSpreadsheet:badevents', ...
            'The events list mismatches for stage %d.', ...
            spreadsheet(firstrow, 1) );
    end
    for k = 1:bintuple
        bins2avg(:,:,k) = spreadsheet( ...
            firstrow : (firstrow + rowsPerStage - 1), ...
            firstbin+4 + k - 1 : bintuple : lastbin+4 ); 
    end
    displaydata = [ displaydata
        reshape(mean(bins2avg,3)', [], 1)' ];
    if any(spreadsheet(firstrow : (firstrow + rowsPerStage - 1), 3) ...
            ~= spreadsheet(firstrow, 3) )
        warning('dg_colorplotSpreadsheet:badunitcount', ...
            'The unit counts for stage %d are not consistent.', ...
            spreadsheet(firstrow, 1) );
    end
    if ttext
        unitcountstr = sprintf('%ssess=%2d units=%3d\n', ...
            unitcountstr, spreadsheet(firstrow, 1), spreadsheet(firstrow, 3) );
    else
        unitcountstr = sprintf('%sStage=%2d #Units=%3d\n', ...
            unitcountstr, spreadsheet(firstrow, 1), spreadsheet(firstrow, 3) );
    end
end
if nrows > 0 && size(displaydata,1) < nrows
    displaydata(end+1:nrows, :) = NaN;
elseif nrows < 0 && size(displaydata,1) < -nrows
    nblanks = -nrows - size(displaydata,1);
    displaydata = [ zeros(nblanks, size(displaydata,2))
        displaydata ];
    displaydata(1:nblanks, :) = NaN;
end

% Set up x-axis ticks & labels
numticks = 2 * length(events) + 1;
ticklabel(numticks) = {''};
tickspacing = ((0:(2*length(events)))/2) * ngroups + 0.5;
if max(events) > length(labeltable)
    labeltable{max(events)} = '';
end
for k = 1:length(labeltable)
    if isempty(labeltable{k})
        labeltable{k} = sprintf('evt %d', k');
    end
end
ticklabel(2:2:numticks) = labeltable(events);

% Figure must exist before getting colormap:
hF = figure;
oldcmap = get(hF, 'ColorMap');
if isempty(cmap)
    cmap = oldcmap;
end

% Calculate color scale so that cmap(1) is not occupied by real data:
if isempty(clim)
    clim = [ min(min(displaydata)) max(max(displaydata)) ];
else
    if min(min(displaydata)) < clim(1)
        displaydata(displaydata < clim(1)) = clim(1);
    end
end
if bipolarcolor
    if nyellow >= 0
        cmap = dg_bipolarcolorscale(clim, size(cmap,1), nyellow);
    else
        cmap = dg_bipolarcolorscale(clim, size(cmap,1));
    end
end
cmap(1,:) = [1 1 1];    % used for plotting NaN
nlevels = size(cmap,1) - 1;
levelwidth = (clim(2) - clim(1)) / nlevels;
clim(1) = clim(1) - (1+10*eps) * levelwidth;

set(hF, 'DefaulttextInterpreter', 'none', 'ColorMap', cmap);
hI = imagesc(displaydata, clim);
if ttext
    hT = text(1.3, 1.1, unitcountstr, ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontName', 'Ariel', 'FontSize', 8);
else
    hT = text(1.3, 1.1, unitcountstr, ...
        'Units', 'normalized', 'VerticalAlignment', 'top');
end
title(sprintf(...
    '%s\ngross bins: %d; bins used:%d; bins per display bin: %d\n\n', ...
    filename, nbinsgross, nbins2use, bintuple ));
hA = get(hI,'Parent');
if rowheight
    axesheight = size(displaydata,1) * rowheight;
    bottom = 0.9 - axesheight;
else
    axesheight = 0.8;
    bottom = 0.1;
end
set(hA, ...
    'Position',[.1 bottom .6 axesheight], ...
    'CLim', clim );
hCB = colorbar('peer', hA);
hCBI = findobj(hCB,'Type','image');
set(hCBI, 'CData', (2:size(cmap,1))');
set(hCB,'Position', [.8 .6 .02 .3], 'YAxisLocation', 'left');
set(hA, 'XTick', tickspacing);
set(hA, 'XTickLabel', ticklabel);
set(hA, 'Tickdir', 'out');

% set the yticks according to session number
if ttext
    set(hA, 'YTick', 1:length(startnewstage))
    set(hA, 'YTickLabel', spreadsheet(startnewstage,1))
end

% Annoyingly, to restore the old color map in subsequent figures without
% changing the color map in figure(hF), it seems we must create a new
% figure:
hF2 = figure;
colormap(oldcmap);
close(hF2);

