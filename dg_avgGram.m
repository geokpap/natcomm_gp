function [hI, hCB, stuff] = dg_avgGram(figlist, varargin)
%[hI, hCB] = dg_avgGram(figlist)
%   Averages the data in a list of currently open or previously saved
%   figure files and displays the result.  Works on output of
%   lfp_makePasteup as well as dg_showGram.  The figure files must contain
%   at least two images, the first of which contains the colorbar and the
%   remainder of which contain the data.  This is the normal result of
%   creating the figure by calling imagesc first and colorbar second.  The
%   axes labellings, 'XData', and 'YData' are copied from the last figure
%   in the list to the newly created figure, along with any line objects
%   that have been drawn on the same axes as the gram.  All files must
%   contain the same number of images.
%
%   Use dg_recolorGram to set color scale.

% INPUTS:
%   figlist - a list of figures to average; if of class 'double', then each
%   element should be the number of an open figure; if of class 'cell'
%   then each element should be a filename string, either absolute or
%   relative to current working directory.  May also be of class
%   'matlab.ui.Figure', which is the figure handle as returned in R2014b
%   and later by the functions 'figure' and 'gcf'. If figures are pasteups
%   e.g. from lfp_makepasteup, then they must all contain the same number
%   of panels for the same alignment events.
% OUTPUTS:
%   [hI, hCB]:  handles to new image(s) and colorbar(s) (see 'std' option).
%       If the input figures are multi-panel pasteups, then hI points to
%       the first panel(s) of the output figure(s).
%   stuff:  a structure containing miscellaneous values returned by
%   	various options.
% OPTIONS:
% 'alignL' - truncates different width grams on the right as needed,
%   keeping left edges in register.
% 'baseline' - treats the first figure as a baseline period to be
%   subtracted from the average of the remaining figures; the average of
%   the first figure is taken in the X direction, and is replicated as many
%   times as needed to create a baseline gram that fits the remaining
%   figures in the X direction.  The number of elements in the Y direction
%   must still be the same as in the following figures.  The resulting
%   time-independent baseline gram is then subtracted from the result of
%   averaging the remaining grams.  If 'log' is specified, then the logs of
%   the baseline gram and the averaged gram are taken first, then the
%   subtraction is done.  You cannot specify 'baseline' and 'norm' at the
%   same time, but if there is a legitimate reason to combine them you can
%   apply one, and then call dg_avgGram again to apply the other (you pick
%   the order).  The baseline fig must NOT be a pasteup, but the remaining
%   figs must all still have the same number of images.  Returns baseline
%   spectrum in stuff.basespectrum.
% 'diff' - treats the first figure as a reference spectrogram to be
%   subtracted from the average of the remaining figures.  This is similar
%   to 'baseline', but without the averaging in the X direction, so each
%   individual time window (i.e. each column of pixels) in the first figure
%   gets subtracted from the corresponding window in the average of the
%   remaining figures.
% 'freqlim' - rescales each gram to the y-range specified by <freqlim(1:2)>
%   and so the resulting average is restricted to this frequency range. Do
%   not know if it works with 'baseline' option!
% 'log' - converts result to dB.
% 'maxlevels', N - when using "'sem', 'auto'", limits the number of levels
%   marked to <N>.
% 'nodisplay' - To control display visibility. By default the figure will
%   be visible.
%  Specifying 'nodisplay' still creates the figure, but with 'Visible'
%  set to 'off'.
% 'norm' - divides result by its maximum value.
% 'sem', level - adds white contour lines to the average figure show where
%   the standard error of the mean is equal to <level>.  If <level>
%   contains more than one element, then contour lines are drawn at the
%   remaining values in progressively darker shades of gray, the darkest of
%   which is [.3 .3 .3].  If <level> is the string 'auto', then contour
%   levels are set (approximately) logarithmically either at integral log10
%   values or at 10^N multiples of 1, 2 and 5; either way, the first
%   element of the series is chosen so that if there were a smaller element
%   before it, the smaller element would be less than the smallest SEM in
%   the gram. Returns stuff.C and stuff.h, which are cell arrays containing
%   the C and H values returned from the 'contour' function for each
%   contour level and each image, in contours X images form.  Also returns
%   stuff.semlevel, which is the vector of contour levels plotted (useful
%   when <level> = 'auto').  Note that in pasteups, the contour lines are
%   drawn separately for each panel, which means there can be open contour
%   lines at splice points in the middle of the pasteup.
% 'semcolor', color - <color> is three element RGB row vector on a 0 to 1
%   scale; default is [1 1 1].
% 'semwidth', width - <width> is a line width in points; default = 1.
% 'std' - creates an additional gram showing the standard deviation across
%   the input figures.  Variance is computed after normalization (if any).
%   When this option is invoked, the return values <hI> and <hCB> are
%   two-element vectors where the first element refers to the average
%   figure and the second refers to the std figure.

% NOTES
% Event markers and other line objects (if any) are simply copied over from
% the last gram in the list. No attempt is made average placement of event
% markers.

%$Rev: 238 $
%$Date: 2016-03-07 17:46:16 -0500 (Mon, 07 Mar 2016) $
%$Author: dgibson $

stuff = [];
freqlim = [];

argnum = 1;
alignLflag = false;
baselineflag = false;
diffflag = false;
display_str = 'on';
logflag = false;
maxlevels = [];
normflag = false;
semcolor = [1 1 1];
semlevel = [];
semwidth = 1;
stdflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'alignL'
            alignLflag = true;
        case {'baseline' 'diff'}
            if isequal(varargin{argnum}, 'baseline')
                baselineflag = true;
            else
                diffflag = true;
            end
            baselinefig = figlist(1);
            figlist = figlist(2:end);
        case 'freqlim'
            argnum = argnum + 1;
            freqlim = varargin{argnum};
        case 'log'
            logflag = true;
        case 'maxlevels'
            argnum = argnum + 1;
            maxlevels = varargin{argnum};
        case 'nodisplay'
            display_str = 'off';
        case 'norm'
            normflag = true;
        case 'sem'
            argnum = argnum + 1;
            semlevel = varargin{argnum};
        case 'semcolor'
            argnum = argnum + 1;
            semcolor = varargin{argnum};
        case 'semwidth'
            argnum = argnum + 1;
            semwidth = varargin{argnum};
        case 'std'
            stdflag = true;
        otherwise
            error('dg_avgGram:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if (baselineflag || diffflag) && normflag
    error('dg_avgGram:badopt', ...
        'Do not specify ''norm'' and ''baseline'' or ''diff'' together.' );
end

switch(class(figlist))
    case {'double' 'matlab.ui.Figure'}
        hF_in = figlist;
        figname = cell(size(figlist));
        for figidx = 1:length(figlist)
            if isa(figlist, 'matlab.ui.Figure')
                figname{figidx} = sprintf('Figure %d', ...
                    get(hF_in(figidx), 'Number')); 
            else
                figname{figidx} = sprintf('Figure %d', hF_in(figidx)); 
            end
        end
        if baselineflag || diffflag
            hF_baseline = baselinefig;
            if verLessThan('matlab','8.4.0')
                basename = sprintf('Figure %d', baselinefig);
            else
                basename = sprintf('Figure %d', baselinefig.Number);
            end
        end
    case 'cell'
        figname = figlist;
        if baselineflag || diffflag
            basename = baselinefig{1};
        end
    otherwise
        error('Bad argument type, must be double or cell');
end

if baselineflag || diffflag
    if isequal(class(figlist), 'cell')
        open(basename);
        hF_baseline = gcf;
    end
    hI_baseline = findobj(hF_baseline, 'Type', 'image');
    if verLessThan('matlab','8.4.0')
        basedata = get(hI_baseline(2), 'CData');
    else
        basedata = get(hI_baseline, 'CData');
    end
    if isequal(class(figlist), 'cell')
        close(hF_baseline);
    end
end

cdataArray = {};
for figidx = 1:length(figlist)
    if isequal(class(figlist), 'cell')
        open(figlist{figidx});
        hF_in(figidx) = gcf;
    end
    if ~isempty(freqlim)
        hLabels = findobj(gcf,'Type','text');
        for k = 1:length(hLabels);
            set(hLabels,'Units','normalized');
        end
        ylim(freqlim);
    end
    hI_in = findobj(hF_in(figidx), 'Type', 'image');
    % iterate over each panel of pasteup:
    if verLessThan('matlab','8.4.0')
        numimg = length(hI_in)-1;
        imgoffset = 1;
    else
        numimg = length(hI_in);
        imgoffset = 0;
    end
    if numimg > 1 && (baselineflag || diffflag)
        error('dg_avgGram:numimgdiff', ...
            'The ''diff'' and ''baseline'' options cannot be used on pasteups.')
    end
    for imagenum = 1:numimg
        cdata = get(hI_in(imagenum+imgoffset), 'CData');
        if length(cdataArray) < imagenum
            cdataArray{imagenum}(:,:,1) = cdata; %#ok<AGROW>
            if baselineflag || diffflag
                if size(cdata,1) ~= size(basedata,1)
                    if isequal(class(figlist), 'cell')
                        close(hF_in);
                    end
                    error('dg_avgGram:base1', ...
                        '%s does not match vertical size of baseline.', ...
                        figname{figidx} );
                end
                if diffflag && size(cdata,2) ~= size(basedata,2)
                    error('dg_avgGram:base2', ...
                        '%s does not match horizontal size of baseline.', ...
                        figname{figidx} );
                end
            end
        else
            if ~isequal(size(cdata), ...
                    [size(cdataArray{imagenum},1) size(cdataArray{imagenum},2)] )
                if alignLflag && (size(cdata, 1) == size(cdataArray{imagenum},1))
                    % truncate the larger array on the right
                    if size(cdata, 2) > size(cdataArray{imagenum},2)
                        warning('dg_avgGram:alignL1', ...
                            'Clipping %s image %d from %d to %d columns', ...
                            figname{figidx}, imagenum, ...
                            size(cdata, 2), size(cdataArray{imagenum},2) );
                        cdata(:, size(cdataArray{imagenum},2)+1 : end) = [];
                    else
                        warning('dg_avgGram:alignL2', ...
                            'At %s, clipping avg from %d to %d columns', ...
                            figname{figidx}, ...
                            size(cdataArray{imagenum}, 2), size(cdata,2) );
                        cdataArray{imagenum}(:, size(cdata,2)+1 : end, :) = []; %#ok<AGROW>
                    end
                else
                    if isequal(class(figlist), 'cell')
                        close(hF_in);
                    end
                    error('dg_avgGram:size', ...
                        ['%s\nimage #%d of size ' ...
                        '%d X %d does not match size of preceding figures ' ...
                        '(%d X %d)'], figname{figidx},  ...
                        imagenum, size(cdata, 1), size(cdata, 2), ...
                        size(cdataArray{imagenum},1), ...
                        size(cdataArray{imagenum},2) );
                end
            end
            cdataArray{imagenum}(:,:,end+1) = cdata;  %#ok<AGROW>
        end
        % Copy image-specific properties for pasteup from the last fig
        xvals = cell(numimg,1);
        yvals = cell(numimg,1);
        if figidx == length(figlist)
            xvals{imagenum} = get(hI_in(imagenum+imgoffset), 'XData');
            % xvals may need to be clipped to match cdataArray when using
            % alignLflag:
            xvals{imagenum} = ...
                xvals{imagenum}(1:size(cdataArray{imagenum}, 2));
            yvals{imagenum} = get(hI_in(imagenum+imgoffset), 'YData');
        end
    end
    % Figure-wide properties:
    if figidx == length(figlist)
        if verLessThan('matlab','8.4.0')
            cbaraxes = get(hI_in(1), 'Parent');
            gramaxes = get(hI_in(2), 'Parent');
        else
            cbaraxes = findobj(hF_in(figidx), 'Type', 'Colorbar');
            gramaxes = findobj(hF_in(figidx), 'Type', 'Axes');
        end
    end
    if isequal(class(figlist), 'cell') && figidx ~= length(figlist)
        close(hF_in(figidx));
    end
end

fignamedisp = '';
for figidx = 1:length(figname)
    fignamedisp = [ fignamedisp sprintf('disp(''%s'');', figname{figidx}) ]; %#ok<AGROW>
end
% imagenum iterates over each panel of pasteup:
newdata = cell(length(cdataArray),1);
newdatastd = cell(length(cdataArray),1);
for imagenum = 1:length(cdataArray)
    newdata{imagenum} = mean(cdataArray{imagenum},3);
    if normflag
        newdata{imagenum} = newdata{imagenum} / max(max(newdata{imagenum}));
    end
    if stdflag || ~isempty(semlevel)
        newdatastd{imagenum} = std(cdataArray{imagenum},0,3);
    end
    if logflag
        if any(any(newdata{imagenum} < 0))
            warning('dg_avgGram:log', ...
                'Ignoring ''log'' option due to non-positive result');
            logflag = false;
        else
            newdata{imagenum} = 10 * log10(newdata{imagenum});
        end
    end
    if baselineflag || diffflag
        if imagenum == 1
            if baselineflag
                stuff.basespectrum = mean(basedata,2);
            else % 'diff'
                stuff.basespectrum = basedata;
            end
            if logflag
                if any(any(stuff.basespectrum(:) < 0))
                    warning('dg_avgGram:baselog', ...
                        'Ignoring ''log'' option on baseline due to non-positive result');
                else
                    stuff.basespectrum = 10 * log10(stuff.basespectrum);
                end
            end
            if baselineflag
                basegram = repmat(stuff.basespectrum, 1, ...
                    size(newdata{imagenum},2));
            else % 'diff'
                basegram = stuff.basespectrum;
            end
        end
        newdata{imagenum} = newdata{imagenum} - basegram;
    end
end
if ~isempty(varargin)
    optionstr = regexprep(dg_thing2str(varargin{1}), '''', '''''');
        if ismember(dg_thing2str(varargin{1}), {'''baseline''', '''diff'''})
            optionstr = sprintf( '%s, %s', optionstr, ...
                regexprep(dg_thing2str(basename), '''', '''''') );
        end
    for argnum = 2:length(varargin)
        optionstr = sprintf( '%s, %s', optionstr, ...
            regexprep(dg_thing2str(varargin{argnum}), '''', '''''') );
        if ismember(dg_thing2str(varargin{argnum}), {'''baseline''', '''diff'''})
            optionstr = sprintf( '%s, %s', optionstr, ...
                regexprep(dg_thing2str(basename), '''', '''''') );
        end
    end
else
    optionstr = '-no options-';
end
fignamedisp = sprintf('%sdisp(''%s'');', ...
    fignamedisp, optionstr);
% start with a copy of the last input fig, then modify with averaged data.
% Note that the image in <gramaxes> might comprise several panels "pasted"
% together.  <hF> is the output figure, <hA> is the output axes.
hF = figure('Visible', display_str);
if verLessThan('matlab','8.4.0')
    hA = copyobj(gramaxes, hF);
    hCB = copyobj(cbaraxes, hF);
    CBpos = get(cbaraxes, 'Position');
    gramPos = get(hA, 'Position');
    CBpos(4) = gramPos(4);
    set(hCB, 'Position', CBpos);
else
    if iscell(gramaxes)
        for panelidx = 1:length(gramaxes)
            h2copy(panelidx) = gramaxes{panelidx}; %#ok<AGROW>
        end
    else
        h2copy = gramaxes;
    end
    h2copy(end+1) = cbaraxes;
    hAllCopies = copyobj(h2copy, hF);
    hCB = hAllCopies(end);
    hA = hAllCopies(1:end-1);
end
hT = title(hA, 'Click here for input figure list & options');
set(hT, 'ButtonDownFcn', fignamedisp);

if isequal(semlevel, 'auto')
    % Do a quick iteration through all images to get max and min std levels
    maxstd = -Inf;
    minstd = Inf;
    for imagenum = 1:length(newdatastd)
        maxstd = max(max(newdatastd{imagenum}(:)), maxstd);
        minstd = min(min(newdatastd{imagenum}(:)), minstd);
    end
    minsem = minstd/sqrt(length(figlist));
    maxsem = maxstd/sqrt(length(figlist));
    logminsem = log10(minsem);
    logmaxsem = log10(maxsem);
    if logmaxsem - logminsem > 2
        % mark integral decades
        semlevel = logspace(ceil(logminsem), floor(logmaxsem), ...
            floor(logmaxsem) - ceil(logminsem) + 1);
        if length(semlevel) > maxlevels
            semlevel = semlevel(1:maxlevels);
        end
    else
        % mark multiples of [1 2 5]
        [mantissa, exp] = dg_findScale([minsem maxsem], ...
            'levels', 0);
        if length(exp) > 2
            mantissa([1 end]) = [];
            exp([1 end]) = [];
            semlevel = mantissa(end:-1:1) .* 10 .^ exp(end:-1:1);
            if length(semlevel) > maxlevels
                semlevel = semlevel(1:maxlevels);
            end
        else
            warning('dg_avgGram:autosem', ...
                'Range of SEM is too narrow to set scale; using midpoint');
            semlevel = mean([minsem maxsem]);
        end
    end
    stuff.semlevel = semlevel;
end
% <hPanels> is the array of image panels belonging to the output axes:
hPanels = reshape(findobj(hA, 'Type', 'image'), [], 1);
hI = hPanels(1);
for imagenum = 1:length(hPanels)
    set(hPanels(imagenum), 'CData', newdata{imagenum});
    if ~isempty(semlevel)
        xi = linspace(xvals{imagenum}(1), xvals{imagenum}(end), ...
            length(xvals{imagenum})*4);
        yi = linspace(yvals{imagenum}(1), yvals{imagenum}(end), ...
            length(yvals{imagenum})*4)';
        semi = interp2(reshape(xvals{imagenum}, 1, []), ...
            reshape(yvals{imagenum}, [], 1), ...
            newdatastd{imagenum}/sqrt(length(figlist)), xi, yi);
        set(hA, 'NextPlot', 'add');
        if numel(semlevel) > 1
            graylevel = linspace(1, .3, numel(semlevel));
        else
            graylevel = 1;
        end
        oldclim = get(hA, 'CLim');
        for levnum = 1:numel(semlevel)
            [stuff.C{levnum, imagenum}, stuff.h{levnum, imagenum}] = ...
                contour(hA, xi, yi, semi, [1 1] * semlevel(levnum), ...
                'Color', semcolor * graylevel(levnum), ...
                'LineWidth', semwidth);
        end
        set(hA, 'CLim', oldclim)
    end
end
% This must come after contour lines for event markers <hL> to remain
% clickable:
hL = findobj(gramaxes, 'Type', 'line');
copyobj(hL, hA);

if stdflag
    hF(2) = figure;
    hAstd = copyobj(gramaxes, hF(2));
    hPanelsStd = findobj(hAstd, 'Type', 'image');
    hI(2) = hPanelsStd(1);
    hCB(2) = colorbar('peer', hAstd);
    labelstr = get(get(cbaraxes, 'YLabel'), 'String');
    set(get(hCB(2), 'YLabel'), 'String', ['SD of ' labelstr]);
    hT = title(hAstd, 'Click here for input figure list & options');
    set(hT, 'ButtonDownFcn', fignamedisp);
    % imagenum iterates over each panel of pasteup:
    for imagenum = 1:length(hPanelsStd)
        set(hPanelsStd(imagenum), 'CData', newdatastd{imagenum});
    end
    copyobj(hL, hAstd);
end

% There is a huge number of variables holding handles to elements of this:
if isequal(class(figlist), 'cell')
    close(hF_in(figidx));
end

