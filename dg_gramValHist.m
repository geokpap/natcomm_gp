function [n, x, nbord, values] = dg_gramValHist(figlist, x, varargin)
%[n, x, nbord, values] = dg_gramValHist(figlist)
%[n, x, nbord, values] = dg_gramValHist(figlist, x)
%   Returns histogram data in the style of Matlab's "hist", describing the
%   maximum values in a list of currently open or previously saved figure
%   files.  The figure files must contain exactly two images, the first of
%   which contains the colorbar and the second of which contains the data.
%   This is the normal result of creating the figure by calling imagesc
%   first and colorbar second.

% INPUTS:
%   figlist - a list of saved figures to average, expressed as a cell
%   string vector of filenames, either absolute or relative to current
%   working directory.
%   x - if given, this should be a list of bin centers as accepted by
%   hist; it can also be [] to invoke the default behavior.
% OUTPUTS:
%   [n, x] -  the frequency counts and the bin locations such that bar(x,n)
%   plots the histogram.
%   nbord - the number of values reported that were extracted from grams
%   where the maximum value occurs within 2 elements of the edge, i.e. in
%   the first, last, next-to-first, or next-to-last column or row.
%   values - the raw list of maximum values.
% OPTIONS:
%   'dir' - <figlist> should be the name of a directory containing figure
%   files; every *.fig file in the directory will be used as input.
%   'draw' - creates a new figure window displaying the bar graph.  The
%   title of the figure is clickable to show a report of the figures that
%   were analyzed, "auto x" (0 = manual x-axis, 1 = automatic x-axis), 
%   nbord, and number of figs.
%   'filefilt', <pattern> - By default, only files that match the <pattern>
%   '*.fig' are included by the 'dir' option; 'filefilt' allows you to
%   specify a different pattern for selecting files.
%   'group' - groups filenames like 'S23Acq01_LFP1ACQ01-T3C1_17.fig' by
%   cluster and event, and uses only the file that contains the highest
%   value of coherence. Specifically, the filename is treated as *LFP<n>*
%   where <n> is a number, and each group consists of all files where only
%   <n> varies.
%   'nobinthresh', nobinthresh - same as 'thresh', except without binning; i.e.
%   clickable details include the count of all values greater than or equal
%   to <nobinthresh> regardless of bin boundaries.
%   'prct', p - clickable details will include the value of the <p>th
%   percentile as computed with interpolation by Matlab's prctile.
%   'thresh', thresh - specify a threshold value <thresh>; the clickable
%   details will include a total of the counts in all bins whose lower
%   limit value is greater than or equal to <thresh>.  It is assumed here
%   that all bins are the same width.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 2
    x = [];
end

argnum = 1;
autox = isempty(x);
dirflag = false;
drawflag = false;
filefilter = '*.fig';
filefiltflag = false;
groupflag = false;
workdir = [];
nobinthresh = [];
prct = [];
thresh = [];
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dir'
            dirflag = true;
        case 'draw'
            drawflag = true;
        case 'filefilt'
            filefiltflag = true;
            argnum = argnum + 1;
            filefilter = varargin{argnum};
        case 'group'
            groupflag = true;
        case 'nobinthresh'
            argnum = argnum + 1;
            nobinthresh = varargin{argnum};
        case 'prct'
            argnum = argnum + 1;
            prct = varargin{argnum};
        case 'thresh'
            argnum = argnum + 1;
            thresh = varargin{argnum};
        otherwise
            error('dg_gramValHist:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
if dirflag
    workdir = figlist;
    files = dir(fullfile(workdir, filefilter));
    figlist = cell(size(files));
    [figlist{:}] = deal(files.name);
    for figidx = 1:length(figlist)
        figlist{figidx} = fullfile(workdir, figlist{figidx});
    end
end

% Ensure that figlist is a column vector:
figlist = reshape(figlist, [], 1);
if isempty(figlist)
    error('dg_gramValHist:emptyfiglist', '<figlist> was empty');
end    

% Sort the file list so that all files of a group are together; package
% each group into a cell so that grouplist becomes a cell vector of
% cell vectors.
if groupflag
    grouplist = {};
    figs2group = figlist;
    while ~isempty(figs2group)
        tokens = regexp(figs2group{1}, ...
            '(.*LFP)([0-9]+)(.*)', 'tokens', 'once');
        if isempty(tokens)
            grouplist = [ grouplist; {figs2group(1)} ];
            figs2group(1) = [];
        else
            figidx2move = [];
            for figidx = 1:length(figs2group)
                if strncmp(figs2group{figidx}, tokens{1}, ...
                        length(tokens{1})) ...
                        && strncmp(...
                        fliplr(figs2group{figidx}), fliplr(tokens{3}), ...
                        length(tokens{3}) )
                    figidx2move = [ figidx2move; figidx ];
                end
            end
            grouplist = [ grouplist; {figs2group(figidx2move)} ];
            figs2group(figidx2move) = [];
        end
    end
else
    % make grouplist with each fig in its own group
    grouplist = cell(size(figlist));
    for figidx = 1:length(figlist)
        grouplist{figidx} = figlist(figidx);
    end
end

nbord = 0;

for groupidx = 1:length(grouplist)
    groupfigs = grouplist{groupidx};
    for figidx = 1:length(groupfigs)
        grpvalues = [];
        open(groupfigs{figidx});
        hF_in = gcf;
        hI_in = findobj(hF_in, 'Type', 'image');
        if numel(hI_in) < 2
            error('dg_gramValHist:badfig', ...
                '"%s" does not contain 2 images', figname{figidx} );
        end
        cdata = get(hI_in(2), 'CData');
        cdata(isinf(cdata)) = NaN;
        close(hF_in);
        [grpvalues(figidx) I] = max(cdata(:));
        [row, col] = ind2sub(size(cdata), I);
        if row < 3 || row > size(cdata,1)-2 ...
                || col < 3 || col > size(cdata,2)-2
            nbord = nbord + 1;
        else
            % check to see if there are any instances of the max value in the
            % "border" regions:
            if any(any(cdata(1:2,:) == grpvalues(figidx))) ...
                    || any(any(cdata(end-1:end,:) == grpvalues(figidx))) ...
                    || any(any(cdata(:,1:2) == grpvalues(figidx))) ...
                    || any(any(cdata(:,end-1:end) == grpvalues(figidx)))
                nbord = nbord + 1;
            end
        end
    end % for figidx
    values(groupidx) = max(grpvalues);
end

if isempty(x)
    [n, x] = hist(values);
else
    [n, x] = hist(values, x);
end

if drawflag
    figure;
    bar(x,n);
    hT = title('Click for details');
    if isempty(workdir)
        detailstr = 'disp(''figlist:'');';
        for figidx = 1:length(figlist)
            detailstr = sprintf('%sdisp(''%s'');', detailstr, figname{figidx});
        end
    else
        detailstr = 'disp(''directory:'');';
        detailstr = sprintf('%sdisp(''%s'');', detailstr, workdir);
        detailstr = sprintf('%sdisp(''filefilter:'');', detailstr);
        detailstr = sprintf('%sdisp(''%s'');', detailstr, filefilter);
    end
    detailstr = sprintf('%sdisp(''auto x: %d'');', detailstr, autox);
    detailstr = sprintf('%sdisp(''nbord: %d'');', detailstr, nbord);
    detailstr = sprintf('%sdisp(''number of figs: %d'');', ...
        detailstr, length(figlist));
    if ~isempty(thresh) && length(x) > 1
        binwidth = x(2) - x(1);
        binedges = x - binwidth/2;
        bins2sum = find(binedges >= thresh);
        binsum = sum(n(bins2sum));
        minbinedge = binedges(bins2sum(1));
        detailstr = sprintf('%sdisp(''actual bin threshold: %d'');', ...
            detailstr, minbinedge );
        detailstr = sprintf('%sdisp(''total counts over thresh: %d'');', ...
            detailstr, binsum);
    end
    if ~isempty(nobinthresh) 
        detailstr = sprintf('%sdisp(''total counts over nobinthresh: %d'');', ...
            detailstr, sum(values >= nobinthresh));
    end
    if ~isempty(prct)
        detailstr = sprintf('%sdisp(''%d percentile: %d'');', ...
            detailstr, prct, prctile(values, prct));
    end
    detailstr = sprintf('%sdisp(''total counts: %d'');', ...
            detailstr, sum(n));
    set(hT, 'ButtonDownFcn', detailstr);
end
