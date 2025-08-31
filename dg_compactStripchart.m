function [hF, hA, hL, offset] = dg_compactStripchart(data, offset, varargin)
% Plots each row of <data> as one line, with a downwards vertical offset of
% <offset> for each successive row of <data>, all in a single axes.
%INPUTS
% data: CSC data in channels/trials x samples format.  "Empty" rows can be
%   created by setting all samples in the row to NaN.
% offset: the distance, in the same units as <data>, by which each
%   successive trace is displaced downwards.  If <offset> is empty, then it
%   is set to the smallest value necessary to prevent successive traces
%   from crossing each other.  If <offset> is a vector, it must be of
%   length size(data, 1); in this case, each value is used as the absolute
%   vertical offset on the plot for each channel/trial (i.e. each row of
%   <data>), in the same units as <data>, with increasing offsets moving
%   downwards.
%OUTPUTS
% hF: figure handle
% hA: axes handle
% hL: column vector of handles to the line objects plotted
% offset: the value of <offset> actually used in plotting
%OPTIONS
% 'axes', hA - existing axes handle into which waveforms are to be plotted.
% 'colors', colors - <colors> is a three-column array where each row
%   specifies the RGB color used for the corresponding line plot.
% 'fig', hF - existing figure handle where new axes object is created.
% 'LineWidth', linewidth - sets line width for all traces.  Default is 1.
% 'plotCL', CL - plots symmetrical confidence limits as shaded backgrounds.
%   <CL> is an array the same size as <data> containing values to be added
%   to <data> to yield upper CL and subtracted from <data> for lower CL.
% 'samplesize', sampsz - <sampsz> is used to calibrate the X axis.  Where
%   samples represent time points, <sampsz> is the sample period.  Default
%   value is 1.
% 'signif', plevel, N - when using 'plotCL', doubles the line thickness at
%   points that are statistically significantly nonzero by two-tailed
%   t-test at p < plevel; values submitted as <CL> for 'plotCL' option are
%   interpreted as representing two times the standard error of the mean
%   for this purpose.  <N> is a vector of length size(data, 1) containing
%   the number of observations that went into each row of <data>.  <N> may
%   also be a scalar, in which case it applies to all rows of <data>.  If
%   <N> is not known, the value <Inf> can be used to invoke a z test
%   instead of a t test.  When using this option, <hL> has a second column
%   containing handles to the highlight lines.  If there is only one color
%   specified, then significant points are plotted in red.
% 'x0', x0 - <x0> is the sample number that should be calibrated as zero in
%   the x values plotted.  Default is 1.

%$Rev: 259 $
%$Date: 2017-08-02 14:37:57 -0400 (Wed, 02 Aug 2017) $
%$Author: dgibson $

colors = [];
CL = [];
hA = [];
hF = [];
linewidth = 1;
N = [];
plevel = [];
sampsz = 1;
x0 = 1;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'axes'
            argnum = argnum + 1;
            hA = varargin{argnum};
            hF = get(hA, 'Parent');
        case 'colors'
            argnum = argnum + 1;
            colors = varargin{argnum};
        case 'fig'
            argnum = argnum + 1;
            hF = varargin{argnum};
        case {'LineWidth' 'linewidth'}
            argnum = argnum + 1;
            linewidth = varargin{argnum};
        case 'plotCL'
            argnum = argnum + 1;
            CL = varargin{argnum};
        case 'samplesize'
            argnum = argnum + 1;
            sampsz = varargin{argnum};
        case 'signif'
            argnum = argnum + 1;
            plevel = varargin{argnum};
            argnum = argnum + 1;
            N = varargin{argnum};
            if numel(N) == 1
                N = repmat(N, size(data,1), 1);
            end
        case 'x0'
            argnum = argnum + 1;
            x0 = varargin{argnum};
        otherwise
            error('dg_compactStripchart:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

numrows = size(data,1);
if isempty(hF)
    hF = figure;
end
if isempty(hA)
    hA = axes('Parent', hF, 'NextPlot', 'add');
else
    set(hA, 'NextPlot', 'add');
end
if isempty(colors)
    colors = get(hA, 'ColorOrder');
end
if isempty(offset)
    offset = max(max(abs(diff(data, 1))));
elseif length(offset) > 1 && length(offset) ~= size(data,1)
    error('dg_compactStripchart:offsetlength', ...
        '<offset> must either be a scalar or of length equal to the number of rows in <data>.');
end
if isempty(CL) && ~isempty(plevel)
    error('dg_compactStripchart:plevel', ...
        'The ''signif'' option cannot be specified without using ''plotCL''.');
end
if ~isempty(plevel) && ( numel(plevel) > 1 || isempty(N) || ...
        (numel(N) > 1) && ~isequal(length(N), size(data, 1)) )
    error('dg_compactStripchart:badsignif', ...
        'See header comments regarding values for ''signif'' option.');
end

xvals = ((0:size(data,2)-1) - x0 + 1) * sampsz;
hL = NaN(size(data,1), 1 + ~isempty(plevel));
if ~isempty(plevel)
    SEM = CL/2;
end

% Do the plotting
for rownum = 1:numrows
    if numel(offset) == 1
        vpos = -(rownum - 1) * offset;
    else
        vpos = -offset(rownum);
    end
    if isempty(CL)
        hL(rownum) = plot(hA, xvals, data(rownum, :) + vpos, 'Color', ...
            colors(mod((rownum - 1), size(colors, 1)) + 1, :) );
    else
        CLdata = [ reshape(xvals, [], 1) ...
            data(rownum, :)' - CL(rownum, :)' + vpos' ...
            data(rownum, :)' + CL(rownum, :)' + vpos' ...
            data(rownum, :)' + vpos'];
        [~, hL(rownum)] = dg_plotShadeCL(hA, CLdata, ...
            'Color', colors( ...
            mod((rownum - 1), size(colors, 1)) + 1, : ));
    end
    set(hL(rownum), 'LineWidth', linewidth);
    if ~isempty(plevel)
        % Do stats for significantly non-zero, and add highlighting.
        % The factor of 2 is because we don't know a priori whether <data> will
        % be positive or negative, so we take abs(data) and look at upper tail.
        % The factor of <SEM> is because 'tcdf' assumes normal distribution
        % with sigma = 1.
        p = 2 * tcdf( abs(data(rownum, :) ./ SEM(rownum, :)), ...
            N(rownum), 'upper' );
        issig = p < plevel;
        data2plot = data(rownum, :);
        data2plot(~issig) = NaN;
        if size(colors,1) > 1
            highlightclr = colors(mod((rownum - 1), size(colors, 1)) + 1, :);
        else
            highlightclr = [0.8 0 0];
        end
        hL(rownum,2) = plot(hA, xvals, data2plot + vpos, ...
            'Color', highlightclr, 'LineWidth', 2*linewidth);
    end
end


end