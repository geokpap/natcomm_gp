function [hF, hA] = lfp_plot(plotdata, result)
% Plots a result in the style to which lfp_lib users should get accustomed.
%INPUTS
% plotdata:  a structure with the following fields, all of which are
%       optional:
%   align - a numeric or cell array specifying the alignment event(s)
%   append - becomes the value to the 'append' option to
%       lfp_createFigTitle; i.e. gets appended to the second line of the
%       title.
%   cbarlabel - see <mode> 'gram'
%   clickstr - as for lfp_createFigTitle, a string containing text to
%       display in the command window when the title is clicked
%   CLstyle - 'line' or 'shade'; defaults to 'line'.
%	CLtype - determines how to plot CL data; possible values are:
%      'symmetric' - plots symmetric CLs at <result(:,3)> above and below
%           <result(:,2)>; note that negative values in <result(:,3)> will
%           produce screwy-looking plots.
%      'nonsymmetric' - plots independent CLs at <result(:,3)> and
%           <result(:,4)>; <result(:,3)> should be the lower CL and
%           <result(:,4)> should be the upper CL when CLstyle is 'shade'.
%   extraline - as for lfp_createFigTitle, a string appended as an extra
%       line to the title
%   figtype - as for lfp_createFigTitle,  a string that identifies what
%       type of analysis was done
%   filenames - string or cell array of strings specifying the
%       channel(s) that was(were) analyzed.  Up to two will be listed after
%       after <sessionstr>; more than two will be appended to <clickstr>.
%   hA - axes handle into which to plot; if empty or missing, a new
%       figure window with a new axes is constructed.
%   mode - a string that determines overall plot and/or <result> format; if
%           empty or missing, 'standard' mode is used.  Possible values
%           are:
%       'standard' - <result> contains one row per sample; col 1 =
%           sample time in seconds, col 2 = sample value, col 3 (optional)
%           and col 4 (optional) contain confidence limit values that are
%           interpreted according to the value of the 'CLtype' field.
%       'gram' - result is a matrix to plot as a pseudocolor image, but in
%           Cartesian (xy) coordinates.  If there are 'xvals' and/or
%           'yvals' fields in <plotdata>, they are used by imagesc to make
%           tick mark labels; if not, then row and column indices are used.
%           If there is a 'cbarlabel' field, it is used as a label for the
%           colorbar; otherwise there is no label.
%   ntrigs - number of "triggers", i.e. number of alignment reference
%       events found (e.g. by lfp_getSamples)
%   plotopts - optional arguments suitable for submission to the Matlab
%       'plot' function when CLstyle is 'line', or to dg_plotShadeCL when
%       CLstyle is 'shade'
%   sessionstr - as for lfp_createFigTitle 'sessionstr' option, a string
%       that overrides the current value of lfp_SessionNames as the source
%       of <sessionstr>, which is the first item on the second line of the
%       title.  <sessionstr> may be '', in which case there is no session
%       label included in the title.
%   trials - as for lfp_createFigTitle, a list of trial numbers or a string
%       describing the trial list
%   win - as for lfp_createFigTitle, the time window relative to
%       lfp_AlignmentRef that was analyzed
%   xlab - a string containing a label for the x axis
%   xval - see <mode> 'gram'
%   ylab - a string containing a label for the y axis
%   yval - see <mode> 'gram'
% result:  data to be plotted; see <plotdata.mode> comments
%OUTPUT
% hF:  handle to the figure window where the plot was made

%$Rev: 408 $
%$Date: 2020-04-20 15:24:17 -0400 (Mon, 20 Apr 2020) $
%$Author: dgibson $

global lfp_CLimAll

appendflag = false;
appendstr = '';

titleopts = {}; % cell row vector

if isfield(plotdata, 'hA') && ~isempty(plotdata.hA)
    hF = get(plotdata.hA, 'Parent');
else
    hF = figure;
    hA = axes('Parent', hF);
end
if isfield(plotdata, 'align') && ~isempty(plotdata.align)
    titleopts = [titleopts, {'alignment', plotdata.align}];
end
if isfield(plotdata, 'sessionstr')
    titleopts = [titleopts, {'sessionstr', plotdata.sessionstr}];
end
if isfield(plotdata, 'append')
    appendflag = true;
    appendstr = [appendstr ' ' plotdata.append];
end
if isfield(plotdata, 'mode') && ~isempty(plotdata.mode)
    mode = plotdata.mode;
else
    mode = 'standard';
end
if isfield(plotdata, 'ntrigs')
    appendflag = true;
    appendstr = sprintf('%s n=%s', ...
        appendstr, dg_thing2str(plotdata.ntrigs) );
end
if isfield(plotdata, 'figtype') && ~isempty(plotdata.figtype)
    figtype = plotdata.figtype;
else
    figtype = 'Mystery Analysis';
end
if isfield(plotdata, 'win')
    win = plotdata.win;
else
    win = NaN;
end
if isfield(plotdata, 'extraline')
    extraline = plotdata.extraline;
else
    extraline = '';
end
if isfield(plotdata, 'clickstr')
    clickstr = plotdata.clickstr;
else
    clickstr = '';
end
if isfield(plotdata, 'CLtype')
    CLtype = plotdata.CLtype;
else
    CLtype = '';
end
if isfield(plotdata, 'CLstyle') && ~isempty(plotdata.CLstyle)
    CLstyle = plotdata.CLstyle;
else
    CLstyle = 'line';
end
if isfield(plotdata, 'filenames')
    if ischar(plotdata.filenames) || length(plotdata.filenames) < 3
        titleopts = [titleopts, {'filenames', plotdata.filenames}];
    else
        titleopts = [titleopts, {'filenames', 'Click for filenames'}];
        for k = 1:length(plotdata.filenames)
            if isempty(clickstr)
                clickstr = 'Filenames:';
            end
            clickstr = sprintf('%s\n%s', ...
                clickstr, plotdata.filenames{k});
        end
    end
end
if isfield(plotdata, 'plotopts')
    plotopts = plotdata.plotopts;
else
    plotopts = {};
end
if isfield(plotdata, 'trials')
    trials = plotdata.trials;
else
    trials = lfp_enabledTrials;
end
if isfield(plotdata, 'xlab')
    xlab = plotdata.xlab;
else
    xlab = 'Mystery Parameter';
end
if isfield(plotdata, 'ylab')
    ylab = plotdata.ylab;
else
    ylab = 'Mystery Parameter';
end

if appendflag
    titleopts = [titleopts, {'append', appendstr}];
end

switch mode
    case 'standard'
        if isempty(CLtype)
            plot(hA, result(:,1), result(:,2), plotopts{:});
        else
            switch CLstyle
                case 'line'
                    switch CLtype
                        case 'symmetric'
                            plot(hA, result(:,1), [ result(:,2) ...
                                result(:,2)-result(:,3) ...
                                result(:,2)+result(:,3) ], plotopts{:});
                        case 'nonsymmetric'
                            plot(hA, result(:,1), [ result(:,2) ...
                                result(:,3:4) ], plotopts{:});
                    end
                case 'shade'
                    switch CLtype
                        case 'symmetric'
                            dg_plotShadeCL(hA, [result(:,1) ...
                                result(:,2)-result(:,3) ...
                                result(:,2)+result(:,3) ...
                                result(:,2) ], plotopts{:});
                        case 'nonsymmetric'
                            dg_plotShadeCL(hA, [result(:,1) ...
                                result(:,3:4) result(:,2) ], plotopts{:});
                    end
                otherwise
                    error('lfp_plot:CLstyle', ...
                        'No such CLstyle as "%s"', dg_thing2str(CLstyle));
            end
        end
        set(hA,'XGrid', 'on', 'YGrid', 'on');
    case 'gram'
        if isfield(plotdata, 'xval') && ~isempty(plotdata.xval)
            xvals = plotdata.xval;
        else
            xvals = 1:size(result,2);
        end
        if isfield(plotdata, 'yval') && ~isempty(plotdata.yval)
            yvals = plotdata.yval;
        else
            yvals = 1:size(result,1);
        end
        if isfield(plotdata, 'cbarlabel') && ~isempty(plotdata.cbarlabel)
            cbarlabel = plotdata.cbarlabel;
        else
            cbarlabel = '';
        end
        [hI, hCB] = dg_showGram(hA, xvals, yvals, result, '', ...
            xlab, ylab, cbarlabel);
        if ~isempty(lfp_CLimAll)
            hCB = dg_recolorGram(hCB, lfp_CLimAll, hI);
        end
    otherwise
        error('lfp_plot:mode', ...
            'Unrecognized mode: %s', mode);
end
lfp_createFigTitle(hA, figtype, trials, win, ...
    extraline, clickstr, titleopts{:});
xlabel(hA, xlab);
ylabel(hA, ylab);

