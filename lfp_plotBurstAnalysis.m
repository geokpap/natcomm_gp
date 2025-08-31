function lfp_plotBurstAnalysis(values, plotdata, ...
    varargin)
%INPUTS
% plotdata:  as returned from lfp_burstAnalysis.  Fields used are:
%   align
%   aligns
%   numchannels
%   offsets
%   options
%   timepts
%   titlestr
%   trials
%   win
% values:  as returned from lfp_burstAnalysis.  Fields explicitly
% referenced are:
%   isactive
%   rawcoact
%   thresh
%OPTIONS
% 'Color', color - dg_plotShadeCL, i.e. sets the color of both the CL patch
%    and the median solid line.  'color' (no caps) will also work.
% 'FaceAlpha', alpha - as for dg_plotShadeCL.
% 'FaceColor', facecolor -  as for dg_plotShadeCL.
% 'hA', hA - plots the result into axes <hA>.  Can only be used when
%   plotting a single result.
% 'numchanonly' - limits the second title line to showing the number of
%   channels.
% 'plotresults', resultlist - plots each result listed in <resultlist> in a
%   separate figure.  <resultlist> is a cell vector of strings containing
%   the names of the fields to plot.
% 'popburst', trialnum - ignores <presetvals> and <presetnum>, and instead
%   computes whether the population of channels as a whole was
%   statisticially significantly coactive, i.e. whether rawcoact > thresh,
%   for a given trialnum.  The proper row of <values.rawcoact> is chosen
%   so that the corresponding element of <plotdata.trials> is the specified
%   <trialnum>.  Does not work with 'smoothing'.
% 'preset', presetnum - chooses one of several different collections of
%   values to plot together.  Default = 1.  Has no effect if combined with
%   'plotresults'.
% 'smoothing', smoothing - smooths each plot with a Hanning win of width
%   <2*smoothing+1>.  In order to prevent NaNs from taking over the world,
%   we linearly interpolate any NaN runs before smoothing, then replace
%   them again with NaNs after smoothing.  Note that the number of data
%   points in the smoothed version is reduced by 2*<smoothing> so as not to
%   include any of the zero padding added to the ends during 'conv'.
% 'ylims', ylims - <ylims> is a 3 X 2 numeric array giving the y limits for
%   one plot on each row, numbered the same as for 'plotresults'.  Use NaN
%   in column 1 for automatic y limits.

%$Rev: 253 $
%$Date: 2012-01-26 16:35:16 -0500 (Thu, 26 Jan 2012) $
%$Author: dgibson $

global lfp_SelectedTrials lfp_OrigTrialNums

presetvals = {
    {'highscore', 'lowscore', 'highratio', 'thresh'}
    {'highscore', 'lowscore', 'combinedscore'}
    {'highratio', 'highratio2'}
    {'combinedscore', 'highratio2', 'thresh'}
    {'combinedscore', 'highratio2'}
    {'combinedscore', 'thresh'}
    {'combinedscore'}
    };

axes_handle = [];
numchanonlyflag = false;
presetnum = 1;
plotopts = {};
resultlist = [];
smoothing = [];
trialnum = [];  % if not empty, then 'popburst' was invoked
ylims = NaN(3,2);
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case {'Color' 'color'}
            plotopts = [plotopts varargin(argnum:argnum+1)];
            argnum = argnum + 1;
        case {'FaceAlpha' 'facealpha'}
            plotopts = [plotopts varargin(argnum:argnum+1)];
            argnum = argnum + 1;
        case {'FaceColor' 'facecolor'}
            plotopts = [plotopts varargin(argnum:argnum+1)];
            argnum = argnum + 1;
        case 'hA'
            argnum = argnum + 1;
            axes_handle = varargin{argnum};
            if ~ishandle(axes_handle)
                error('<hA> must be a handle');
            end
        case 'numchanonly'
            numchanonlyflag = true;
        case 'plotresults'
            argnum = argnum + 1;
            resultlist = varargin{argnum};
            if ~iscell(resultlist)
                error('<resultlist> must be a cell array');
            end
        case 'popburst'
            argnum = argnum + 1;
            trialnum = varargin{argnum};
        case 'preset'
            argnum = argnum + 1;
            presetnum = varargin{argnum};
        case 'smoothing'
            argnum = argnum + 1;
            smoothing = varargin{argnum};
        case 'ylims'
            argnum = argnum + 1;
            ylims = varargin{argnum};
        otherwise
            error('lfp_plotBurstAnalysis:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~isempty(resultlist)
    presetvals = {resultlist};
    presetnum = 1;
end

if ~isempty(axes_handle)
    if length(presetvals{presetnum}) > 1
        error('lfp_plotBurstAnalysis:hA', ...
            'The ''hA'' option can only be used to plot a single result.');
    end
end

if isempty(plotdata.trials) && isempty(plotdata.titlestr)
    warning('lfp_plotBurstAnalysis:trials', ...
        '<plotdata.trials> was empty, assuming all currently enabled trials');
    plotdata.trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
end

names = fieldnames(plotdata);
ispasteup = ismember('aligns', names);
resultstr = presetvals{presetnum};

if ~isempty(trialnum)
    % 'popburst' was invoked
    if ~isempty(smoothing)
        % This is an error because allowing it would make the code
        % complicated and is of dubious value (rawcoact, thresh, and
        % isactive would all have to be smoothed, and what would the
        % results mean anyway???) 
        error('lfp_plotBurstAnalysis:badcombo', ...
            '''smoothing'' is not available for ''popburst''');
    end
    resultstr = {'rawcoact' 'isactive'};
    rowidx = find(plotdata.trials == trialnum);
    if isempty(rowidx)
        error('Trialnum %d is not a member of plotdata.trials.', trialnum);
    end
    values.isactive = values.rawcoact(rowidx,:) > values.thresh(1,:);
    if isempty(lfp_OrigTrialNums)
        trialstr = sprintf('trialnum %d', trialnum);
    else
        trialstr = lfp_getTrialID(trialnum);
    end
    plotdata.titlestr = sprintf('%s %s', plotdata.titlestr, trialstr);
end


if ~isempty(smoothing)
    sw = hanning(2*smoothing+1);
    for fieldnum = 1:length(resultstr)
        if isfield(values, resultstr{fieldnum})
            for rownum = 1:size(values.(resultstr{fieldnum}),1)
                values.(resultstr{fieldnum})(rownum,:) = ...
                    dg_naninterpSmooth( ...
                    values.(resultstr{fieldnum})(rownum,:), sw );
            end
            values.(resultstr{fieldnum}) = values.(resultstr{fieldnum})...
                (:, smoothing + 1 : end - smoothing);
        end
    end
    plotdata.timepts = plotdata.timepts(smoothing + 1 : end - smoothing);
end

if isempty(axes_handle)
    hA = NaN(length(presetvals{presetnum}),1);
    hF = figure;
end
if isempty(plotdata.titlestr)
    trialsval = plotdata.trials;
else
    trialsval = plotdata.titlestr;
    if any(cellfun( @isequal, ...
            repmat({'onset'}, size(plotdata.options)), plotdata.options ))
        trialsval = sprintf('%s onset', trialsval);
    elseif any(cellfun( @isequal, ...
            repmat({'offset'}, size(plotdata.options)), plotdata.options ))
        trialsval = sprintf('%s offset', trialsval);
    else
        trialsval = sprintf('%s sustained', trialsval);
    end
end
if ~isfield(plotdata, 'aligns') || isempty(plotdata.aligns)
    alignval = plotdata.align;
else
    alignval = plotdata.aligns;
end
if isfield(plotdata, 'numchannels')
    numchstr = sprintf('numchannels=%d', plotdata.numchannels);
else
    numchstr = 'numchannels not available';
end
for resultidx = 1:length(resultstr)
    if isempty(axes_handle)
        hA(resultidx) = subplot(length(resultstr), 1, resultidx, ...
            'Parent', hF);
    else
        hA(resultidx) = axes_handle;
    end
end
for axnum = 1:length(resultstr)
    if isequal(resultstr{axnum}, 'rawcoact')
        % This is the weird one: we want rawcoact(rowidx,:)
        % superimposed on thresh.
        plot(hA(axnum), plotdata.timepts, values.rawcoact(rowidx,:));
        set(hA(axnum), 'NextPlot', 'add');
        plot(hA(axnum), plotdata.timepts, values.thresh(1,:), 'r');
    else
        if isfield(values, resultstr{axnum})
            if size(values.(resultstr{axnum}), 1) < 3
                plot(hA(axnum), plotdata.timepts, ...
                    values.(resultstr{axnum})(1,:));
            else
                dg_plotShadeCL(hA(axnum), ...
                    [ plotdata.timepts' ...
                    values.(resultstr{axnum})([3 2 1],:)' ], ...
                    plotopts{:});
            end
        end
    end
    grid(hA(axnum), 'on');
    ylabel(hA(axnum), resultstr{axnum});
    xlabel(hA(axnum), 'Time, s');
    if axnum == 1
        if numchanonlyflag
            lfp_createFigTitle(hA(axnum), '', ...
                trialsval, [], '', '', 'append', numchstr);
        else
            lfp_createFigTitle(hA(axnum), 'coactivation', ...
                trialsval, plotdata.win, '', '', ...
                'alignment', alignval, 'append', numchstr);
        end
    end
    if axnum <= size(ylims,1) && ~isnan(ylims(axnum,1))
        set(hA(axnum), 'YLim', (ylims(axnum,:)));
    end
    if ispasteup
        offsets = [0 reshape(plotdata.offsets, 1, [])];
        if length(plotdata.aligns) ~= length(offsets)
            warning('lfp_plotBurstAnalysis:offsets', ...
                'There are different numbers of offsets and aligns.');
        end
        numevts2plot = min(length(offsets), length(plotdata.aligns));
        evtIDs = zeros(numevts2plot,1);
        for ididx = 1:numevts2plot
            if ~isempty(plotdata.aligns{ididx})
                evtIDs(ididx) = plotdata.aligns{ididx}(1);
            end
        end
        evtmatrix = [ offsets(1:numevts2plot)' evtIDs ];
        lfp_plotEvtMarkers(hA(axnum), evtmatrix);
    end
end

end


