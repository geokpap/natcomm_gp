function lfp_plotBurstAnalysis(highscore, lowscore, highratio, plotdata, ...
    varargin)
%INPUTS
% plotdata:  as returned from lfp_burstAnalysis.
% highscore, lowscore, highratio:  as returned from dg_computeCoactivation.
%OPTIONS
% 'smoothing', smoothing - smooths each plot with a Hanning win of width
%   <2*smoothing+1>.  In order to prevent NaNs from taking over the world,
%   we linearly interpolate any NaN runs before smoothing, then replace
%   them again with NaNs after smoothing.
% 'plotresults', resultlist - plots each result listed in <resultlist> in a
%   separate figure.  <resultlist> is a list of integers 1, 2 and/or 3,
%   where 1=highscore, 2=lowscore, 3=highratio.
% 'ylims', ylims - <ylims> is a 3 X 2 numeric array giving the y limits for
%   one plot on each row, numbered the same as for 'plotresults'.  Use NaN
%   in column 1 for automatic y limits.

%$Rev:  $
%$Date:  $
%$Author:  $

global lfp_SelectedTrials;

resultlist = [];
smoothing = [];
ylims = NaN(3,2);
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'plotresults'
            argnum = argnum + 1;
            resultlist = varargin{argnum};
            if ~all(ismember(resultlist, 1:3))
                error('<resultlist> must contain only integers 1:3');
            end
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

if isempty(plotdata.trials)
    warning('lfp_plotBurstAnalysis:trials', ...
        '<plotdata.trials> was empty, assuming all currently enabled trials');
    plotdata.trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
end

if ~isempty(smoothing)
    sw = hanning(2*smoothing+1);
    swgain = sum(sw);
    for rownum = 1:3
        highscore(rownum,:) = naninterpsmooth( ...
            highscore(rownum,:), sw, swgain, smoothing);
        lowscore(rownum,:) = naninterpsmooth( ...
            lowscore(rownum,:), sw, swgain, smoothing);
        highratio(rownum,:) = naninterpsmooth( ...
            highratio(rownum,:), sw, swgain, smoothing);
    end
end

hA = NaN(3,1);
hF = NaN(3,1);
resultstr = {'highscore', 'lowscore', 'highratio'};
if isempty(resultlist)
    hF = figure;
    for resultidx = 1:3
        hA(resultidx) = subplot(3, 1, resultidx, 'Parent', hF);
    end
else
    for resultidx = reshape(resultlist, 1, [])
        hF(resultidx) = figure;
        hA(resultidx) = axes('Parent', hF(resultidx));
        lfp_createFigTitle(hA(resultidx), 'coactivation', ...
            plotdata.trials, plotdata.win, '', '');
    end
end
for axnum = 1:3
    if isempty(resultlist) || ismember(axnum, resultlist)
        eval(sprintf( ...
            'plot(hA(axnum), plotdata.timepts, %s'')', resultstr{axnum}));
        grid(hA(axnum), 'on');
        ylabel(hA(axnum), resultstr{axnum});
        xlabel(hA(axnum), 'Time, s');
        if ~isempty(resultlist) || axnum == 1
            lfp_createFigTitle(hA(axnum), 'coactivation', ...
                plotdata.trials, plotdata.win, '', '');
        end
        if axnum <= size(ylims,1) && ~isnan(ylims(axnum,1))
            set(hA(axnum), 'YLim', (ylims(axnum,:)));
        end
    end
end

end


function result = naninterpsmooth(x, sw, swgain, smoothing)
[nanidx nonnanidx] = dg_findNaNruns(x);
if ~isempty(nanidx) && nanidx(1) == 1
    x(1:nonnanidx(1)) = x(nonnanidx(1));
    headnans = 1:nonnanidx(1)-1;
    nanidx(1) = [];
    if ~isempty(nonnanidx)
        nonnanidx(1) = [];
    end
else
    headnans = [];
end
if length(nanidx) > length(nonnanidx)
    x(nanidx(end):end) = x(nanidx(end)-1);
    tailnans = nanidx(end):length(x);
    nanidx(end) = [];
else
    tailnans = [];
end
for runidx = 1:length(nanidx)
    x(nanidx(runidx):nonnanidx(runidx)-1) = interp1( ...
        [nanidx(runidx)-1, nonnanidx(runidx)], ...
        x([nanidx(runidx)-1, nonnanidx(runidx)]), ...
        nanidx(runidx):nonnanidx(runidx)-1);
end
s1 = conv(x, sw);
result = s1(smoothing + 1 : end - smoothing) / swgain;
for runidx = 1:length(nanidx)
    result(nanidx(runidx):nonnanidx(runidx)-1) = NaN;
end
result(headnans) = NaN;
result(tailnans) = NaN;

end

