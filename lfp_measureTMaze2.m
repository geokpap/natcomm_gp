function mazespec2 = lfp_measureTMaze2(x, y, varargin)
% Constructs a lookup table based on the
% distributions of video position coordinates that correspond to certain
% key trial events.  The state of lfp_SelectedTrials and lfp_BadTrials is
% used in the usual way to select trials that contain legitimate data. Only
% the first instance of each event per trial is used. The position for each
% trial is computed as the mean of the five samples centered on the event.
% After calculating the median for a given event, any data points that are
% more than 20 pixels from the median in either dimension are thrown out,
% and the calculations are repeated until all points are within 20 pixels
% of the median.  The outlier points are logged and a warning is issued,
% and they are plotted as asterisks but not included in the results
% returned in <mazespec2>.
%INPUTS
% <x>, <y> are filenums for channels containing video x-position and
%   y-position traces respectively.
%OUTPUTS
% <mazespec2> is a lfp_linearizeRodentTracker2-compatible T-maze layout
%   specification.  Fields: 'evtIDs', 'means' (mean position over trials),
%   'medians' (median position over trials), 'stds' (standard deviation of
%   position over trials), 'n'; all have one row for each evtID whose
%   position was measured, <means> and <stds> have columns [x y].
%OPTIONS
% 'plot' - show Katy-style tracker-with-events plot.  Video tracker traces
%   and the maze outline are drawn by lfp_disp.  StimOn events are
%   plotted in blue, Goal events are plotted in red.  mazespec2.RGoal,
%   mazespec2.StimOn, and mazespec2.LGoal are drawn as magenta lines;
%   mazespec2.X1, mazespec2.X2, mazespec2.Y0, and mazespec2.Y1 as green
%   lines. 
% 'outline' - includes outline of the date-matched calibration in the plot,
%   but suppresses display of figure legend (which gets complicated)
% 'evtIDs', evtIDs - cell vector specifying a list of event IDs to measure;
%   default is the standard set of photobeam events (see 'evtIDs' in code).

%$Rev: 136 $
%$Date: 2010-07-02 23:08:22 -0400 (Fri, 02 Jul 2010) $
%$Author: dgibson $

lfp_declareGlobals;

% The standard set of photobeam events - note that L/R pairs must be kept
% separate:
evtIDs = {13 6 [31 38] 9 14 15 16 7 8 17 18};

outlineflag = false;
plotflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'evtIDs'
            argnum = argnum + 1;
            evtIDs = varargin{argnum};
        case 'plot'
            plotflag = true;
        case 'outline'
            outlineflag = true;
        otherwise
            error('lfp_measureTMaze2:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

[xposn, yposn, outliersx, outliersy, mazespec2.medians] = ...
    lfp_getEvtPosns(x, y, evtIDs);
mazespec2.evtIDs = evtIDs;
[mazespec2.means(:, 1), mazespec2.n] = dg_nanTolerantMean(xposn, 2);
mazespec2.means(:, 2) = dg_nanTolerantMean(yposn, 2);
mazespec2.stds(:, 1) = dg_nanTolerantStd(xposn, 0, 2);
mazespec2.stds(:, 2) = dg_nanTolerantStd(yposn, 0, 2);

if plotflag
    options = {'eye'};
    if ~outlineflag
        options(end+1) = {'trackeronly'};
    end
    hF = lfp_disp(lfp_enabledTrials(find(lfp_SelectedTrials)), [x y], ...
        [], options{:});
    hA = get(hF, 'CurrentAxes');
    xlim(hA, [0 650]);
    ylim(hA, [0 400]);
    axis(hA, 'equal', 'ij');
    grid(hA, 'on');
    title(hA, lfp_SessionNames{1});
    lfp_plotTMazeEvents(hA, evtIDs, ...
        xposn, yposn, outliersx, outliersy, ~outlineflag);
end


