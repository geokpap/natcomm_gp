function hL = dg_plotMazespec(hA, mazespec, varargin)
% Superimposes blue outline of T maze on the plot pointed to by axes <hA>,
% with markers for X1 (blue, front of start box) and StimOn, LGoal, RGoal
% photobeams (red).
%INPUTS
% mazespec - a dg_measureTMaze/lfp_linearizeRodentTracker compatible
%   mazespec.  Required Fields: 'X1', 'Y0', 'Y1', 'cmperpixel'.
%OUTPUT
% hL - a column vector of line handles for all the lines plotted by this
%   function.
%OPTIONS
% 'mazespec2' - treats <mazespec> as a mazespec2 (see lfp_measureTMaze2)
%   and plots a green Start marker in place of the StimOn marker.
%   <mazespec.evtIDs> must contain events 17 (RGoal), 18 (LGoal), and 13
%   (Start), and <mazespec.medians> must contain values for those events.

%$Rev: 49 $
%$Date: 2010-03-26 17:54:40 -0400 (Fri, 26 Mar 2010) $
%$Author: dgibson $

XStimOn = [];
XStart = [];
mazespec2 = [];
hL = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'mazespec2'
            mazespec2 = mazespec;
        otherwise
            error('dg_plotMazespec:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

oldNextPlot = get(hA, 'NextPlot');
set(hA, 'NextPlot', 'add');

overhang = 10;   % number of pixels by which markers exceed track width
x0_x1 = 20.2;    % X0 position to front of start box, cm
x0_2 = 119.5;   % distance from X0 position to X2 position, cm
x0_StimOn = 20.2 + 48.5;   % X0 position to StimOn photobeam, cm
goal_goal = 66;    % distance between goal photobeams, cm
y4_y3 = 73.8;     % distance from Y4 position to Y3 position, cm
x1_start = 13.7;    % distance from X1 position to Start photobeam, cm
trackwidth = 7.7;   % cm
endbit = 3.5;   % distance from Goal photobeam to end of track, cm

% Extract or calculate key parameters from mazespec:
if isempty(mazespec2)
    cmperpixel = mazespec.cmperpixel;
    names = fieldnames(mazespec);
    if ismember('X0', names)
        X0 = mazespec.X0;
    else
        X0 = mazespec.X1 - x0_x1/cmperpixel;
    end
    X1 = mazespec.X1;
    Y0 = mazespec.Y0;
    Y1 = mazespec.Y1;
    
    X2 = X0 + x0_2/cmperpixel;
    X3 = X2 + (Y1 - Y0);    % assume the track is the same width in X and Y
    XStimOn = X0 + x0_StimOn/cmperpixel;
    Y3 = (Y0 + Y1)/2 + y4_y3/(2*cmperpixel);
    Y4 = (Y0 + Y1)/2 - y4_y3/(2*cmperpixel);
    YRGoal = (Y0 + Y1)/2 + goal_goal/(2*cmperpixel);
    YLGoal = (Y0 + Y1)/2 - goal_goal/(2*cmperpixel);
else
    RGoalidx = cellfun(@isequal, ...
        mazespec2.evtIDs, repmat({17}, size(mazespec2.evtIDs)) );
    LGoalidx = cellfun(@isequal, ...
        mazespec2.evtIDs, repmat({18}, size(mazespec2.evtIDs)) );
    Startidx = cellfun(@isequal, ...
        mazespec2.evtIDs, repmat({13}, size(mazespec2.evtIDs)) );
    cmperpixel = goal_goal / ...
        (mazespec2.medians(RGoalidx,2) - mazespec2.medians(LGoalidx,2));
    XStart = mazespec2.medians(Startidx,1);
    X0 = XStart - (x0_x1 + x1_start)/cmperpixel;
    X1 = XStart - x1_start/cmperpixel;
    Y0 = mazespec2.medians(Startidx,2) - trackwidth/(2*cmperpixel);
    Y1 = mazespec2.medians(Startidx,2) + trackwidth/(2*cmperpixel);
    X2 = mazespec2.medians(RGoalidx,1) - trackwidth/(2*cmperpixel);
    X3 = mazespec2.medians(RGoalidx,1) + trackwidth/(2*cmperpixel);
    YRGoal = mazespec2.medians(RGoalidx,2);
    YLGoal = mazespec2.medians(LGoalidx,2);
    Y3 = YRGoal + endbit/cmperpixel;
    Y4 = YLGoal - endbit/cmperpixel;
end

% T maze track outline:
hL(1,1) = plot([X0 X2 X2 X3 X3 X2 X2 X0 X0], ...
    [Y0 Y0 Y4 Y4 Y3 Y3 Y1 Y1 Y0]);
% Front of start box:
hL(end+1,1) = plot([X1 X1], [Y0 Y1]);
% Start photobeam:
if ~isempty(XStart)
    hL(end+1,1) = plot([XStart XStart], [Y0-overhang Y1+overhang], 'g');
end
% StimOn, LGoal, RGoal photobeams:
if ~isempty(XStimOn)
    hL(end+1,1) = plot([XStimOn XStimOn], [Y0-overhang Y1+overhang], 'r');
end
hL(end+1,1) = plot([X2-overhang X3+overhang], [YLGoal YLGoal], 'r');
hL(end+1,1) = plot([X2-overhang X3+overhang], [YRGoal YRGoal], 'r');

set(hA, 'NextPlot', oldNextPlot);