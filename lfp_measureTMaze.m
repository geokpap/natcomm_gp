function mazespec = lfp_measureTMaze(x, y, varargin)
% Measurement of T-maze coordinates and scale is done by examining the
% distributions of video position coordinates that correspond to certain
% key trial events.  The state of lfp_SelectedTrials and lfp_BadTrials is
% used in the usual way to select trials that contain legitimate data.  It
% is assumed that there is only one instance of each event per trial.
%INPUTS
% <x>, <y> are filenums for channels containing video x-position and
%   y-position traces respectively.
%OUTPUTS
% <mazespec> is a lfp_linearizeRodentTracker-compatible T-maze layout
%   specification.  Fields: 'X1', 'X2', 'Y0', 'Y1', 'cmperpixel', 'StimOn'
%   (x coordinate), 'RGoal' (y coordinate), 'LGoal', 'trackwidth' (width of
%   of central (1-2*tailwidth) of actual trace).
%OPTIONS
% 'goalspacing', <goalspacing> - physical distance (cm) from one goal
%   photobeam to the other
% 'plot' - show Katy-style tracker-with-events plot.  Video tracker traces
%   and the maze outline are drawn by lfp_disp.  StimOn events are
%   plotted in blue, Goal events are plotted in red.  mazespec.RGoal,
%   mazespec.StimOn, and mazespec.LGoal are drawn as magenta lines;
%   mazespec.X1, mazespec.X2, mazespec.Y0, and mazespec.Y1 as green lines.
% 'start2stimon', <start2stimon> - physical distance (cm) from front edge
%   of start box to Stim On photobeam
% 'stimon', <StimOn> - <StimOn> is a list of event IDs that will be treated
%   as Stimulus Onset events.  The default value is [21 22 31 38].  For
%   Jianbin data, use "'stimon', 23".
% 'stimon2crossbar', <stimon2crossbar> - physical distance (cm) from the
%   Stim On photobeam to the centerline of the crossbar of the T (the track
%   connecting the two goals)
%
% The default values for all of the physical distances are as measured by
%   Hu Dan and Ledia 06-Nov-2008.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

start2stimon = 48.5;   % centimeters
goalspacing = 66;   % centimeters
stimon2crossbar = 53.85;
StimOn = [21 22 31 38];

plotflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'goalspacing'
            argnum = argnum + 1;
            goalspacing = varargin{argnum};
        case 'plot'
            plotflag = true;
        case 'start2stimon'
            argnum = argnum + 1;
            start2stimon = varargin{argnum};
        case 'stimon'
            argnum = argnum + 1;
            StimOn = varargin{argnum};
        case 'stimon2crossbar'
            argnum = argnum + 1;
            stimon2crossbar = varargin{argnum};
        otherwise
            error('lfp_measureTMaze:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

tailsize = .05;

RGoal = 17;
LGoal = 18;
EventList = {StimOn RGoal LGoal};
maxerrcm = 2;
stimonevts = [];
rgoalevts = [];
lgoalevts = [];

% Collect relevant event indices from trial periods
for trial = lfp_enabledTrials(find(lfp_SelectedTrials))
    trialevents = lfp_Events( lfp_TrialIndex(trial,1) ...
        : lfp_TrialIndex(trial,2), : );
    evtidx = find(ismember(trialevents(:,2), StimOn));
    if ~isempty(evtidx)
        stimonevts(end+1, 1) = evtidx(1) + lfp_TrialIndex(trial,1) - 1;
    end
    evtidx = find(ismember(trialevents(:,2), RGoal));
    if ~isempty(evtidx)
        rgoalevts(end+1, 1) = evtidx(1) + lfp_TrialIndex(trial,1) - 1;
    end
    evtidx = find(ismember(trialevents(:,2), LGoal));
    if ~isempty(evtidx)
        lgoalevts(end+1, 1) = evtidx(1) + lfp_TrialIndex(trial,1) - 1;
    end
end

stimons = lfp_time2index(lfp_Events(stimonevts, 1));
rgoals = lfp_time2index(lfp_Events(rgoalevts, 1));
lgoals = lfp_time2index(lfp_Events(lgoalevts, 1));

% Remove any trials where the data are NaN:
rgoals(isnan(lfp_Samples{y}(rgoals))) = [];
lgoals(isnan(lfp_Samples{y}(lgoals))) = [];


% Expand the numbers of samples by 5 additional samples on each side of
% each event:
stimons = expandpoints(stimons, 5);
rgoals = expandpoints(rgoals, 5);
lgoals = expandpoints(lgoals, 5);

numtailpts = floor(tailsize * length(stimons)) + 1;
if 2*numtailpts >= length(stimons)
    error('lfp_measureTMaze:stimon_x', ...
        'There are not enough stim on points to measure.' );
end
stimon_x = sort(reshape(lfp_Samples{x}(stimons), [], 1));
stimon_x = stimon_x(numtailpts + 1 : end - numtailpts);
stimon_y = sort(reshape(lfp_Samples{y}(stimons), [], 1));
stimon_y = stimon_y(numtailpts + 1 : end - numtailpts);

numtailpts = floor(tailsize * length(rgoals)) + 1;
if 2*numtailpts >= length(rgoals)
    error('lfp_measureTMaze:rgoal_x', ...
        'There are not enough R Goal points to measure.' );
end
rgoal_x = sort(reshape(lfp_Samples{x}(rgoals), [], 1));
rgoal_x = rgoal_x(numtailpts + 1 : end - numtailpts);
rgoal_y = sort(reshape(lfp_Samples{y}(rgoals), [], 1));
rgoal_y = rgoal_y(numtailpts + 1 : end - numtailpts);

numtailpts = floor(tailsize * length(lgoals)) + 1;
if 2*numtailpts >= length(lgoals)
    error('lfp_measureTMaze:lgoal_x', ...
        'There are not enough L Goal points to measure.' );
end
lgoal_x = sort(reshape(lfp_Samples{x}(lgoals), [], 1));
lgoal_x = lgoal_x(numtailpts + 1 : end - numtailpts);
lgoal_y = sort(reshape(lfp_Samples{y}(lgoals), [], 1));
lgoal_y = lgoal_y(numtailpts + 1 : end - numtailpts);

avgrgoal_y = mean(rgoal_y);
avglgoal_y = mean(lgoal_y);
mazespec.RGoal = avgrgoal_y;
mazespec.LGoal = avglgoal_y;
goalspacingpixels = avgrgoal_y - avglgoal_y;
cmperpixel_y = goalspacing / goalspacingpixels;

avgstimon_x = mean(stimon_x);

allgoal_x = [rgoal_x; lgoal_x];
cmperpixel_x = stimon2crossbar / (mean(allgoal_x) - avgstimon_x);
if (max(cmperpixel_x, cmperpixel_y) ...
        / min(cmperpixel_x, cmperpixel_y)) > 1.05
    warning('lfp_measureTMaze:cmperpixel', ...
        'Horizontal (%.2f) and vertical (%.2f) cmperpixel differ by factor of %.2f', ...
        cmperpixel_x, cmperpixel_y, ...
        max(cmperpixel_x, cmperpixel_y) ...
        / min(cmperpixel_x, cmperpixel_y) );
end
mazespec.cmperpixel = (cmperpixel_x + cmperpixel_y) / 2;

maxerrpixels = maxerrcm / mazespec.cmperpixel;
if abs(lgoal_x(1) - rgoal_x(1)) > maxerrpixels
    warning('lfp_measureTMaze:goal_x', ...
        'Estimates of goal X position differ by %.0f pixels (%.1f cm)', ...
        abs(lgoal_x(1) - rgoal_x(1)), ...
        abs(lgoal_x(1) - rgoal_x(1)) * mazespec.cmperpixel );
end

mazespec.StimOn = avgstimon_x;
trackwidth_y = stimon_y(end) - stimon_y(1);
trackwidth_x = max(allgoal_x) - min(allgoal_x);
mazespec.trackwidth = (trackwidth_x + trackwidth_y) / 2;
if abs(trackwidth_y - trackwidth_x) > maxerrpixels
    warning('lfp_measureTMaze:trackwidth', ...
        'Estimates of track width differ by %.0f pixels (%.1f cm)', ...
        abs(trackwidth_y - trackwidth_x), ...
        abs(trackwidth_y - trackwidth_x) * mazespec.cmperpixel );
end
start2StimOnpixels = start2stimon / mazespec.cmperpixel;
mazespec.X1 = avgstimon_x - start2StimOnpixels;
mazespec.X2 = (lgoal_x(1) + rgoal_x(1)) / 2  - mazespec.trackwidth/2;
mazespec.Y0 = stimon_y(1) - mazespec.trackwidth/2;
mazespec.Y1 = stimon_y(end) + mazespec.trackwidth/2;

if plotflag
    hF = lfp_disp(lfp_enabledTrials(find(lfp_SelectedTrials)), [x y], ...
        [], 'eye', 'pass', {'mazespec', mazespec});
    hA = get(hF, 'CurrentAxes');
    hold(hA, 'on');
    xlim(hA, [0 650]);
    ylim(hA, [0 400]);
    axis(hA, 'equal', 'ij');
    grid(hA, 'on');
    title(hA, lfp_SessionNames{1});
    plot(hA, lfp_Samples{x}(stimons), lfp_Samples{y}(stimons), ...
        'b.', 'MarkerSize', 18);
    plot(hA, lfp_Samples{x}(rgoals), lfp_Samples{y}(rgoals), ...
        'r.', 'MarkerSize', 18);
    plot(hA, lfp_Samples{x}(lgoals), lfp_Samples{y}(lgoals), ...
        'r.', 'MarkerSize', 18);
    plot( hA, ...
        [mean(rgoal_x) - 2*mazespec.trackwidth mean(rgoal_x) + ...
        2*mazespec.trackwidth], ...
        [mazespec.RGoal mazespec.RGoal], 'm' );
    plot( hA, ...
        [mean(lgoal_x) - 2*mazespec.trackwidth mean(lgoal_x) + ...
        2*mazespec.trackwidth], ...
        [mazespec.LGoal mazespec.LGoal], 'm' );
    plot( hA, ...
        [mazespec.StimOn mazespec.StimOn], ...
        [mean(stimon_y) - 2*mazespec.trackwidth mean(stimon_y) + ...
        2*mazespec.trackwidth], ...
        'm' );
    plot( hA, ...
        [mazespec.X1 mazespec.X1], ...
        [mean(stimon_y) - 2*mazespec.trackwidth mean(stimon_y) + ...
        2*mazespec.trackwidth], ...
        'g' );
    plot( hA, ...
        [mazespec.X2 mazespec.X2], ...
        [mean(stimon_y) - 2*mazespec.trackwidth mean(stimon_y) + ...
        2*mazespec.trackwidth], ...
        'g' );
    plot( hA, ...
        [mean(rgoal_x) + 2*mazespec.trackwidth mean(rgoal_x) + ...
        4*mazespec.trackwidth], ...
        [mazespec.Y0 mazespec.Y0], 'g' );
    plot( hA, ...
        [mean(rgoal_x) + 2*mazespec.trackwidth mean(rgoal_x) + ...
        4*mazespec.trackwidth], ...
        [mazespec.Y1 mazespec.Y1], 'g' );
end


function result = expandpoints(rawpts, n)
result = reshape( ...
    repmat(rawpts, 1, 2*n + 1) + repmat(-n:n, length(rawpts), 1), ...
    [], 1 );


