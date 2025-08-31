function data = lfp_plot_eye_rodent(filenums, trials, startsample, endsample, ...
    varargin)
%OPTIONS
% 'autocolor'
% 'evtIDs', evtIDs - cell vector specifying a list of event IDs to measure;
%   default is the standard set of photobeam events (see 'evtIDs' in code).
%   Used only with 'manualcalib'.  Interacts with 'midT' in a way that
%   cannot be described other than in the code.
% 'mazespec', mazespec - same syntax and data structure accepted by
%   lfp_linearizeRodentTracker; overrides value looked up based on file
%   dates.  If mazespec.X0 is provided, then it overrides the value
%   normally calculated from X1 and cmperpixel.  Required fields in
%   <mazespec> are as for dg_plotMazespec.
% 'manualcalib' - plots data without T maze outline and provides GUI
%   tools for calibrating the tracker and realigning slipped video.  The
%   buttons labeled "Set Stem" through "Set Start" are used to set cursors
%   that specify the locations of the corresponding T maze landmarks.
%   "Outline" plots an outline of the T maze with goal and start markers to
%   fit the cursor positions as well as possible.  "Events" toggles on and
%   off a display of the positions recorded on the video tracker at the
%   times of the photobeam events.  "Realign" invokes an algorithm to
%   estimate video timestamp slippage and realign the timestamps to
%   compensate.  If the "Use Original File" box is checked, then the
%   original vt1.nvt or vt1.dat file is re-read from lfp_DataDir to get the
%   raw timestamps and coordinates.  If the box is not checked, then a new
%   series of timestamps is created and the position data are read from the
%   same CSC channels that were used to create the tracker plot.  This
%   makes it possible to use processed (e.g. smoothed) video data in place
%   of the raw data, but note that there will not be a strict
%   correspondence between the new timestamps and the old timestamps.  The
%   "Tol" slider controls how many frames of slippage will be tolerated
%   without making adjustments.  It is recommended that this be set no
%   lower than 2.0 when making initial adjustments based on a manual video
%   calibration, since a rat traveling at 75 cm/s will spend about 3 frames
%   traversing a distance (usually about 8 cm) that corresponds to 25
%   pixels, and 20 - 25 pixels is about as accurate as the photobeam-video
%   reconciliation gets. 
%     The typical workflow is to click buttons in order from left to right,
%   but the only requirements are that all cursors be set before plotting
%   the outline, and that the outline be plotted before realigning the
%   video timestamps.
% 'midT', midT - specifies a set of events to use as the mid-T photobeam
%   marker; default is [31 38].  Used only with 'manualcalib'.

%$Rev: 333 $
%$Date: 2014-10-16 18:56:27 -0400 (Thu, 16 Oct 2014) $
%$Author: dgibson $

hF = gcf;
hA = axes('Box', 'on');

lfp_declareGlobals;

if ~isequal(size(startsample), size(endsample)) || size(startsample,2) ~= 1
    error('lfp_plot_eye_rodent:samples', ...
        'startsample, endsample must be column vectors of same length');
end

data = [];
autocolorflag = false;
evtIDs = [];
manualcalibflag = false;
mazespec = [];
midT = [31 38];
numtraces = size(startsample,1);
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'autocolor'
            autocolorflag = true;
        case 'evtIDs'
            argnum = argnum + 1;
            evtIDs = varargin{argnum};
        case 'manualcalib'
            manualcalibflag = true;
        case 'mazespec'
            argnum = argnum + 1;
            mazespec = varargin{argnum};
        case 'midT'
            argnum = argnum + 1;
            midT = varargin{argnum};
        otherwise
            error('lfp_plot_eye_rodent:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
if isempty(evtIDs)
    evtIDs = {13 6 midT 9 14 15 16 7 8 17 18};
end

if manualcalibflag
    sessiondatestr = 'manual-calib';
    mazespec.calibstr = 'manual-calib';
    if ~isempty(lfp_XLimAll)
        warning('lfp_plot_eye_rodent:xlimall', ...
            'lfp_XLimAll is not empty, which is usually not good for calibrating tracker data.');
    end
else
    if isempty(mazespec)
        % Find the mazespec based on session date and calib dates
        [mazespec, sessiondatestr] = lfp_findMazespec;
    else
        sessiondatestr = 'user-mazespec';
        names = fieldnames(mazespec);
        if ~ismember('calibstr', names)
            mazespec.calibstr = 'user-mazespec';
        end
    end
end

set(hA, 'NextPlot', 'add');
buttXmarg = 5;
buttXposn = 10;
buttYposn = 10;
buttHeight = 20;
if manualcalibflag
    set(hA, 'ButtonDownFcn', @lfp_plot_eye_rodent_manualcalib);
    buttWidth = 55;
    uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Set Stem', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_stem);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 80;
    uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Set Crossbar', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_xbar);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 60;
    uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Set Goals', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_goal);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 55;
    uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Set Start', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_start);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 40;
    uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Outline', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_outline);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 40;
    uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Events', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_events);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 50;
    uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Realign', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_realign);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 15;
    uicontrol(hF, ...
        'Style', 'checkbox', ...
        'Value', 1, ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Tag', 'lfp_plot_eye_rodent_origfile_checkbox', ...
        'Callback', @lfp_plot_eye_rodent_origfile);
    uicontrol(hF, ...
        'Style', 'text', ...
        'String', 'Use Original File', ...
        'Position', [buttXposn-30 buttYposn+20 80 15]);
    buttXposn = buttXposn + buttWidth + buttXmarg;
    buttWidth = 100;
    uicontrol(hF, ...
        'Style', 'slider', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight/2], ...
        'Tag', 'lfp_plot_eye_rodent_tol_slider', ...
        'Max', [5], 'Min', [0], ...
        'SliderStep', [0.1 0.2], ...
        'Value', 2, ...
        'Callback', @lfp_plot_eye_rodent_tol);
    uicontrol(hF, ...
        'Style', 'text', ...
        'String', 'Tol = 2.0', ...
        'Tag', 'lfp_plot_eye_rodent_tol_text', ...
        'Position', [buttXposn+30 buttYposn+12 50 10]);
    figposn = get(hF, 'Position');
    buttXposn = 10;
    buttHeight = 30;
    buttYposn = figposn(4) - buttHeight;
    buttWidth = 100;
    huic = uicontrol(hF, ...
        'Style', 'pushbutton', ...
        'String', 'Save Outline', ...
        'Position', [buttXposn buttYposn buttWidth buttHeight], ...
        'Callback', @lfp_plot_eye_rodent_savemazespec);
    set(huic, 'Units', 'normalized');
    mydata.stem = NaN;
    mydata.xbar = NaN;
    mydata.LGoal = NaN;
    mydata.RGoal = NaN;
    mydata.start = NaN;
    mydata.hstem = NaN;
    mydata.hxbar = NaN;
    mydata.hLGoal = NaN;
    mydata.hRGoal = NaN;
    mydata.hstart = NaN;
    mydata.hevents = NaN;
    mydata.houtline = NaN;
    mydata.filenums = filenums;
    mydata.evtIDs = evtIDs;
    mydata.midT = midT;
    mydata.trials = trials;
    mydata.ms2 = [];
    if exist(fullfile(lfp_DataDir, 'mazespec2.mat'), 'file')
        load(fullfile(lfp_DataDir, 'mazespec2.mat'));
        mydata.ms2 = mazespec2;
        mydata.houtline = dg_plotMazespec(hA, mydata.ms2, 'mazespec2');
    end
    % IMPORTANT REVISION NOTE: the following line must go last
    set(hA, 'UserData', mydata);
else
    dg_plotMazespec(hA, mazespec);
end

% set plot properties
axis equal;
xlim([0 750]);
ylim([0 450]);
axis ij;

smoothlength=5;
sampleskip=12;
if autocolorflag
    tracecolortable = get(gca, 'ColorOrder');
end
for tracenum = 1:numtraces
    calibrated_samples1c = smooth(lfp_Samples{filenums(1)}( ...
        startsample(tracenum):sampleskip:endsample(tracenum)),smoothlength);
    calibrated_samples2c = smooth(lfp_Samples{filenums(2)}( ...
        startsample(tracenum):sampleskip:endsample(tracenum)),smoothlength);
    if autocolorflag
        tracecolor = tracecolortable( ...
            mod(tracenum-1, size(tracecolortable,1)) + 1, : );
    else
        tracecolor = 'k';
    end
    hL = plot(hA, calibrated_samples1c, calibrated_samples2c, '.', ...
        'Color', tracecolor, 'MarkerSize', 4);
    if manualcalibflag
        set(hL, 'ButtonDownFcn', @lfp_plot_eye_rodent_line);
    end
end

trialstyle = lfp_TrialStyle;
sessionstr = lfp_SessionNames{1};
if numtraces > 1
    titlestr = [ sprintf('%s align=%s session: %s\n ', ...
        sessionstr, mat2str(lfp_AlignmentRef), sessiondatestr) ...
        lfp_getTrialsLabel(trials, trialstyle) ...
        sprintf(' n=%d calib: %s', length(lfp_enabledTrials(trials))), ...
        mazespec.calibstr ];
    title(titlestr, 'Interpreter', 'none');
else
    titlestr = sprintf('%s align=%s trial %s', ...
        sessionstr, mat2str(lfp_AlignmentRef), ...
        lfp_getTrialID(trials) );
    title(titlestr, 'Interpreter', 'none');
end
ylabel(lfp_FileNames{2});
xlabel(lfp_FileNames{1});

end


function lfp_plot_eye_rodent_manualcalib(hA,evnt)
% The 'UserData' prop of the fig window contains the
% identity of the most recently clicked cursor pushbutton.
% Other interaction data is stored in the 'UserData' prop of the axes.
cursorcolor = [.7 .7 .7];
cp = get(hA,'CurrentPoint');    % x=cp(1,1); y=cp(1,2);
hF = get(hA, 'Parent');
mydata = get(hA, 'UserData');
mymode = get(hF, 'UserData');
if ~isempty(mymode)
    switch mymode
        case 'stem'
            mydata.stem = cp(1,2);
            if isnan(mydata.hstem)
                mydata.hstem = plot(hA, get(hA, 'XLim'), ...
                    [mydata.stem mydata.stem], ...
                    'Color', cursorcolor);
                set(mydata.hstem, ...
                    'ButtonDownFcn', @lfp_plot_eye_rodent_line);
            else
                set(mydata.hstem, 'YData', [mydata.stem mydata.stem]);
            end
        case 'xbar'
            mydata.xbar = cp(1,1);
            if isnan(mydata.hxbar)
                mydata.hxbar = plot(hA, [mydata.xbar mydata.xbar], ...
                    get(hA, 'YLim'), ...
                    'Color', cursorcolor);
                set(mydata.hxbar, ...
                    'ButtonDownFcn', @lfp_plot_eye_rodent_line);
            else
                set(mydata.hxbar, 'XData', [mydata.xbar mydata.xbar]);
            end
        case 'goal'
            if cp(1,2) < mydata.stem
                mydata.LGoal = cp(1,2);
                if isnan(mydata.hLGoal)
                    mydata.hLGoal = plot(hA, get(hA, 'XLim'), ...
                        [mydata.LGoal mydata.LGoal], ...
                        'Color', cursorcolor);
                    set(mydata.hLGoal, ...
                        'ButtonDownFcn', @lfp_plot_eye_rodent_line);
                else
                    set(mydata.hLGoal, 'YData', [mydata.LGoal mydata.LGoal]);
                end
            else
                mydata.RGoal = cp(1,2);
                if isnan(mydata.hRGoal)
                    mydata.hRGoal = plot(hA, get(hA, 'XLim'), ...
                        [mydata.RGoal mydata.RGoal], ...
                        'Color', cursorcolor);
                    set(mydata.hRGoal, ...
                        'ButtonDownFcn', @lfp_plot_eye_rodent_line);
                else
                    set(mydata.hRGoal, 'YData', [mydata.RGoal mydata.RGoal]);
                end
            end
        case 'start'
            mydata.start = cp(1,1);
            if isnan(mydata.hstart)
                mydata.hstart = plot(hA, [mydata.start mydata.start], ...
                    get(hA, 'YLim'), ...
                    'Color', cursorcolor);
                set(mydata.hstart, ...
                        'ButtonDownFcn', @lfp_plot_eye_rodent_line);
            else
                set(mydata.hstart, 'XData', [mydata.start mydata.start]);
            end
        otherwise
            error('Burma!');
    end
    set(hA, 'UserData', mydata);
end
end


function lfp_plot_eye_rodent_stem(hObject,evnt)
set(get(hObject, 'Parent'), 'UserData', 'stem');
end


function lfp_plot_eye_rodent_xbar(hObject,evnt)
set(get(hObject, 'Parent'), 'UserData', 'xbar');
end


function lfp_plot_eye_rodent_goal(hObject,evnt)
set(get(hObject, 'Parent'), 'UserData', 'goal');
end


function lfp_plot_eye_rodent_line(hObject,evnt)
hA = get(hObject, 'Parent');
lfp_plot_eye_rodent_manualcalib(hA,evnt);
end


function lfp_plot_eye_rodent_events(hObject,evnt)
hF = get(hObject, 'Parent');
hA = findobj(hF, 'Type', 'axes');
hA = hA(end);   % legend lives in a new axes, we want the original
mydata = get(hA, 'UserData');
[xposn, yposn, outliersx, outliersy] = ...
    lfp_getEvtPosns( mydata.filenums(1), mydata.filenums(2), ...
    mydata.evtIDs, 'trials', mydata.trials);
if ~isnan(mydata.hevents)
    delete(mydata.hevents);
    delete(mydata.hlegend);
    mydata.hevents = NaN;
else
    lfp_declareGlobals;
    mydata.hevents = lfp_plotTMazeEvents(hA, mydata.evtIDs, ...
        xposn, yposn, outliersx, outliersy, false);
    nonoutliers = true(size(mydata.hevents));
    for k = 1:length(mydata.hevents)
        set(mydata.hevents(k), ...
            'ButtonDownFcn', @lfp_plot_eye_rodent_line);
        if isequal(get(mydata.hevents(k), 'Marker'), '*')
            nonoutliers(k) = false;
        end
    end
    legendstr = cell(size(mydata.evtIDs));
    for k = 1:length(mydata.evtIDs)
        if length(mydata.evtIDs{k}) == 0 || length(mydata.evtIDs{k}) > 1
            legendstr{k} = mat2str(mydata.evtIDs{k});
        else
            if isempty(lfp_EventNames{mydata.evtIDs{k}})
                legendstr{k} = sprintf('evt %d', mydata.evtIDs{k});
            else
                legendstr(k) = lfp_EventNames(mydata.evtIDs{k});
            end
        end
    end
    mydata.hlegend = legend(mydata.hevents(nonoutliers), legendstr, ...
        'Location', 'NorthWest');
end
set(hA, 'UserData', mydata);
end


function lfp_plot_eye_rodent_start(hObject,evnt)
set(get(hObject, 'Parent'), 'UserData', 'start');
end


function lfp_plot_eye_rodent_outline(hObject,evnt)
hF = get(hObject, 'Parent');
hA = findobj(hF, 'Type', 'axes');
hA = hA(end);   % legend lives in a new axes, we want the original
mydata = get(hA, 'UserData');
canonicalms2.medians = [
    -1        0
    -0.563    0
    -0.138    0
    0         0.349
    0        -0.349
    0         1
    0        -1
    ];
canonicalms2.evtIDs = { 13 mydata.midT 14 15 16 17 18 };
if ~isnan(mydata.LGoal) ...
        && ~isnan(mydata.RGoal) ...
        && ~isnan(mydata.xbar) ...
        && ~isnan(mydata.stem) ...
        && ~isnan(mydata.start)
    % Fit canonicalms2 to UserData
    xscale = (mydata.xbar - mydata.start);
    yscale = (mydata.RGoal - mydata.LGoal)/2;
    fittedms2.evtIDs = canonicalms2.evtIDs;
    fittedms2.medians = canonicalms2.medians .* ...
        repmat([xscale yscale], length(fittedms2.evtIDs), 1) + ...
        repmat([mydata.xbar mydata.stem], length(fittedms2.evtIDs), 1);
    if ~isnan(mydata.houtline)
        delete(mydata.houtline);
    end
    mydata.houtline = dg_plotMazespec(hA, fittedms2, 'mazespec2');
    mydata.ms2 = fittedms2;
    for k = 1:length(mydata.houtline)
        set(mydata.houtline(k), ...
            'ButtonDownFcn', @lfp_plot_eye_rodent_line);
    end
else
    disp('Please set all of stem, xbar, LGoal, Rgoal, start cursors first.');
end
set(hA, 'UserData', mydata);
end


function lfp_plot_eye_rodent_realign(hObject,evnt)
lfp_declareGlobals;
hF = get(hObject, 'Parent');
hA = findobj(hF, 'Type', 'axes');
hA = hA(end);   % legend lives in a new axes, we want the original
mydata = get(hA, 'UserData');
if isnan(mydata.houtline)
    disp('Please plot outline first.');
    return
end
hCheckbox = findobj(hF, 'Tag', 'lfp_plot_eye_rodent_origfile_checkbox');
useorigfile = get(hCheckbox, 'Value');
if useorigfile
    vtfile = '';
    vtfilenames = {'vt1.nvt' 'vt1.dat'};
    for k = 1:length(vtfilenames)
        if exist(fullfile(lfp_DataDir, vtfilenames{k}), 'file')
            vtfile = vtfilenames{k};
            break
        end
    end
    if isempty(vtfile)
        msg = sprintf('Session %s contains no video tracker file.', ...
            lfp_DataDir);
        lfp_log(msg);
        warning('lfp_realignVideo:novt', '%s', msg);
        return
    end
    disp('Reading original video file...');
    [VT_timestamps, VT_X, VT_Y] = Nlx2MatVT_v3( ...
        fullfile(lfp_DataDir, vtfile), [1 1 1 0 0 0], 0, 1,[] );
    VT_timestamps = VT_timestamps * 1e-6;
else
    framedur = 0.035151999999925;
    numframes = ( lfp_TimeStamps(end) + ...
        (lfp_SamplesPerFrame - 1) * lfp_SamplePeriod - lfp_TimeStamps(1) ...
        ) / framedur;
    VT_timestamps = lfp_TimeStamps(1) + (0:numframes-1) * framedur;
    frameidx = lfp_time2index(VT_timestamps);
    % eliminate repetitions caused by recording gaps:
    repeatedframe = diff(frameidx)==0;
    frameidx(repeatedframe) = [];
    VT_timestamps(repeatedframe) = [];
    VT_X = lfp_Samples{mydata.filenums(1)}(frameidx);
    VT_Y = lfp_Samples{mydata.filenums(2)}(frameidx);
end
myevtids = {[13]  mydata.midT  [14]  [15]  [16]  [17]  [18]};
mytestcols = [1      1       1     2     2     2     2];
gotevtid = false(size(myevtids));
for k = 1:length(myevtids)
    gotevtid(k) = any(cellfun(@isequal, ...
        mydata.ms2.evtIDs, repmat(myevtids(k), size(mydata.ms2.evtIDs)) ));
end
hSlider = findobj(hF, 'Tag', 'lfp_plot_eye_rodent_tol_slider');
frametol = get(hSlider, 'Value');
disp('Computing realignment...');
[new_VT_timestamps, gutreportstr, scorecard] = ...
    lfp_realignVideoGuts(myevtids(gotevtid), mytestcols(gotevtid), ...
    VT_timestamps, ...
    25, mydata.ms2, mydata.filenums, frametol);
disp(gutreportstr);
for k = 1:length(scorecard)
    disp(sprintf('trial %2d: %s', k, scorecard{k}));
end
[FileName, PathName] = uiputfile( ...
    '*.nvt', 'Saving realigned video tracker data', ...
    fullfile(lfp_DataDir, 'vt_realigned.nvt') );
if ~isequal(FileName, 0)
    disp('Saving realigned video tracker data...');
    new_VT_timestamps = round(new_VT_timestamps * 1e6);
    if ~ispc
        error('lfp_plot_eye_rodent:notPC', ...
            'Neuralynx format files can only be created on Windows machines.');
    end
    Mat2NlxVT_411(fullfile(PathName, FileName), 0, 1, 1, ...
        length(VT_X), [1 1 1 0 0 0 0], new_VT_timestamps, VT_X, VT_Y );
end
end


function lfp_plot_eye_rodent_origfile(hObject,evnt)
% Do nothing!
end


function lfp_plot_eye_rodent_tol(hObject,evnt)
hF = get(hObject, 'Parent');
value = get(hObject, 'Value');
step = get(hObject, 'SliderStep');
minval = get(hObject, 'Min');
maxval = get(hObject, 'Max');
value = minval + step(1) * (maxval-minval) * ...
    round((value-minval)/(step(1)*(maxval-minval)));
set(hObject, 'Value', value);
hText = findobj(hF, 'Tag', 'lfp_plot_eye_rodent_tol_text');
set(hText, 'String', sprintf('Tol = %03.1f', value));
end


function lfp_plot_eye_rodent_savemazespec(hObject,evnt)
lfp_declareGlobals;
hF = get(hObject, 'Parent');
hA = findobj(hF, 'Type', 'axes');
hA = hA(end);   % legend lives in a new axes, we want the original
mydata = get(hA, 'UserData');
if isnan(mydata.houtline)
    disp('Please plot outline first.');
    return
end
[FileName, PathName] = uiputfile( ...
    '*.mat', 'Saving T maze specifications', ...
    fullfile(lfp_DataDir, 'mazespec2.mat') );
mazespec2 = mydata.ms2;
if ~isequal(FileName, 0)
    save(fullfile(PathName, FileName), 'mazespec2');
end
end
