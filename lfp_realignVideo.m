function scorecard = lfp_realignVideo(sessiondir, varargin)
%scorecard = lfp_realignVideo(sessiondir)
% Meant to work with all types of mouse and rat sessions. Takes raw
% Neuralynx video tracker file or two CSC channels as input, uses
% 'mazespec2.mat' file in <sessiondir> as video calibration. Produces a
% video file with altered timestamps, which by default is a Yasuo-format
% (VT output) file saved to the animal directory. Note that the trial
% numbers in the output file may not agree with those in some other
% existing Yasuo file, but they will agree with the trial numbers assigned
% by lfp_lib.  Displays and logs a report after all processing is done, and
% returns a trial-by-trial summary scorecard, each of whose elements is one
% of
%   'good' (no adjustment required)
%   'bad' (no adjustment possible)
%   'well-adjusted' (an easy-to-determine time offset)
%   'adjusted with difficulty'
% OPTIONS
%   'convertcalib' - look up the saved 'mazespec' (note that this is
%       different from 'mazespec2') and convert it to a mazespec2.  Note
%       that the 'mazespec' video image model is intrinsically inaccurate,
%       as it fails to address the issues of keystoning and foreshortening.
%   'filenums', filenums - uses the traces recorded in CSC channels
%       <filenums> as video tracker data instead of raw Neuralynx file.
%       <filenums> must be a 2-element array.  Default is to read from raw
%       Neuralynx video tracker file.
%   'frametol', frametol - number of frames of timestamp slippage that are
%       tolerated as unavoidable error.  Default is 1.
%   'maxpixerr', maxpixerr - <maxpixerr> is a threshold used to determine
%       whether each event is consistently aligned between video and events
%       file.  Default value is 20 pixels.
%   'mazespec2', ms2 - uses <ms2> as the mazespec2 for the session instead
%       of reading the one in the session directory or creating one.
%   'offsetplot' - creates a figure containing plots of offset vs. trial
%       number and offset vs. time.  For the offset vs. trial plot, the
%       offset at the beginning of the trial is used.
%   'output', outfilename - writes output to <outfilename>.  If
%       <outfilename> has no directory component, then the file is written
%       to <sessiondir>.  If <outfilename> has a '.nvt' extension, then a
%       Neuralynx-format file is written.
%   'selfcalib' - create a mazespec2 using the actual session data; but
%       note that this will fail to correct a systematic offset that
%       affects the entire session, and will most likely fail disastrously
%       if many trials have significant offsets.
%   'trialfrac', trialfrac - the fraction of trials that must be "good" to
%       do the realignment.  Default = 0.5.  Trials are considered "bad" if
%       they have breaks (NaNs) anywhere in either video trace during the
%       common trial time aligned on lfp_NominalTrialStart or [15 16 17
%       18].  (This definition is a computational convenience which
%       substitutes for finding trials that contain any breaks between
%       lfp_NominalTrialStart and lfp_NominalTrialEnd.)
% NOTES
% Assumes that sessiondir is indeed a session directory, contained within
% an animal directory.

%$Rev: 333 $
%$Date: 2014-10-16 18:56:27 -0400 (Thu, 16 Oct 2014) $
%$Author: dgibson $

scorecard = [];
maxpixerr = 20;    % maximum position error permitted, in pixels
frametol = 1;
trialfrac = 0.5;    % fraction of trials that must be "good" to continue

argnum = 1;
convertcalibflag = false;
filenums = [];
ms2 = [];
outfilename = '';
offsetplotflag = false;
selfcalibflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'convertcalib'
            convertcalibflag = true;
        case 'filenums'
            argnum = argnum + 1;
            filenums = varargin{argnum};
        case 'frametol'
            argnum = argnum + 1;
            frametol = varargin{argnum};
        case 'maxpixerr'
            argnum = argnum + 1;
            maxpixerr = varargin{argnum};
        case 'mazespec2'
            argnum = argnum + 1;
            ms2 = varargin{argnum};
        case 'offsetplot'
            offsetplotflag = true;
        case 'output'
            argnum = argnum + 1;
            outfilename = varargin{argnum};
        case 'selfcalib'
            selfcalibflag = true;
        case 'trialfrac'
            argnum = argnum + 1;
            trialfrac = varargin{argnum};
        otherwise
            error('lfp_realignVideo:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if (selfcalibflag + ~isempty(ms2) + convertcalibflag) > 1
    error('lfp_realignVideo:options', ...
        'You may only specify one of ''mazespec2'', ''convertcalib'', or ''selfcalib''');
end
if ~isempty(filenums) && numel(filenums) ~= 2
    error('lfp_realignVideo:filenums', ...
        'There must be exactly two elements in <filenums>');
end

% Find events file:
evfile = '';
evtfilenames = {'events_90.evtsav' 'events.evtsav' 'events.nev' 'events.dat'};
for k = 1:length(evtfilenames)
    if exist(fullfile(sessiondir, evtfilenames{k}), 'file')
        evfile = evtfilenames{k};
        break
    end
end
if isempty(evfile)
    msg = sprintf('Session %s contains no events file.', sessiondir);
    lfp_log(msg);
    error('lfp_realignVideo:noevts', '%s', msg);
end
if isempty(filenums)
    % Find video tracker file:
    vtfile = '';
    vtfilenames = {'vt1.nvt' 'vt1.dat'};
    for k = 1:length(vtfilenames)
        if exist(fullfile(sessiondir, vtfilenames{k}), 'file')
            vtfile = vtfilenames{k};
            break
        end
    end
    if isempty(vtfile)
        msg = sprintf('Session %s contains no video tracker file.', ...
            sessiondir);
        lfp_log(msg);
        error('lfp_realignVideo:novt', '%s', msg);
    end
end

% Figure out what kind of data we have:
[TS, TTL] = Nlx2MatEV(fullfile(sessiondir, evfile), [1, 0, 1, 0, 0], 0, 1);
strobed = TTL < 0;
TTL = TTL(strobed) + 2^15;
eventIDs = unique(TTL);
eventcounts = NaN(size(eventIDs));
for k = 1:length(eventIDs)
    eventcounts(k) = sum(TTL == eventIDs(k));
end
% Assume that every session is either rat or mouse; thus ~israt implies
% that it is a mouse session:
if ismember(90, eventIDs) && ismember(91, eventIDs) ...
        && ismember(92, eventIDs)
    israt = eventcounts(eventIDs == 90) > ...
        (eventcounts(eventIDs == 91) +  eventcounts(eventIDs == 92));
else
    israt = ismember(90, eventIDs);
end
if israt
    numtrials = eventcounts(eventIDs == 90);
else
    numtrials = eventcounts(eventIDs == 91) +  eventcounts(eventIDs == 92);
end
% isauditory and istactile are not mutually exclusive, so a completely
% intact auditory-tactile session will have tactile cue events in half the
% trials and auditory cue events in the other half:
isauditory = sum(ismember(TTL, [31 38])) > numtrials/3;
istactile = israt && (sum(ismember(TTL, [21 22])) > numtrials/3) || ...
    ~israt && sum(TTL==20) > numtrials/3;
gateidx = find(TTL==11);
if isempty(gateidx)
    msg= sprintf('Session %s contains no gate events.', sessiondir);
    lfp_log(msg);
    error('lfp_realignVideo:nogate', '%s', msg);
end
iscuefirst = any(ismember(TTL(1:gateidx(1)), [31 38])) && ...
    ~any(ismember(TTL(gateidx(end):end), [31 38]));
isstandard = ~any(ismember(TTL(1:gateidx(1)), [31 38])) && ...
    any(ismember(TTL(gateidx(end):end), [21 22 31 38]));
reportstr = sprintf('Session %s is ', sessiondir);

% Choose appropriate midT events and setup, and start report:
midT = [];
if israt
    if iscuefirst
        midT = 23;
        setupname = 'jianbin';
        reportstr = [reportstr 'cue-first '];
    else
        if isauditory
            midT = [midT 31 38];
        end
        if istactile
            if israt
                midT = [midT 21 22];
            end
        end
        if isstandard
            reportstr = [reportstr 'standard '];
        else
            reportstr = [reportstr 'unknown '];
        end
        setupname = 'katy'; % has tactile and auditory cue values
    end
    reportstr = [reportstr 'rat'];
else
    setupname = 'mouse';
    reportstr = [reportstr 'mouse'];
    midT = [midT 20];
end
if isauditory && istactile
    reportstr = sprintf('%s tactile-auditory.\n', reportstr);
    if isempty(midT)
        midT = [31 38 21 22];
    end
end
if isauditory && ~istactile
    reportstr = sprintf('%s auditory.\n', reportstr);
end
if ~isauditory && istactile
    reportstr = sprintf('%s tactile.\n', reportstr);
end
if ~isauditory && ~istactile
    reportstr = sprintf('%s unknown modality.\n', reportstr);
end
lfp_changeSetup(setupname);

lfp_declareGlobals;

%load session
lfp_read2('preset', sessiondir, {evfile ''});
if isempty(filenums)
    lfp_add('preset', sessiondir, {vtfile}, ...
        'Neuralynx Video (*.NVT, *.DAT)', 0);
    XYchan = lfp_ActiveFilenums(end-1:end);
else
    XYchan = filenums;
end
lfp_XLimAll = [];
lfp_AlignmentRef = lfp_NominalTrialStart;
lfp_disp([],XYchan,[],'avg','marknanbad','noplot');
lfp_AlignmentRef = [15 16 17 18];
lfp_disp([],XYchan,[],'avg','marknanbad','noplot');
if length(lfp_BadTrials) > length(lfp_SelectedTrials)*(1-trialfrac)
    msg = sprintf( ...
        'Session %s contains badly corrupted video data, skipping.', ...
        sessiondir);
    lfp_log(msg);
    error('lfp_realignVideo:badvideo', '%s', msg);
end
lfp_BadTrials = [];
evtIDs = {13 6 midT 9 14 15 16 7 8 17 18};
if selfcalibflag
    % Get coordinates from this very session:
    ms2 = lfp_measureTMaze2(XYchan(1), XYchan(2), ...
        'evtIDs', evtIDs, 'plot', 'outline');
    msg = sprintf('Using self-calibration');
    lfp_log([sessiondir ': ' msg]);
    warning('lfp_realignVideo:nocalib', '%s', msg);
end
if convertcalibflag
    [mazespec, sessiondatestr] = lfp_findMazespec;
    ms2 = dg_convertMazespec(mazespec, evtIDs);
end
if isempty(ms2)
    ms2file = fullfile(sessiondir, 'mazespec2.mat');
    if exist(ms2file, 'file')
        S = load('-mat', ms2file);
        myfields = fieldnames(S);
        if ismember('mazespec2', myfields)
            ms2 = S.mazespec2;
        else
            error('lfp_realignVideo:mazespec2', ...
                '%s does not contain a mazespec2', ms2file);
        end
    else
        error('lfp_realignVideo:noms2', ...
            'There is no mazespec2 file.');
    end
end
if isempty(filenums)
    [VT_timestamps, VT_X, VT_Y] = Nlx2MatVT_v3( ...
        fullfile(sessiondir, vtfile), [1 1 1 0 0 0], 0, 1,[] );
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

% For each trial, find spatial coordinates of all events that exist among
% {13 midT 14 [15 16] [17 18]}, and if any are "outliers", correct slippage.
myevtids =  {13 midT  14  15 16 17 18};
mytestcols = [1 1     1   2  2  2  2 ];
gotevtid = false(size(myevtids));
for k = 1:length(myevtids)
    gotevtid(k) = any(cellfun(@isequal, ...
        ms2.evtIDs, repmat(myevtids(k), size(ms2.evtIDs)) ));
end
[new_VT_timestamps, gutreportstr, scorecard] = ...
    lfp_realignVideoGuts(myevtids(gotevtid), mytestcols(gotevtid), ...
    VT_timestamps, maxpixerr, ms2, XYchan, frametol);
reportstr = sprintf('%s\n%s', reportstr, gutreportstr);
lfp_log(reportstr);
disp(reportstr);

if offsetplotflag
    offsets = new_VT_timestamps - VT_timestamps;
    trialstartTS = lfp_Events(lfp_TrialIndex(:,1), 1);
    trialoffsets = NaN(size(lfp_SelectedTrials));
    for k = 1:length(trialoffsets)
        [d, videoframenum] = min(abs(new_VT_timestamps - trialstartTS(k)));
        trialoffsets(k) = offsets(videoframenum);
    end
    hOffsetFig = figure;
    hA = subplot(2, 1, 1);
    plot(new_VT_timestamps, offsets);
    title('Video Timestamp Offset vs. Time');
    xlabel('Time, s');
    ylabel('offset, s');
    figure(hOffsetFig);
    hA = subplot(2, 1, 2);
    plot(1:length(trialoffsets), trialoffsets);
    title('Video Timestamp Offset vs. Trial Number');
    xlabel('Trial');
    ylabel('offset, s');
end

if isempty(outfilename)
    [animaldir, sessionstr] = fileparts(sessiondir);
    outfilename = ['v2' sessionstr '.dat'];
    if exist(fullfile(animaldir, outfilename), 'file')
        n = 1;
        while exist(fullfile(animaldir, [outfilename '.bak' num2str(n)]), ...
                'file')
            n = n + 1;
        end
        warning('lfp_saveVTfile:mvfile', ...
            'Moving file %s to %s', fullfile(animaldir, outfilename), ...
            fullfile(animaldir, [outfilename '.bak' num2str(n)]) );
        movefile(fullfile(animaldir, outfilename), ...
            fullfile(animaldir, [outfilename '.bak' num2str(n)]));
    end
    outfilename = fullfile(animaldir, outfilename);
end
[outdir, outname, outext] = fileparts(outfilename);
if isequal(lower(outext), '.nvt')
    if ~ispc
        error('lfp_realignVideo:notPC', ...
            'Neuralynx format files can only be created on Windows machines.');
    end
    Mat2NlxVT_411(outfilename, 0, 1, 1, ...
        length(VT_X), [1 1 1 0 0 0 0], new_VT_timestamps, VT_X, VT_Y );
else
    lfp_saveVTfile( outfilename, new_VT_timestamps, VT_X, VT_Y );
end

end

