function [windows, winsamples] = lfp_collectDataWindows(channel, trials, ...
    evtbounds, offsets, window, noclip, nodataOKflag, varargin)
%[windows, winsamples] = lfp_collectDataWindows(channel, trials, ...
%    evtbounds, offsets, window)
%INPUTS
%  <trials>: Each enabled trial (i.e. lfp_SelectedTrials(trial) is 'true'
%  and <trial> is not a member of lfp_BadTrials) in the integer vector
%  <trials> is used to compute results. <trials> can be a row vector or a
%  column vector, or it can be empty, in which case all selected trials are
%  used (equivalent to specifying <trials> as
%  lfp_enabledTrials(find(lfp_SelectedTrials))). If <trials> is a string,
%  then it is interpreted as containing Unique Trial IDs, i.e. the
%  combination of a session name and a trial number as it was originally
%  numbered in that session.  The session name component is optional (see
%  lfp_parseTrialStr and lfp_getTrialNum for further details). 
%
%  <evtbounds> is a 1x2 cell array.  The first cell should contain
%  a list of event IDs, and the event in the list that occurs earliest in
%  the trial is used to specify the start time. The second cell works similarly to
%  specify end time.  
%
%  <offsets> is a 1x2 numeric array.  Each element represents an offset in
%  seconds that is added to the time of the corresponding event in
%  <evtbounds> to compute the final start and end times that delimit the
%  data to use to compute results.  The end time is truncated so that the
%  interval from start time to end time is an integral multiple of
%  <window>.
%
%  <window> specifies the width of the moving time window used to select
%  data.  The step size for the moving window is exactly one window, so
%  successive windows are disjoint and therefore strictly independent.  If
%  the value is 0, then the entire time period specified by <evtbounds> and
%  <offsets> is used.
%
%  <noclip> - if 'auto', then the maximum and minimum
%  values found in <channel> are treated as clipping values.
%  Otherwise, it should be a 2-element array containing the lower
%  and upper clipping values in any order.  Any windows containing
%  more than one point with either of those values will be excluded.
%
%  <nodataOKflag> - if 'true', issues a warning and keeps going instead of
%  raising an error when a trial has some fatal flaw in the data. 
%OUTPUTS
% <windows> is a two-column matrix with starting sample indices in column 1
%   and ending sample indices in col 2.  Each row represents one window.
%   Although the windows are collected in a trial-oriented manner, trial
%   boundaries are not indicated in the result.
% <winsamples> is the number of samples in one data window.
%OPTIONS
% 'margintime', margintime - adds <margintime> seconds of data before and
%   after each entry in <windows>.  This does not alter the step size, so
%   the end result is that each window overlaps with its neighbor by
%   <margintime> seconds.

%$Rev: 389 $
%$Date: 2017-06-12 15:07:03 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

margintime = 0;
argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'margintime'
            argnum = argnum + 1;
            margintime = varargin{argnum};
        otherwise
            error('funcname:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end


global lfp_ActiveFilenums lfp_Events lfp_Samples lfp_SamplePeriod ...
    lfp_SelectedTrials lfp_TrialIndex lfp_TrialRec

if isempty(trials)
    trials = 1:length(lfp_SelectedTrials);
end

if ischar(trials)
    trials = lfp_parseTrialStr(trials, session);
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

if ~ismember(channel, lfp_ActiveFilenums)
    error('lfp_collectDataWindows:noSuchFilenum', ...
        'You requested non-existent file number: %d', channel);
end
if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_collectDataWindows:badTrials1', ...
            '<trials> must be an integer vector.');
    end
    trials = trials';
end
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_collectDataWindows:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
        num2str(trials( ...
        trials > size(lfp_TrialIndex,1) | trials < 1 )) ]);
end
if ~isa(trials, 'double') || ~all(fix(trials) == trials)
    error('lfp_collectDataWindows:badTrials2', ...
        '<trials> must be an integer vector.');
end

if ~iscell(evtbounds)
    error('lfp_collectDataWindows:badevtbounds', ...
        '<evtbounds> must be a 1 X 2 cell array.');
elseif ~isequal(size(evtbounds), [1 2])
    error('lfp_collectDataWindows:badevtbounds2', ...
        '<evtbounds> must be a 1 X 2 cell array.');
end

if ~isnumeric(offsets)
    error('lfp_collectDataWindows:offsets', ...
        '<offsets> must be a 1 X 2 numeric array.');
elseif ~isequal(size(evtbounds), [1 2])
    error('lfp_collectDataWindows:offsets2', ...
        '<offsets> must be a 1 X 2 numeric array.');
end

if ~isnumeric(window)
    error('lfp_collectDataWindows:window', ...
        '<window> must be a scalar number.');
elseif ~isequal(size(window), [1 1])
    error('lfp_collectDataWindows:window2', ...
        '<window> must be a scalar number.');
end

if ~isnumeric(margintime)
    error('lfp_collectDataWindows:margintime', ...
        '<margintime> must be a scalar number.');
elseif ~isequal(size(margintime), [1 1])
    error('lfp_collectDataWindows:margintime2', ...
        '<margintime> must be a scalar number.');
end

if isequal(noclip, 'auto')
    noclip = [min(lfp_Samples{channel}(:)) max(lfp_Samples{channel}(:))];
elseif ~isempty(noclip) && ( .../
        ~(isequal(size(noclip), [1 2]) || isequal(size(noclip), [2 1])) ...
        || ~isnumeric(noclip) )
    error('lfp_collectDataWindows:noclip', ...
        '<noclip> must be a 2-element numeric array.');
end

if window > 0
    winsamples = floor(window/lfp_SamplePeriod);
    % in the unlikely case that window is an integral number of samples, we do
    % not want to include the last sample:
    if winsamples == window/lfp_SamplePeriod
        winsamples = winsamples - 1;
    end
else
    winsamples = 0;
end
windows = cell(0);

for trialidx = 1:length(trials)
    trial = trials(trialidx);
    % find start time and endtime for data:
    startevtidx = lfp_TrialIndex(trial,1);
    endevtidx = lfp_TrialIndex(trial,2);
    startevtidx = find(...
        ismember(...
        lfp_Events(startevtidx:endevtidx, 2), evtbounds{1} )) ...
        + startevtidx - 1;
    if isempty(startevtidx)
        errid = 'lfp_collectDataWindows:evtbounds1';
        errmsg = sprintf(...
            'Trial %d has no ''evtbounds'' start event', trial );
        if nodataOKflag
            warning(errid, '%s', errmsg);
            continue
        else
            error(errid, '%s', errmsg);
        end
    else
        startevtidx = startevtidx(1);
    end
    endevtidx = find(...
        ismember(...
        lfp_Events(startevtidx:endevtidx, 2), evtbounds{2} )) ...
        + startevtidx - 1;
    if isempty(endevtidx)
        errid = 'lfp_collectDataWindows:evtbounds2';
        errmsg = sprintf(...
            'Trial %d has no ''evtbounds'' end event', trial );
        if nodataOKflag
            warning(errid, '%s', errmsg);
            continue
        else
            error(errid, '%s', errmsg);
        end
    else
        endevtidx = endevtidx(1);
    end

    % Verify that the time range of data to process is within the bounds of
    % the trial's recorded time segment.
    timeinterval(1) = lfp_Events(startevtidx,1) + offsets(1) - margintime;
    timeinterval(2) = lfp_Events(endevtidx,1) + offsets(2) + margintime;
    startsample = lfp_time2index(timeinterval(1));
    endsample = lfp_time2index(timeinterval(2));
    if startsample < lfp_TrialRec(trial,1)
        errid = 'lfp_collectDataWindows:startsample';
        errmsg = sprintf(...
            'Trial %d start time precedes recorded time', trial );
        if nodataOKflag
            warning(errid, '%s', errmsg);
            continue
        else
            error(errid, '%s', errmsg);
        end
    elseif endsample > lfp_TrialRec(trial,2)
        errid = 'lfp_collectDataWindows:endsample';
        errmsg = sprintf(...
            'Trial %d end time exceeds recorded time', trial );
        if nodataOKflag
            warning(errid, '%s', errmsg);
            continue
        else
            error(errid, '%s', errmsg);
        end
    end

    if winsamples > 0
        % Compute the start sample end sample of
        % each window to be analyzed in this trial. <winsamples> does not
        % include <margintime>, but <startsample> and <endsample> do.
        % Note that <windows> is
        % a cell array at this point, so there is one allocation here for
        % each trial, plus one more allocation by cell2mat after the total
        % size is known.
        marginpts = round(margintime/lfp_SamplePeriod);
        nwindows = floor( (endsample - startsample + 1 - 2 * marginpts) ...
            / winsamples );
        windows{trialidx,1} = zeros(nwindows,2);
        windows{trialidx,1}(:,1) = startsample + ...
            (0:nwindows-1)' * winsamples;
        windows{trialidx,1}(:,2) = windows{trialidx}(:,1) + winsamples - 1;
    else
        windows{trialidx,1} = [startsample endsample];
    end
end
windows = cell2mat(windows);
isbadwin = false(size(windows,1), 1);
if ~isempty(noclip)
    for k = 1:size(windows,1)
        if sum( lfp_Samples{channel}(windows(k,1):windows(k,2)) == ...
                noclip(1) | ...
                lfp_Samples{channel}(windows(k,1):windows(k,2)) == ...
                noclip(2) ) > 2
            isbadwin(k) = true;
        end
    end
    msg = sprintf( ...
        'lfp_collectDataWindows eliminated %d clipped windows, of %d', ...
        sum(isbadwin), length(isbadwin) );
    disp(msg);
    lfp_log(msg);
    windows(isbadwin, :) = [];
end

