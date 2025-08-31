function lfp_replaceEyeFixIDs(fixations, varargin)
%lfp_replaceEyeFixIDs(fixations)
%lfp_replaceEyeFixIDs(fixations, 'insertnew')
%   <fixations> is a fixations table as returned by lfp_EyeTabulation.
%   Replaces all eye fixation start event IDs in lfp_Events with an ID code
%   whose decimal representation is of the form PTMMIIINNN, where:
%     P is the index number (i.e. the row number) of the task period as
%       defined by lfp_getTaskPeriods.
%     T is the temporal ordinal ID of the fixation target (since this is
%       a single digit, this code cannot be used for tasks that have more
%       than 9 temporally ordered sets of targets - Joey's has 3, Theresa's
%       has 1; T=0 means fixation was not on a relevant target);
%     MM (2 digits) is the ordinal number of the fixation of this spatial  
%       target in the current task period, or 0 for non-target fixations;
%     III (3 digits) is the spatial ID of fixation target (000 =
%       non-target);
%     NNN (3 digits) is the ordinal number of the fixation in the current 
%       task period (e.g., 1st, 2nd, ..., nth fixation in the current
%       trial;  starting at 001, so 000 is not used).
%   The leading non-zero P constrains all eye fixation start event IDs to
%   have values in the range > 1000000000.  Note that these extremely large
%   values make it impractical to extend event-related arrays such as
%   lfp_EventNames and lfp_EventColors to handle these ID codes.
%
%   It is assumed that lfp_SelectedTrials has the same value it did when
%   lfp_EyeTabulation was run; if its length is different, then an error is
%   signalled, but it could have the same length and still be different.

%
%   Requires that FixStart be assigned by lfp_getEvtIDs.
%
%   The specifics of the task are encapsulated in the functions
%   lfp_getTaskParams, lfp_getTaskPeriods, and lfp_relevantTargets, which
%   each user should modify according to their own needs.
%
% OPTIONS
%   'countT' - causes MM to count the number of fixations of the temporal
%       ordinal ID <T> rather than the spatial ID <III>.
%   'insertnew' - inserts new events with the new codes rather than
%       replacing the IDs in the existing events with the new codes.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;
allfixindices = find(lfp_Events(:,2) == FixStart);
trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if length(trials) ~= length(fixations)
    error('lfp_replaceEyeFixIDs:badSelect', ...
        'Trial selection state has changed since generating fixations table' );
end

countTflag = false;
insertnewflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'countT'
            countTflag = true;
        case 'insertnew'
            insertnewflag = true;
            newevents = zeros(length(allfixindices), 2);
            neweventidx = 0;
        otherwise
            error('lfp_replaceEyeFixIDs:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

taskperiods = lfp_getTaskPeriods;
for trial = trials
    fixtable = fixations{trial};
    if isempty(fixtable)
        warning('lfp_replaceEyeFixIDs:emptytable', ...
            'The fixations table for trial %d is empty; skipping.', ...
            trial );
        continue
    end
    startevtidx = lfp_TrialIndex(trial,1);
    endevtidx = lfp_TrialIndex(trial,2);
    % fixindices is an index into lfp_Events for this trial's fixations
    fixindices = allfixindices( ...
        allfixindices >= startevtidx & allfixindices <= endevtidx );
    % taskperiodtimes contains start (col 1) and end (col 2) timestamps in
    % this trial for each row of taskperiods
    taskperiodtimes = zeros(size(taskperiods,1),2);
    for pd = 1:size(taskperiods,1)
        startidx = find(ismember( ...
            lfp_Events(startevtidx:endevtidx,2), taskperiods{pd,2} )) ...
            + startevtidx - 1 ;
        endidx = find(ismember( ...
            lfp_Events(startevtidx:endevtidx,2), taskperiods{pd,3} )) ...
            + startevtidx - 1 ;
        if length(startidx) > 1 || length(endidx) > 1
            warning('lfp_replaceEyeFixIDs:multibound', ...
                'Multiple boundary events for period %s in trial %d', ...
                taskperiods{pd,1}, trial );
            startidx = startidx(1);
            endidx = endidx(1);
        end
        taskperiodtimes(pd, :) = [ ...
            lfp_Events(startidx,1) lfp_Events(endidx,1) ];
    end
    prevtaskperiodidx = 0;
    % fixidx is an index into lfp_Events
    for fixidx = fixindices'
        fixTS = lfp_Events(fixidx,1);
        % fixrow is a row index into <fixtable>
        fixrowidx = find(fixtable(:,3) > fixTS - 5e-7 ...
            & fixtable(:,3) < fixTS + 5e-7 );
        if isempty(fixrowidx)
            warning('lfp_replaceEyeFixIDs:nofixrow', ...
                'Fixation at %d is not in fixation table for trial %d; skipping.', ...
                fixTS, trial );
            continue
        end
        if size(fixrowidx,1) > 1
            warning('lfp_replaceEyeFixIDs:multiTS', ...
                'Multiple timestamp matches to fixation at %d; using first.', ...
                fixTS );
            fixrowidx(2:end,:) = [];
        end
        taskperiodidx = find( fixTS >= taskperiodtimes(:,1) ...
            & fixTS < taskperiodtimes(:,2) );
        if isempty(taskperiodidx)
            warning('lfp_replaceEyeFixIDs:noperiod', ...
                'Fixation at %d is not in any task period; skipping.', fixTS );
            continue
        end
        if taskperiodidx ~= prevtaskperiodidx
            % Re-zero counters at start of each task period
            mm = zeros(999,1);    % vector of MM values for this task period
            fixnum = 0; % within this task period
            prevtaskperiodidx = taskperiodidx;
        end
        % targID contains III value
        targID = lfp_xy2targID( ...
            fixtable(fixrowidx,1), fixtable(fixrowidx,2), ...
            lfp_getTargets(trial, taskperiods{taskperiodidx}) ); 
        MM = 0;
        [targetIDs] = ...
            lfp_relevantTargets(trial, taskperiods{taskperiodidx});
        T = find(ismember(targetIDs, targID));
        if isempty(T)
            T = 0;
        end
        if targID > 0
            if countTflag
                idx = T;
            else
                idx = targID;
            end
            if idx > 0
                mm(idx) = mm(idx) + 1;
                MM = mm(idx);
            end
        end
        fixnum = fixnum + 1;
        newEvtID = fixnum + 1000*targID + 1e6*MM ...
            + 1e8*T + taskperiodidx*1e9;
        if insertnewflag
            neweventidx = neweventidx + 1;
            newevents(neweventidx,:) = [ fixTS newEvtID ];
        else
            lfp_Events(fixidx,2) = newEvtID;
        end
    end
end

if insertnewflag
    lfp_Events = [lfp_Events
        newevents ];
    lfp_Events = sortrows(lfp_Events);
    disp('Recalculating lfp_TrialIndex...');
    lfp_TrialIndex = lfp_createTrialIndex(lfp_Events);
    if ~isempty(lfp_SamplePeriod)
        [startSampleIndex, endSampleIndex] = ...
            lfp_addTIcols3n4(lfp_TrialIndex(:,1), lfp_TrialIndex(:,2), lfp_Events);
        lfp_TrialIndex = [lfp_TrialIndex startSampleIndex endSampleIndex];
    end
end

