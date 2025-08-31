function lfp_cleanEyeEvents(minsaccdur, minfixdur, minsaccint, minfixint)
%lfp_cleanEyeEvents(minsaccdur, minfixdur, minsaccint, minfixint)
%   Edits lfp_Events to eliminate non-paired starts and ends of saccades,
%   and likewise for fixations (the closer starts and ends are kept).
%   Edits lfp_Events to delete events that are too brief, and then to fuse
%   repeated events that are sufficiently close in time.  Any edit whose
%   minimum time parameter is set to zero is skipped.
%   minsaccdur:  min saccade duration in seconds
%   minfixdur:  min fixation duration in seconds
%   minsaccint:  min interval between successive saccades in seconds
%   minfixint:  min interval between successive fixations in seconds

% NOTE:  This code depends on these eye event constants being defined in
% the appropriate lfp_getEvtIDs_* script:  SaccStart, SaccEnd, FixStart,
% FixEnd, EyeEvents.

% OTHER NOTE:  If this function ever needs to be substantially modified
% again, this duplication should be eliminated:  the same formal operations
% are being performed on fixations as on saccades, and those operations
% should be called as functions twice rather than copied-and-pasted as they
% are now.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;

saccstarts = find(lfp_Events(:,2) == SaccStart);
saccends = find(lfp_Events(:,2) == SaccEnd);
[pairs, extrastart, extraend] = dg_zip(saccstarts, saccends);
disp(sprintf(...
    'Removing %d extra saccade starts, %d extra saccade ends', ...
    length(extrastart), length(extraend) ));
lfp_Events([saccstarts(extrastart); saccends(extraend)],:) = [];
fixstarts = find(lfp_Events(:,2) == FixStart);
fixends = find(lfp_Events(:,2) == FixEnd);
[pairs, extrastart, extraend] = dg_zip(fixstarts, fixends);
disp(sprintf(...
    'Removing %d extra fixation starts, %d extra fixation ends', ...
    length(extrastart), length(extraend) ));
lfp_Events([fixstarts(extrastart); fixends(extraend)],:) = [];

% Remove saccades that are too short to be real:
if minsaccdur > 0
    saccstarts = find(lfp_Events(:,2) == SaccStart);
    saccends = find(lfp_Events(:,2) == SaccEnd);
    durations = lfp_Events(saccends,1) - lfp_Events(saccstarts,1);
    remove = find(durations < minsaccdur);
    disp(sprintf(...
        'Removing %d short saccades', ...
        length(remove) ));
    lfp_Events([saccstarts(remove); saccends(remove)],:) = [];
end

% Remove fixations that are too short to be real:
if minfixdur > 0
    fixstarts = find(lfp_Events(:,2) == FixStart);
    fixends = find(lfp_Events(:,2) == FixEnd);
    durations = lfp_Events(fixends,1) - lfp_Events(fixstarts,1);
    remove = find(durations < minfixdur);
    disp(sprintf(...
        'Removing %d short fixations', ...
        length(remove) ));
    lfp_Events([fixstarts(remove); fixends(remove)],:) = [];
end

% Fuse saccades that are separated by too little time to be
% distinct:
if minsaccint > 0
    saccstarts = find(lfp_Events(:,2) == SaccStart);
    saccends = find(lfp_Events(:,2) == SaccEnd);
    intervals = lfp_Events(saccstarts(2:end),1) ...
        - lfp_Events(saccends(1:end-1),1);
    remove = find(intervals < minsaccint);
    % But don't do it if any other eye events intervene!  (Non-eye events
    % are fine, because they could be absolutely anything.)
    if any(saccstarts(remove+1) - saccends(remove) > 1)
        % Yes, there some intervening events; must examine
        dontremove = [];
        for rmidx = 1:length(remove)
            eventrange = ...
                saccends(remove(rmidx)) + 1 : saccstarts(remove(rmidx)+1) - 1;
            if any(ismember(lfp_Events(eventrange, 2), EyeEvents))
                dontremove(end+1) = rmidx;
            end
        end
        disp(sprintf(...
            'Of %d short inter-saccade intervals, %d were interrupted', ...
            length(remove), length(dontremove) ));
        remove(dontremove) = [];
    end
    disp(sprintf(...
        'Removing %d short inter-saccade intervals', ...
        length(remove) ));
    lfp_Events([saccends(remove); saccstarts(remove+1)],:) = [];
end

% Fuse fixations that are separated by too little time to be
% distinct:
if minfixint > 0
    fixstarts = find(lfp_Events(:,2) == FixStart);
    fixends = find(lfp_Events(:,2) == FixEnd);
    intervals = lfp_Events(fixstarts(2:end),1) ...
        - lfp_Events(fixends(1:end-1),1);
    remove = find(intervals < minfixint);
    % But don't do it if any other eye events intervene!
    if any(fixstarts(remove+1) - fixends(remove) > 1)
        dontremove = [];
        for rmidx = 1:length(remove)
            eventrange = ...
                fixends(remove(rmidx)) + 1 : fixstarts(remove(rmidx)+1) - 1;
            if any(ismember(lfp_Events(eventrange, 2), EyeEvents))
                dontremove(end+1) = rmidx;
            end
        end
        disp(sprintf(...
            'Of %d short inter-fixation intervals, %d were interrupted', ...
            length(remove), length(dontremove) ));
        remove(dontremove) = [];
    end
    disp(sprintf(...
        'Removing %d short inter-fixation intervals', ...
        length(remove) ));
    lfp_Events([fixends(remove) fixstarts(remove+1)],:) = [];
end

lfp_Events = sortrows(lfp_Events);
disp('Recalculating lfp_TrialIndex...');
lfp_TrialIndex = lfp_createTrialIndex(lfp_Events);
[startSampleIndex, endSampleIndex] = ...
    lfp_addTIcols3n4(lfp_TrialIndex(:,1), lfp_TrialIndex(:,2), lfp_Events);
lfp_TrialIndex = [lfp_TrialIndex startSampleIndex endSampleIndex];

lfp_log(sprintf('lfp_cleanEyeEvents %d %d %d %d', ...
    minsaccdur, minfixdur, minsaccint, minfixint ));
