function lfp_cleanEyeEvents2(minsaccdur, minfixdur, minsaccint, minfixint)
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

% Modified by DG 8-Apr-2008 from a copy of lfp_cleanEyeEvents; saccade
% fusion has been replaced by incorporation into the next fixation.
% Specifically, for any series of saccades that violate minsaccint, the
% first in the series will be kept as is; if the last is separated from the
% following fixation by less than minfixint, then the second through last
% are ALL fused into the fixation.

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
    shortint = find(intervals < minsaccint);
    % But don't do it if any other eye events intervene!  (Non-eye events
    % are fine, because they could be absolutely anything.)
    if any(saccstarts(shortint+1) - saccends(shortint) > 1)
        % Yes, there some intervening events; must examine
        dontremove = [];
        for rmidx = 1:length(shortint)
            eventrange = ...
                saccends(shortint(rmidx)) + 1 : saccstarts(shortint(rmidx)+1) - 1;
            if any(ismember(lfp_Events(eventrange, 2), EyeEvents))
                dontremove(end+1) = rmidx;
            end
        end
        msg = sprintf(...
            'Of %d short inter-saccade intervals, %d were interrupted', ...
            length(shortint), length(dontremove) );
        lfp_log(msg);
        disp(msg);
        if length(dontremove)
            warning('lfp_cleanEyeEvents2:bogus', ...
                'Non-interruption assumption below may be false!' );
        end
        shortint(dontremove) = [];
    end
    % <shortint> points into the saccade lists where saccstarts(shortint+1)
    % is within minsaccint of saccends(shortint), so the new (9-Apr-08)  
    % requirement is to find runs of consecutive numbers in <shortint>, and
    % for each single (nonconsecutive) entry or run of consecutive entries
    % in <shortint>, examine the interval following the last saccend in the
    % run, and if it's less than minfixint, change the second saccstart in
    % the run to a fixstart, remove all the remaining saccade events in the
    % run, and remove the original fixstart.
    events2delete = [];
    numfused = 0;
    ptr = 1;
    while ptr < length(shortint)
        runstart = ptr; % first element in the run; points into shortint
        while ptr < length(shortint) && ...
                (shortint(ptr+1) == shortint(ptr) + 1)
            ptr = ptr + 1;
        end
        % ptr now points to the last shortint in the run:
        runend = ptr;
        % At this point, I rely on the assumption that if the interval to
        % the next fix start is less than minfixint, then there will be no
        % intervening eye events (see warning('lfp_cleanEyeEvents2:bogus'
        % ... )).  Therefore, we can stop searching when we find ANY eye
        % event, and if the interval is less than minfixint AND the event
        % is a fix start, we have a hit.
        evtix = saccends(shortint(runend) + 1) + 1;
        saccendTS = lfp_Events(evtix,1);
        while evtix <= size(lfp_Events,1) && ...
                ~ismember(lfp_Events(evtix,2), EyeEvents)
            evtix = evtix + 1;
        end
        if evtix <= size(lfp_Events,1)
            % evtix now points to next eye event in lfp_Events
            if (lfp_Events(evtix,1) - saccendTS < minfixint) ...
                    && (lfp_Events(evtix,2) == FixStart)
                % It's a fix start within minfixint of previous sacc end
                newfixstart = saccstarts(shortint(runstart)+1);
                lfp_Events(newfixstart, 2) ...
                    = FixStart;
                msg = sprintf('Moved fix %.6f to %.6f', ...
                    lfp_Events(evtix,1), ...
                    lfp_Events(newfixstart, 1) );
                del = [ saccstarts( ...
                    shortint(runstart) + 2 : shortint(runend) + 1 )
                    saccends(shortint(runstart) + 1 : shortint(runend) + 1)
                    evtix ];
                events2delete = [ events2delete
                    del
                    ];
                msg = sprintf('%s; deleted %s', ...
                    msg, mat2str(lfp_Events(del,1)) );
                lfp_log(msg);
                numfused = numfused + 1;
            end
        end
        % resume with next shortint after run
        ptr = runend + 1;
    end
    lfp_Events(events2delete,:) = [];
end
msg = sprintf('Fused %d saccade-fixation clusters', numfused);
disp(msg);
lfp_log(msg);

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
        msg = sprintf(...
            'Of %d short inter-fixation intervals, %d were interrupted', ...
            length(remove), length(dontremove) );
        disp(msg);
        lfp_log(msg);
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
