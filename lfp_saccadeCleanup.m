function lfp_saccadeCleanup(x, y, velo, threshes, maxvsaccstart)

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIds;

saccEvtIDs = [SaccStart SaccEnd
    CplxSaccStart CplxSaccEnd
    BTSaccStart BTSaccEnd];
toremove = false(size(lfp_Events,1), 1);
fix2mergeevts = zeros(0, 4); % start, end to remove, start to remove, end
saccstartevts = find(ismember(lfp_Events(:,2), saccEvtIDs(:,1)));
saccendevts = find(ismember(lfp_Events(:,2), saccEvtIDs(:,2)));

if length(saccstartevts) ~= length(saccendevts)
    msg = sprintf('different numbers of saccade starts and ends');
    warning('lfp_saccadeCleanup:mismatch', ...
        '%s', ...
        msg );
    lfp_log(['lfp_saccadeCleanup: ' msg]);
end
minlen = min(length(saccstartevts), length(saccendevts));
if any(saccstartevts(1:minlen) > saccendevts(1:minlen)) ...
        || any(saccendevts(1:minlen-1) > saccstartevts(2:minlen))
    msg = 'There are disordered saccade events';
    lfp_log(['lfp_saccadeCleanup: ' msg]);
    warning('lfp_saccadeCleanup:disordered1', '%s', msg);
end

saccstarttimes = lfp_Events(saccstartevts, 1);
diffsst = diff(saccstarttimes);
shortdifftime = 0.120;
is_short = diffsst < shortdifftime;
msglist = {sprintf('saccade starts <.120s apart:\n')};
% i.	Implement new rules for whole trial for saccade starts <.120s apart
for saccnum = find(is_short)'
    sacc1start = lfp_time2index(lfp_Events(saccstartevts(saccnum), 1));
    evtidx = saccstartevts(saccnum) + 1;
    while evtidx <= size(lfp_Events, 1) && ~ismember( ...
            lfp_Events(evtidx, 2), saccEvtIDs(:,2) )
        evtidx = evtidx + 1;
    end
    if ~ismember(lfp_Events(evtidx, 2), saccEvtIDs(:,2))
        error('lfp_saccadeCleanup:nosaccend', ...
            'The first sacc has no end.' );
    end
    sacc1endevt = evtidx;
    sacc1end = lfp_time2index(lfp_Events(sacc1endevt, 1));
    sacc2start = lfp_time2index(lfp_Events(saccstartevts(saccnum+1), 1));
    % start searching for sacc2end at the later of sacc1end or sacc2start:
    evtidx = max(sacc1endevt, saccstartevts(saccnum+1)) + 1;
    while evtidx <= size(lfp_Events, 1) && ~ismember( ...
            lfp_Events(evtidx, 2), saccEvtIDs(:,2) )
        evtidx = evtidx + 1;
    end
    if ~ismember(lfp_Events(evtidx, 2), saccEvtIDs(:,2))
        error('lfp_saccadeCleanup:nosaccend', ...
            'The second sacc has no end.' );
    end
    sacc2endevt = evtidx;
    sacc2end = lfp_time2index(lfp_Events(sacc2endevt, 1));
    saccv1 = max(lfp_Samples{velo}(sacc1start:sacc1end));
    saccv2 = max(lfp_Samples{velo}(sacc2start:sacc2end));
    if saccv1 > threshes.blinkv1 && ...
            saccv2 > threshes.blinkv1
        % Case 1 (Note that Case 1b requires no action)
        if ~(saccstartevts(saccnum+1) > sacc1endevt )
            % Case 1a:
            msglist{end+1} = sprintf(...
                'lfp_saccadeCleanup:disordered2. Removing disordered nested start-end pair at %.6f - %.6f\n', ...
                lfp_Events(saccstartevts(saccnum+1), 1), ...
                lfp_Events(sacc1endevt, 1) );
            toremove([saccstartevts(saccnum+1) sacc1endevt]) = true;
        elseif mean(lfp_Samples{velo}( sacc1end : sacc2start )) ...
                > threshes.sacc
            % Case 1c:
            msglist{end+1} = sprintf( ...
                'lfp_saccadeCleanup:merge. Merging saccade end-start pair at %.6f - %.6f\n', ...
                lfp_Events(sacc1endevt, 1), ...
                lfp_Events(saccstartevts(saccnum+1), 1) );
            toremove([saccstartevts(saccnum+1) sacc1endevt]) = true;
            % also remove the fixation(s) if there are any between
            if any(ismember( ...
                    lfp_Events(saccstartevts(saccnum+1):sacc1endevt, 2), ...
                    [FixStart FixEnd BadFixStart BadFixEnd] ))
                fixstartevts = find(ismember( ...
                    lfp_Events(saccstartevts(saccnum+1):sacc1endevt, 2), ...
                    [FixStart BadFixStart] )) + saccstartevts(saccnum+1) - 1;
                fixendevts = find(ismember( ...
                    lfp_Events(saccstartevts(saccnum+1):sacc1endevt, 2), ...
                    [FixEnd BadFixEnd] )) + saccstartevts(saccnum+1) - 1;
                if length(startevts) ~= length(endevts)
                    error('lfp_saccadeCleanup:mismatchedfix', ...
                        'Fixation starts and ends do not match: %s', ...
                        dg_thing2str(lfp_Events(ismember( ...
                    lfp_Events(saccstartevts(saccnum+1):sacc1endevt, 2), ...
                    [FixStart FixEnd BadFixStart BadFixEnd] )), 1 ));
                end
                toremove([fixstartevts fixendevts]) = true;
            end
        end
    elseif saccv1 > threshes.blinkv1
        % Case 2: sacc2 is low-v
        if sqrt( (lfp_Samples{x}(sacc2end) ...
                - lfp_Samples{x}(sacc2start))^2 ...
                + (lfp_Samples{y}(sacc2end) ...
                - lfp_Samples{y}(sacc2start))^2 ) < 30
            % too small
            msglist{end+1} = sprintf( ...
                'lfp_saccadeCleanup:case2.  Removing sacc at %.6f - %.6f\n', ...
                lfp_Events(saccstartevts(saccnum+1), 1), ...
                lfp_Events(sacc2endevt, 1) );
            toremove([saccstartevts(saccnum+1) sacc2endevt]) = true;
            % if there is a fixation before and after the
            % saccade, with no intervening blinks or rec breaks, merge them
            % into one.  
            result = findfix2mergeevts(saccstartevts, ...
                saccnum+1, sacc2endevt);
            if ~isempty(result)
                fix2mergeevts(end+1, :) = result;
            end
        end
    elseif saccv2 > threshes.blinkv1
        % other Case 2: sacc1 is low-v
        if sqrt( (lfp_Samples{x}(sacc1end) ...
                - lfp_Samples{x}(sacc1start))^2 ...
                + (lfp_Samples{y}(sacc1end) ...
                - lfp_Samples{y}(sacc1start))^2 ) < 30
            % too small
            msglist{end+1} = sprintf( ...
                'lfp_saccadeCleanup:otherCase2. Removing sacc at %.6f - %.6f\n', ...
                lfp_Events(saccstartevts(saccnum), 1), ...
                lfp_Events(sacc1endevt, 1) );
            toremove([saccstartevts(saccnum) sacc1endevt]) = true;
            result = findfix2mergeevts(saccstartevts, ...
                saccnum, sacc1endevt);
            if ~isempty(result)
                fix2mergeevts(end+1, :) = result;
            end
        end
    else
        % Case 3a.	Tentative rule - remove all in this instance… ???
        msglist{end+1} = sprintf( ...
            'lfp_saccadeCleanup:case3a. Removing saccs at %.6f - %.6f and  %.6f - %.6f\n', ...
            lfp_Events(saccstartevts(saccnum), 1), ...
            lfp_Events(sacc1endevt, 1), ...
            lfp_Events(saccstartevts(saccnum+1), 1), ...
            lfp_Events(sacc2endevt, 1) );
        toremove([ saccstartevts(saccnum) sacc1endevt ...
            saccstartevts(saccnum+1) sacc2endevt ]) = true;
        % Note that the following usually (or perhaps always) results in
        % the middle fixation's being listed as the following fixation in
        % one row of fix2mergeevts, and then again as the preceding
        % fixation in the next row.
        result = findfix2mergeevts(saccstartevts, ...
            saccnum, sacc1endevt);
        if ~isempty(result)
            fix2mergeevts(end+1, :) = result;
        end
        result = findfix2mergeevts(saccstartevts, ...
            saccnum+1, sacc2endevt);
        if ~isempty(result)
            fix2mergeevts(end+1, :) = result;
        end
    end
end
lfp_log(cat(2, msglist{:}));

% Do the actual removals before proceeding to the next phase of cleanup.
% Doing the removals makes it necessary to recompute all the event indices.
% Worse than that, it also entails doing the fixation merges BEFORE doing
% the removals so that their event indices don't get invalidated.
%
% Iterate through fix2mergeevts marking events to remove and changing IDs
% of events that aren't toremove; this should automatically handle the note
% to Case 3a gracefully.
for k = 1:size(fix2mergeevts, 1)
    lfp_Events(fix2mergeevts(k, 1), 2) = MergedFixStart;
    lfp_Events(fix2mergeevts(k, 4), 2) = MergedFixEnd;
    toremove(fix2mergeevts(k, 2:3)) = true;
end
lfp_Events(toremove, :) = [];
saccstartevts = find(ismember(lfp_Events(:,2), saccEvtIDs(:,1)));
saccendevts = find(ismember(lfp_Events(:,2), saccEvtIDs(:,2)));

% Re-initialize lists whose indices were invalidated by the removal:
toremove(:) = false;
fix2mergeevts = zeros(0, 4); % start, end to remove, start to remove, end

%ii. Remove the whole saccade and associated markers if the start is >
%  maxvsaccstart and it comes after a blink or a rec break without any
%  intervening EyeEvents.
saccstarts = lfp_time2index(lfp_Events(saccstartevts, 1));
is_highv = lfp_Samples{velo}(saccstarts) > maxvsaccstart ;
msglist = {sprintf('suspiciously timed high-v saccstarts:\n')};
for saccnum = reshape(find(is_highv), 1, [])
    % Search forwards from sacc start for sacc end:
    evtidx = saccstartevts(saccnum);
    while evtidx <= size(lfp_Events, 1) && ~ismember( ...
            lfp_Events(evtidx, 2), saccEvtIDs(:,2) )
        evtidx = evtidx + 1;
    end
    if ~ismember(lfp_Events(evtidx, 2), saccEvtIDs(:,2))
        msg = sprintf( 'The high-v sacc at %.6f has no end.', ...
            lfp_index2time(saccstartevts(saccnum)) );
        lfp_log(['lfp_saccadeCleanup: ' msg]);
        warning('lfp_saccadeCleanup:nosaccend2', ...
             '%s', msg );
         continue
    end
    sacc1endevt = evtidx;
    evtidx = saccstartevts(saccnum);
    wrongrecseg = lfp_RecSegments( ...
        lfp_Events(evtidx,1) >= lfp_RecSegments(:,1) & ...
        lfp_Events(evtidx,1) <= lfp_RecSegments(:,2), :);
    recseg = lfp_RecSegments( ...
        lfp_time2index(lfp_Events(evtidx,1)) >= lfp_RecSegments(:,1) & ...
        lfp_time2index(lfp_Events(evtidx,1)) <= lfp_RecSegments(:,2), :);
    if isempty(recseg)
        % The event does not belong to any recorded segment, so we treat
        % the event as if it belongs to the preceding record segment
        recidx = find(lfp_Events(evtidx,1) >= lfp_RecSegments(:,1));
        if isempty(recidx)
            % There isn't any preceding segment, either.  Treat as if it
            % belongs to first segment.
            msg = sprintf('The high-v saccade at %.6f precedes all recorded segments', ...
                lfp_Events(evtidx,1) );
            lfp_log(['lfp_saccadeCleanup: ' msg]);
            warning('lfp_saccadeCleanup:norecseg', ...
                '%s', msg );
            recseg = lfp_RecSegments(1,:);
        else
            msg = sprintf( 'The high-v saccade at %.6f is not in any recorded segment', ...
                lfp_Events(evtidx,1) );
            lfp_log(['lfp_saccadeCleanup: ' msg]);
            warning('lfp_saccadeCleanup:norecseg2', ...
                '%s', msg );
            recseg = lfp_RecSegments(recidx(end),:);
        end
    end
    % search backwards to find the previous EyeEvent; skip the current sacc
    % start first:
    evtidx = evtidx - 1;
    while evtidx > 0 && ...
            ~ismember(lfp_Events(evtidx,2), EyeEvents)
        evtidx = evtidx - 1;
    end
    % If there is no such event (evtidx == 0), the sacc is the first event
    % in the fragment and should be removed; otherwise, if it's in the
    % previous rec segment, or if it's a blink end, remove the saccade.
    removethis = evtidx == 0 || ...
        ( lfp_time2index(lfp_Events(evtidx,1)) < recseg(1) || ...
        lfp_Events(evtidx,2) == BlinkEnd );
    if removethis
        msglist{end+1} = sprintf( ...
                'lfp_saccadeCleanup:highv.  Removing sacc at %.6f - %.6f\n', ...
                lfp_Events(saccstartevts(saccnum), 1), ...
                lfp_Events(sacc1endevt, 1) );
        toremove([saccstartevts(saccnum) sacc1endevt]) = true;
        result = findfix2mergeevts(saccstartevts, ...
            saccnum, sacc1endevt);
        if ~isempty(result)
            fix2mergeevts(end+1, :) = result;
        end
    end
    if isempty(wrongrecseg)
        msglist{end+1} = sprintf( ...
            'Bug 20081109 would have crashed on saccade at  %.6f - %.6f\n', ...
            lfp_Events(saccstartevts(saccnum), 1), ...
            lfp_Events(sacc1endevt, 1) );
    else
        wrongremove = evtidx == 0 || ...
            ( lfp_time2index(lfp_Events(evtidx,1)) < wrongrecseg(1) || ...
            lfp_Events(evtidx,2) == BlinkEnd );
        if removethis == wrongremove
            msglist{end+1} = sprintf( ...
                'Bug 20081109 made no difference to saccade at  %.6f - %.6f\n', ...
                lfp_Events(saccstartevts(saccnum), 1), ...
                lfp_Events(sacc1endevt, 1) );
        elseif wrongremove
            msglist{end+1} = sprintf( ...
                'Bug 20081109 would have removed saccade at  %.6f - %.6f\n', ...
                lfp_Events(saccstartevts(saccnum), 1), ...
                lfp_Events(sacc1endevt, 1) );
        else
            msglist{end+1} = sprintf( ...
                'Bug 20081109 have failed to remove saccade at  %.6f - %.6f\n', ...
                lfp_Events(saccstartevts(saccnum), 1), ...
                lfp_Events(sacc1endevt, 1) );
        end
    end
end
lfp_log(cat(2, msglist{:}));

% Iterate through fix2mergeevts marking events to remove and changing IDs
% of events that aren't toremove; this should automatically handle the note
% to Case 3a gracefully.
for k = 1:size(fix2mergeevts, 1)
    lfp_Events(fix2mergeevts(k, 1), 2) = MergedFixStart;
    lfp_Events(fix2mergeevts(k, 4), 2) = MergedFixEnd;
    toremove(fix2mergeevts(k, 2:3)) = true;
end

lfp_Events(toremove, :) = [];
% easy way to invoke lfp_createTrialIndex, etc:
lfp_createEvents(@(a) a, 666, []);


function result = findfix2mergeevts(saccstartevts, saccnum, saccendevt)
% Finds fixations preceding saccstartevts(saccnum) and following
% saccendevt without any intervening physical eye events (see
% <stopsearchIDs>). 
lfp_declareGlobals;
lfp_getEvtIDs;
result = [];
stopsearchIDs = [
    BlinkStart
    BlinkEnd
    SaccStart
    SaccEnd
    CplxSaccStart
    CplxSaccEnd
    BadFixStart
    BadFixEnd
    BTSaccStart
    BTSaccEnd    ];
fixstartIDs = [FixStart BlinkFixStart MergedFixStart];
fixendIDs = [FixEnd BlinkFixEnd MergedFixEnd];
% First find the fixation before:
evtidx = saccstartevts(saccnum) - 1;
while evtidx > 1 ...
        && ~ismember(lfp_Events(evtidx,2), fixendIDs) ...
        && ~ismember(lfp_Events(evtidx,2), stopsearchIDs)
    evtidx = evtidx - 1;
end
if evtidx == 0 || ...
        evtidx == 1 && ~ismember(lfp_Events(evtidx,2), fixendIDs)
    % Hit beginning of events without finding a fix end:
    return
end
if evtidx && ~ismember(lfp_Events(evtidx,2), stopsearchIDs)
    fix1endevt = evtidx;
    while evtidx > 1 ...
            && ~ismember(lfp_Events(evtidx,2), fixstartIDs) ...
            && ~ismember(lfp_Events(evtidx,2), stopsearchIDs)
        evtidx = evtidx - 1;
    end
    if evtidx && ismember(lfp_Events(evtidx,2), stopsearchIDs) ...
            || evtidx==1 && ~ismember(lfp_Events(evtidx,2), fixstartIDs)
        % found the fix end, but not the fix start:
        return
    elseif evtidx
        fix1startevt = evtidx;
        % Second, find fixation after:
        evtidx = saccendevt + 1;
        while evtidx < size(lfp_Events,1) ...
                && ~ismember(lfp_Events(evtidx,2), fixstartIDs) ...
                && ~ismember(lfp_Events(evtidx,2), stopsearchIDs)
            evtidx = evtidx + 1;
        end
        if evtidx == size(lfp_Events,1) ...
                && ~ismember(lfp_Events(evtidx,2), fixstartIDs)
            % Hit end of events without finding fix start
            return
        end
        if evtidx < size(lfp_Events,1) ...
                && ~ismember(lfp_Events(evtidx,2), stopsearchIDs)
            fix2startevt = evtidx;
            while evtidx < size(lfp_Events,1) ...
                    && ~ismember(lfp_Events(evtidx,2), fixendIDs) ...
                    && ~ismember(lfp_Events(evtidx,2), stopsearchIDs)
                evtidx = evtidx + 1;
            end
            if evtidx < size(lfp_Events,1) ...
                    && ismember(lfp_Events(evtidx,2), stopsearchIDs) ...
                    || evtidx == size(lfp_Events,1) ...
                    && ~ismember(lfp_Events(evtidx,2), fixendIDs)
                % found the fix start, but not the fix end; just abort
                return
            elseif evtidx < size(lfp_Events,1)
                fix2endevt = evtidx;
                result = ...
                    [fix1startevt fix1endevt fix2startevt fix2endevt];
            end
        end
    end
end
