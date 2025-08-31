function lfp_removeEyeJunk(velo, evtcount, evts2remove, threshes, vpeakthresh)
% Post-finally, a scan is done to detect and remove post-blink junk, which
% is defined as follows.  The post-blink interval is considered to begin at
% the first sample after blink end, and to end at the last sample before 
% the first saccade that crosses 2*blinkv1, or blink start, following the
% blink end.  Event frequency is read from the CSC channel <evtcount>. If 
% the event frequency >= <threshes.evtfreq> events/s anywhere in the
% post-blink interval, then ALL eye events in the post-blink interval are
% considered junk, and are replaced by a Junk Start on the first sample of
% the post-blink interval and Junk End on the last sample.  <velo> is CSC
% channel containing velocity trace.
% <threshes> is as for lfp_parseEyeTraces6, but in addition it must have an
% <evtfreq> field.
% 01-Oct-2008 DG added another category of junk identified by Theresa,
% namely cases where the first saccade after the blink is fitted to a
% velocity trace that doesn't have a real maximum, as defined by
% lfp_findmax (warning 'lfp_findmax:endpoint' is suppressed during this
% test).  This also entailed getting more picky about unmatched events;
% each category of event (sacc, cplxsacc, blink) is independently zipped to
% ensure that the final lists of all starts and ends contain exclusively
% paired events.  Unpaired events are ignored.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;

maxfixv = 1.75 * threshes.sacc;
minsaccint = 0.12;  % no relation to minsaccint in lfp_parseEyeTraces7
switch vpeakthresh
    case 'blinkv1'
        vpeakthreshval = threshes.blinkv1;
    case 'maxfixv'
        vpeakthreshval = maxfixv;
    otherwise
        error('lfp_removeEyeJunk:badvpeakthresh', ...
            ['The option "' vpeakthresh '" is not recognized.'] );
end

% Find eye events
saccstartevtidx = find(lfp_Events(:,2) == SaccStart);
cplxsaccstartevtidx = find(lfp_Events(:,2) == CplxSaccStart);
blinkstartevtidx = find(lfp_Events(:,2) == BlinkStart);
saccendevtidx = find(lfp_Events(:,2) == SaccEnd);
cplxsaccendevtidx = find(lfp_Events(:,2) == CplxSaccEnd);
blinkendevtidx = find(lfp_Events(:,2) == BlinkEnd);

if length(saccstartevtidx) ~= length(saccendevtidx)
    warning('lfp_removeEyeJunk:unpairedsacc', ...
        'SaccStarts: %d, SaccEnds: %d', ...
        length(saccstartevtidx), length(saccendevtidx) );
end
if length(cplxsaccstartevtidx) ~= length(cplxsaccendevtidx)
    warning('lfp_removeEyeJunk:unpairedcplxsacc', ...
        'SaccStarts: %d, SaccEnds: %d', ...
        length(cplxsaccstartevtidx), length(cplxsaccendevtidx) );
end
if length(blinkstartevtidx) ~= length(blinkendevtidx)
    warning('lfp_removeEyeJunk:unpairedblink', ...
        'SaccStarts: %d, SaccEnds: %d', ...
        length(blinkstartevtidx), length(blinkendevtidx) );
end

allstartevtidx = sort([saccstartevtidx
    cplxsaccstartevtidx 
    blinkstartevtidx]);
allendevtidx = sort([saccendevtidx 
    cplxsaccendevtidx 
    blinkendevtidx]);
blinkstarts = lfp_time2index(lfp_Events(blinkstartevtidx,1));
blinkends = lfp_time2index(lfp_Events(blinkendevtidx,1));
allstarts = lfp_time2index(lfp_Events(allstartevtidx, 1));
allends = lfp_time2index(lfp_Events(allendevtidx, 1));
s = warning('off', 'lfp_findmax:endpoint');
saccvpeaks = lfp_findmax(velo, allstarts, allends);
warning(s.state, 'lfp_findmax:endpoint');

% Find post-blink junk and replace it with Junk Start / Junk End markers.
junkstarts = zeros(0,1);
junkends = zeros(0,1);
h = waitbar(0,'Finding junk periods...');
recsegidx = 1;
% <firststarts> is sample #s of first sacc/blink start after each blink:
firststarts = zeros(size(blinkends)); 
for blinknum = 1:length(blinkends)
    while blinkends(blinknum) > lfp_RecSegments(recsegidx,2)
        if recsegidx < size(lfp_RecSegments,1)
            recsegidx = recsegidx + 1;
        else
            error('lfp_removeEyeJunk:bug', ...
                'Bug: there is a blink end at sample %d, which is beyond the end of the last rec segment', ...
                blinkends(blinknum) );
        end
    end
    postblinkstart = blinkends(blinknum) + 1;
    postblinkend = NaN;
    postblinkstartsidx = find(allstarts > postblinkstart);
    if isempty(postblinkstartsidx)
        % There are no more sacc/blink starts, so we're done.
        break
    end
    firststarts(blinknum) = postblinkstartsidx(1);
    for saccnum = postblinkstartsidx(1):length(allstarts)
        if allstarts(saccnum) > lfp_RecSegments(recsegidx,2)
            % That's it for this recording segment; pack it in
            postblinkend = lfp_RecSegments(recsegidx,2);
            break
        end
        if ismember(allstarts(saccnum), blinkstarts) ...
                || any(lfp_Samples{velo}( ...
                allstarts(saccnum):allends(saccnum) ) ...
                > 2*vpeakthreshval )
            % Found the first big-enough saccade
            postblinkend = allstarts(saccnum) - 1;
            break
        end
    end
    if isnan(postblinkend)
        % All saccades were too low-v to count, so we're done.
        break
    else
        if any( lfp_Samples{evtcount}(postblinkstart:postblinkend) ...
                >= threshes.evtfreq )
            % It's junk
            junkstarts(end+1, 1) = postblinkstart;
            junkends(end+1, 1) = postblinkend;
        end
    end
    waitbar(blinknum/length(blinkends), h );
end
close(h);
% Remove the junk events:
toremove = false(size(lfp_Events,1), 1);
junkstartTS = zeros(size(junkstarts));
junkendTS = zeros(size(junkstarts));
isAnyStartEvt = ismember(lfp_Events(:,2), [BlinkStart
    SaccStart
    FixStart
    MrgStart
    FitStart
    OldFixStart
    CplxSaccStart
    BadFixStart
    BlinkFixStart
    BTSaccStart ]);
isAnyEndEvt = ismember(lfp_Events(:,2), [BlinkEnd
    SaccEnd
    FixEnd
    MrgEnd
    FitEnd
    OldFixEnd
    CplxSaccEnd
    BadFixEnd
    BlinkFixEnd
    BTSaccEnd ]);
isMidFix = lfp_Events(:,2) == MidFix;
msgparts{1} = sprintf('Marked the following runs of events as junk:\n');
h = waitbar(0,'Marking junk periods...');
for k = 1:length(junkstarts)
    junkstartTS(k) = lfp_index2time(junkstarts(k));
    junkendTS(k) = lfp_index2time(junkends(k));
    % Truncation error can be a problem here, so we fudge by 0.5 us:
    isjunk = lfp_Events(:,1) >= junkstartTS(k) - 0.5e-6 & ...
        lfp_Events(:,1) <= junkendTS(k) + 0.5e-6 & ...
        ismember(lfp_Events(:,2), evts2remove);
    % There can be cases where events that should belong to the junk don't
    % get included because a series of unfortunate incidents has placed
    % them on the same sample as the post-blink-defining blink end or sacc
    % start.  These must also be added to the junk list.
    blinkendTS = junkstartTS(k) - lfp_SamplePeriod;
    saccstartTS = junkendTS(k) + lfp_SamplePeriod;
    isprejunkevt = abs(lfp_Events(:,1) - blinkendTS) < lfp_SamplePeriod/2;
    ispostjunkevt = abs(lfp_Events(:,1) - saccstartTS) < lfp_SamplePeriod/2;
    isjunk(isprejunkevt & isAnyStartEvt) = true;
    isjunk(ispostjunkevt & isAnyEndEvt) = true;
    if any((isprejunkevt | ispostjunkevt) & isMidFix)
        error('lfp_removeEyeJunk:midfix', ...
            'There is a midfix on a post-blink defining blink end or sacc start');
    end
    toremove(isjunk) = true;
    junkidx = find(isjunk);
    if isempty(junkidx)
        msgparts{end+1} = sprintf( ...
            'No events to delete in %.6f - %.6f interval\n===', ...
            junkstartTS(k), junkendTS(k) );
    else
        msg = '';
        for eix = reshape(junkidx, 1, [])
            msg = sprintf('%s%.6f %8d %s\n', msg, lfp_Events(eix,1), ...
                lfp_Events(eix,2), lfp_EventNames{lfp_Events(eix,2)} );
        end
        msgparts{end+1} = sprintf('%s===\n', msg);
    end
    waitbar(k/length(junkstarts), h );
end
close(h);
msg = cat(2, msgparts{:});
lfp_log(msg);

% find the firststarts that are NOT in any junk period and test them for
% lfp_findmax badness; if bad, add the saccstart and saccend to <toremove>.
h = waitbar(0,'Finding saccs with bad max...');
for blinknum = 1:length(blinkends)
    saccidx = firststarts(blinknum);
    % saccidx can be zero if there were no post-blink saccades
    if saccidx && ( isempty(junkstarts) ...
            || (~any( allstarts(saccidx) >= junkstarts ...
            & allstarts(saccidx) <= junkends )) )
        % allstarts and allends "should" be paired
        if allends(saccidx) < allstarts(saccidx) ...
                || (saccidx > 1 && allends(saccidx - 1) > allstarts(saccidx))
            warning('lfp_removeEyeJunk:unpairedevt', ...
                'The event at t = %.6f is matched to an earlier end or is earlier than the preceding end', ...
                lfp_index2time(allstarts(saccidx)) );
        end
        % check for lfp_findmax badness:
        if saccvpeaks(saccidx) == allstarts(saccidx)
            toremove(allstartevtidx(saccidx)) = true;
            toremove(allendevtidx(saccidx)) = true;
        end
    end
    waitbar(blinknum/length(blinkends), h );
end
close(h);

% The call to lfp_createEvents does the index recalculation for the
% deletions and the insertions all at once:
lfp_Events(toremove, :) = [];
lfp_createEvents(@(a) a, [JunkStart, JunkEnd], {junkstartTS, junkendTS});
