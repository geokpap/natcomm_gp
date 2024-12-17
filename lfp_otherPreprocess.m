function [params, events, badtrials] = lfp_otherPreprocess(events)
% Timestamps in <events> are still in microseconds at this point.
% For 'ken' data:
%   Put reward values into <params> col. 1 and aversive stim values into
% <params> col. 2; returns <events> with those parameter values removed.
% The values are NaN for trials in which the parameter is not recorded.
% Legitimate events are considered to start with the first PRE_TRIAL
% marker; preceding events are deleted with a warning and log message.
%
% For 'georgios' data:
%   Consider events [LEFT_CROSS RIGHT_CROSS] to serve double-duty as escape
% markers that are always followed by exactly two trial params.  For
% standard conflict trials those are respectively reward and aversion
% magnitudes, and they go in <params> col. 1 and col. 2 respectively.
%   Georgios does not use anything that can be used as trial end events.
% Therefore, these have to be inserted.
%   2-May-2019: Georgios has revised his event codes s.t. all event with
% TTL codes 0 - 200 represent outcome magnitudes, and all other events
% have codes >200.
%   7-Feb-2022 updated 6-Aug-2022: <badtrials> has now been repurposed, and
% <isbadtrial> no longer exists.  The previous meaning of <badtrials>
% became obsolete and incorrect when I added the trial repair code for
% those trials.
%   6-Aug-2022: <badtrials> is passed back to 'lfp_readEvents', which
% passes it back its caller.  As of this writing, the only caller of
% 'lfp_readEvents' that uses it is 'lfp_read2', which uses it as the
% initial value of <lfp_BadTrials>. The other callers ('lfp_add',
% 'lfp_fragmentFiles') just raise an error if it is not empty.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;

badtrials = [];

switch lfp_SetupType
    case 'ken'
        % As of Oct 4 2012, Ken does not use 0 as a TTL code for anything,
        % not even reward or aversion values.  Therefore all events with
        % TTL value 0 are accidental triggerings.
        iszero = events(:,2)==0;
        events(iszero,:) = [];
        % Ken in effect uses two escaped single-element BD sections (reward
        % value and aversion value) that precede the trial, together with a
        % bunch of BD data emitted at 200 ms intervals during the
        % "pre-trial" period.  However, the non-escaped BD values are
        % unique event codes, so there is nothing to gain from unpacking
        % them into <params>.  CONTROL_TASK specifies both ENCODE_REWARD
        % and ENCODE_AVERSI, so presumably the latter gives the size of the
        % "alternate reward" and it is delivered when the monk chooses
        % the AVOID action.
        %   params: col 1 = reward value; col. 2 = aversion or alt.rwd.
        pretrials = find(events(:,2) == PRE_TRIAL);
        % It seems to be routine that the recording doesn't start until the
        % middle of a trial, so these are junk events:
        if pretrials(1) ~= 1
            numbogus = pretrials(1) - 1;
            msg = sprintf(...
                'Deleting %d events before the first PRE_TRIAL marker; see log for event details', ...
                numbogus);
            warning( 'lfp_otherPreprocess:ken3', '%s', msg);
            for k = 1:numbogus
                msg = sprintf( '%s\n%12.0f %3d', ...
                    msg, events(k,1), events(k,2) );
            end
            lfp_log(msg);
            events(1:numbogus, :) = [];
            pretrials = find(events(:,2) == PRE_TRIAL);
        end
        endpretrials = find(events(:,2) == START_TRIAL);
        intertrials = find(events(:,2) == START_INTER_TRIAL);
        if endpretrials(1) < pretrials(1)
            warning( 'lfp_otherPreprocess:ken4', ...
                'The first END_PRE_TRIAL is before the first PRE_TRIAL' );
            endpretrials(endpretrials(1) < pretrials(1)) = [];
        end
        if pretrials(end) > intertrials(end)
            % These are trailing junk events:
            warning( 'lfp_otherPreprocess:ken5', ...
                'The last PRE_TRIAL is after the last ITI, deleting trailing events' );
            pretrials(pretrials > intertrials(end)) = [];
            endpretrials(endpretrials > intertrials(end)) = [];
            events(intertrials(end)+1 : end, :) = [];
        end
        [pairs, extraL, extraR] = dg_zip(pretrials, endpretrials);
        if ~isempty(extraL) || ~isempty(extraR)
            msg = 'There are mismatched pre_trials/end_pre_trials:';
            for k = 1:length(extraL)
                msg = sprintf( '%s\n  extra PRE_TRIAL at t=%.6f sec', ...
                    msg, events(pretrials(extraL(k)),1)*1e-6 );
            end
            for k = 1:length(extraR)
                msg = sprintf( '%s\n  extra END_PRE_TRIAL at t=%.6f sec', ...
                    msg, events(endpretrials(extraR(k)),1)*1e-6 );
            end
            warning('lfp_otherPreprocess:ken1', msg);
        end
        pretrials = pretrials(pairs(:,1));
        endpretrials = endpretrials(pairs(:,2));
        % The event following evt 4 is the reward value, and the event
        % following evt 5 is aversive value.
        rewardvalidx = find(events(:,2) == 4) + 1;
        aversvalidx = find(events(:,2) == 5) + 1;
        params = cell(1, length(pretrials));
        evts2delete = [];
        for t = 1:length(pretrials)
            rwdidx = rewardvalidx( ...
                rewardvalidx > pretrials(t) & rewardvalidx < endpretrials(t) );
            if isempty(rwdidx)
                params{t}(1, 1) = NaN;
            else
                if length(rwdidx) > 1
                    warning('lfp_otherPreprocess:ken2', ...
                        'Trial %d has %d reward value events, using first', ...
                        t, length(rwdidx));
                    rwdidx = rwdidx(1);
                end
                params{t}(1, 1) = events(rwdidx,2);
                evts2delete(end+1) = rwdidx; %#ok<AGROW>
            end
            avridx = aversvalidx( ...
                aversvalidx > pretrials(t) & aversvalidx < endpretrials(t) );
            if isempty(avridx)
                params{t}(1, 1) = NaN;
            else
                if length(avridx) > 1
                    warning('lfp_otherPreprocess:ken2', ...
                        'Trial %d has %d aversion value events, using first', ...
                        t, length(avridx));
                    avridx = avridx(1);
                end
                params{t}(1, 2) = events(avridx,2);
                evts2delete(end+1) = avridx; %#ok<AGROW>
            end
        end
        events(evts2delete, :) = [];
    case 'georgios_ApAv_rev_only'
        error('lfp_otherPreprocess:nocode', ...
            'Need to write code for lfp_SetupType ''%s''.', ...
            lfp_SetupType);
    case {'georgios_ApAv' 'georgios_ApAvApAp' 'georgios_ApAp'}
        if events(1, 2) == FixPtOn
            % This is a trial fragment, and the fact that it starts with
            % FixPtOn will break subsequent code, so we just delete the
            % first event:
            events(1, :) = [];
        end

        % Pull out trial params.
        %   Updated 4-May-2022: There are two types of trials, "free"
        % choice which have two outcome params, one for airpuff magnitude
        % and one for reward magnitude, ["in that order" was wrong, they
        % are in REVERSE order!!!  I.e., same order as they are SPECIFIED
        % in a choice trial. -DG 20240905]; FORCED trials that have only
        % one outcome param.  On FORCED trials, the irrelevant outcome
        % has either magnitude 0 or the same magnitude repeated.
        % However, because 'georgios_ApAv_rev_only' does not <usestrobe>,
        % any repeated events fail to record, so we have to assume that
        % FORCED trials can have either 2 or 1 outcome magnitude.  The
        % possible event sequences are diagrammed in
        % "EvtTransitionNetwork.pptx", but see "Just to make it harder"
        % below. The inevitable Monkey Logic events 9 and 18 are trial
        % start and trial end respectively.  All events up to FixPtOn are
        % bookkeeping events, and the monkey does not need to do anything
        % until FixPtOn, so FixPtOn is always the first physical event of
        % every trial.  There are exactly 1200 trials in a completed
        % session. The old trial start event 255 is actually a manually
        % inserted event, so it may be missing or inserted at the wrong
        % time.
        %   Just to make it harder:
        % 1. There are not necessarily always two outcome magnitudes,
        % either.  ML sends the same outcome magnitude twice in a row for
        % FORCED trials, and may also do that for the occasional "choice"
        % trial.  The second events fails to trigger TTL recording by
        % Cheetah.
        % 2. Not all events 9 or 18 are necessarily trail start or end
        % events, because those are possible outcome magnitudes.  That
        % means that FixPtOn is actually the ONLY reliable trial event.
        % Also, for unclear reasons, events 9 or 18 are sometimes missing.
        % 3. Event 255 "TrialStart" is recorded manually and may also be
        % missing.

        if isequal(lfp_SetupType, 'georgios_ApAvApAp')
            % Before we can even get started, we need to divide the session
            % into ApAv blocks and ApAp blocks because event 213 is
            % double-booked: FixPtOnApAp and RWD_ON_FORCED.  Therefore we
            % cannot simply search for events [FixPtOn FixPtOnApAp].
            allApApevts = [FixPtOnApAp RWD_OFF_FORCEDRedApAp ...
                RWD_OFF_FORCEDYelApAp ...
                choiceCueOnApAp crossAcqApAp ...
                crossOffRedApAp forcedRedOffApAp forcedRedOnApAp ...
                forcedYelOffApAp forcedYelOnApAp rewardOffRedApAp ...
                rewardOffYelApAp rewardOnRedApAp rewardOnYelApAp ...
                squareAcqYelApAp targetsOffYelApAp targetsOnApAp];
            isApApevt = ismember(events(:,2), allApApevts);
            if sum(isApApevt) == 0
                error('lfp_otherPreprocess:notApAvApAp', ...
                    'This session contains no ApAp events.');
            end
            smooApApevt = conv(double(isApApevt), ones(51,1)/51, 'same');
            % <isApAv> knows nothing of actual trial boundaries, but the
            % transitions between <true> and <false>, i.e. <endApAv> and
            % <startApAv>, occur either late in the last trial of the block
            % that is ending or early in the first trial of the block that
            % is beginning, so they are close enough to trial boundaries to
            % use them to distinguish event 213s that mean
            % 'RWD_ON_FORCEDApAv' from event 213s that mean 'FixPtOnApAp'.
            isApAv = smooApApevt < 1/50;
            % Check to make sure blocks are not too short (which would
            % indicate a failure in the criterion for <isApAv>). <endApAv>
            % is last ApAv evt in run, <endApAp> is last ApAp evt in run.
            endApAv = reshape( find(isApAv(1:end-1) & ~isApAv(2:end)), ...
                [], 1 ); 
            endApAp = reshape( find(~isApAv(1:end-1) & isApAv(2:end)), ...
                [], 1 ); 
            % We assume the the following means that there are only two
            % blocks in the session, ApAv followed by ApAp:
            if isempty(endApAp)
                endApAp = length(isApAv);
            end
            numblox = length(endApAv) + length(endApAp);
            if numblox < 4
                warning('lfp_otherPreprocess:numblox', ...
                    'There are only %d blocks in this session.', ...
                    numblox);
            end
            % First evts in runs:
            startApAv = endApAp + 1;
            startApAp = endApAv + 1;
            if isempty(startApAv) || startApAv(1) > endApAv(1)
                startApAv = [1; startApAv];
            end
            if startApAv(end) > endApAv(end)
                endApAv(end+1, 1) = length(isApAv);
            end
            durApAv = endApAv - startApAv + 1;
            durApAv(durApAv==0) = [];
            if any(durApAv < 20)
                error('lfp_otherPreprocess:durApAv', ...
                    'ApAv block is too short.');
            end
            if startApAp(1) > endApAp(1)
                startApAp = [1; startApAp];
            end
            if startApAp(end) > endApAp(end)
                endApAp(end+1, 1) = length(isApAv);
            end
            durApAp = endApAp - startApAp + 1;
            if any(durApAp < 20)
                error('lfp_otherPreprocess:durApAp', ...
                    'durApAp block is too short.');
            end
            % change all event 213s that are in ApAv blocks (formerly known
            % as 'RWD_ON_FORCEDApAv' to 210s (RWD_ON_FORCED):
            is213 = events(:, 2) == 213;
            events(is213 & isApAv, 2) = 210;
            % replace all event 250s (formerly known as
            % 'RWD_ON_FORCEDRedApAp') with 210s (RWD_ON_FORCED):
            events(events(:,2) == 250, 2) = 210;
            % Now we search for fixation point onsets as [FixPtOn
            % FixPtOnApAp]:
            fixptonidx = find(ismember(events(:, 2), ...
                [FixPtOn FixPtOnApAp]));
        else
            % <fixptonidx> is by definition the variable that identifies
            % trials. It depends on <events>.  Therefore, we defer any
            % edits to <events> until after we are finished with this
            % incarnation of <fixptonidx>.
            fixptonidx = find(events(:, 2) == FixPtOn);
        end
        
        % There can be a partially recorded trial at the beginning of the
        % session, so we need to check the events preceding the first
        % <fixptonidx> to make sure they represent intact trials. The
        % minimum possible number of non-zero TTL events preceding FixPtOn
        % is 2, and they must be "targposn / FORCED" and an outcome
        % magnitude in that order.  We work backwards to verify and if
        % necessary amputate.  See code following
        % "warning('lfp_otherPreprocess:missing'" for handling of a trial
        % fragment in the middle or at the end of the session.
        if isempty(fixptonidx)
            error('lfp_otherPreprocess:fixptonidx2', ...
                'Session contains no fixptonidx events.');
        end
        hasmagnitude = events(fixptonidx(1) - 1, 2) < 200;
        hastargposn = false;
        for evtidx = fixptonidx(1) - 2 : -1 : 1
            if ismember(events(evtidx, 2), [CROSS_LEFT CROSS_RIGHT FORCED])
                hastargposn = true;
            end
        end
        % Save this result for future reference when editing events:
        isbadfixpt1 = ~(hasmagnitude && hastargposn);
        badfixpt1idx = [];
        if isbadfixpt1
            % remove trial fragment from <fixptonidx>, and save its index
            % in the events list so that it can be deleted from <events>
            % later on:
            badfixpt1idx = fixptonidx(1);
            fixptonidx(1) = [];
        end
        numtrials = length(fixptonidx);
        if numtrials ~= 1200
            warning('lfp_otherPreprocess:numtrials2', ...
                'Found %d FixPtOn events (usually 1200).', numtrials);
        end
        
        % Extract offer magnitudes.  This must be done based on
        % <fixptonidx>.  By definition, a well-formed trial must have one
        % of [CROSS_LEFT CROSS_RIGHT FORCED] as either the second or third
        % event before FixPtOn, and the one or two events in between
        % [CROSS_LEFT CROSS_RIGHT FORCED] and FixPtOn must be in the 0-200
        % range that represents outcome magnitudes.  Any ill-formed trials
        % are marked as bad and added to <lfp_BadTrials> via <badtrials>,
        % with NaN outcome magnitudes in <lfp_TrialParams>.  For
        % well-formed trials, RWD magnitudes go in
        % <lfp_TrialParams{k}(:,1)> and AVE magnitudes in
        % <lfp_TrialParams{k}(:,2)>.
        isgood = ismember( ...
            events(fixptonidx - 2, 2), [CROSS_LEFT CROSS_RIGHT FORCED] ) ...
            | ismember( ...
            events(fixptonidx - 3, 2), [CROSS_LEFT CROSS_RIGHT FORCED] );
        badtrials = find(~isgood);
        morebadtrials = [];
        % At this point, <badtrials> is an index into <fixptonidx>, and is
        % therefore a suitable value to use to initialize <lfp_BadTrials>
        % in lfp_read2.  This means we are now committed to keeping every
        % FixPtOn that currently exists in <events> when we get around to
        % deleting the outcome magnitude events later on.  We must also
        % treat <fixptonidx> as immutable until then. Trying to keep track
        % of all that while still doing the computations in parallel
        % exceeds my working memory capacity, and there are only O(1000)
        % trials anyway, so we iterate over trials here.
        params = NaN(length(fixptonidx), 2);
        evts2del = zeros(0, 1);
        for trialnum = 1:length(fixptonidx)
            if ismember(trialnum, badtrials)
                % do nothing, <params> is already initialized to NaN.
            else
                if events(fixptonidx(trialnum) - 1, 2) < 0 || ...
                        events(fixptonidx(trialnum) - 1, 2) > 200
                    % What should be the second outcome magnitude is not
                    % an outcome magnitude, so this is another bad trial.
                    morebadtrials(end+1) = trialnum; %#ok<AGROW>
                else
                    params(trialnum, 2) = ...
                        events(fixptonidx(trialnum) - 1, 2);
                    evts2del(end+1, 1) = fixptonidx(trialnum) - 1; %#ok<AGROW>
                    if events(fixptonidx(trialnum) - 2, 2) < 0 || ...
                            events(fixptonidx(trialnum) - 2, 2) > 200
                        % There is only one outcome magnitude event, so we
                        % infer the first outcome magnitude was the same as
                        % the second.
                        params(trialnum, 1) = params(trialnum, 2);
                    else
                        % There were two different outcome magnitude events
                        params(trialnum, 1) = ...
                            events(fixptonidx(trialnum) - 2, 2);
                        evts2del(end+1, 1) = fixptonidx(trialnum) - 2; %#ok<AGROW>
                    end
                end
            end
        end
        params = mat2cell(params, ones(length(fixptonidx),1), 2);
        
        % Now we delete events (notably the param events) from the events
        % list. NaNs represent non-magnitude events that were in the
        % position of the airpuff magnitude.  We also delete the first
        % <fixptonidx> from <events> if it was a trial fragment.
        %   NOTE: this deletion invalidates <fixptonidx>, etc.  (It also
        % might invalidate <numtrials>, but we use that value to determine
        % whether the number of trials has in fact changed.)
        evts2del = [badfixpt1idx; evts2del];
        events(evts2del, :) = [];
        clear fixptonidx badfixpt1idx hasmagnitude
        % At this point all events that denote outcome magnitudes have been
        % moved into <params>.  Cheetah inserts events with zero for TTL at
        % the very beginning and end of the file, so those can be safely
        % ignored.  Any others are defects:
        event0s = find(events(:, 2) == 0);
        defects = setdiff(event0s, [1 2 size(events, 1)]);
        if isempty(defects)
            % All Cheetah, simply delete them:
            events(event0s, :) = [];
        else
            warning('lfp_otherPreprocess:georgios2', ...
                'Replacing event 0s with event %d.', event_zero);
            % <event_zero> must be defined by call to lfp_getEvtIDs:
            events(event0s, 2) = event_zero;
        end
        % We can now reliably find MonkeyLogic trial starts and ends:
        trialstarts = find(events(:, 2) == MLTrialStart);
        trialends = find(events(:, 2) == MLTrialEnd);
        % Make sure that recalculating <fixptonidx> doesn't mess up the
        % trial numbering expressed in <badtrials>.  This functions as an
        % "assert":
        if isequal(lfp_SetupType, 'georgios_ApAvApAp')
            if sum(ismember(events(:, 2), [FixPtOn FixPtOnApAp])) ...
                    ~= numtrials
                error('lfp_otherPreprocess:fixptonidx2', ...
                    'Events editing has altered the number of trials.');
            else
                fixptonidx = find(ismember( events(:, 2), ...
                    [FixPtOn FixPtOnApAp] ));
            end
        else
            if sum(events(:, 2) == FixPtOn) ~= numtrials
                error('lfp_otherPreprocess:fixptonidx', ...
                    'Events editing has altered the number of trials.');
            else
                fixptonidx = find(events(:, 2) == FixPtOn);
            end
        end
        % <trialstarts>, <fixptonidx>, and <trialends> are all indices into
        % <events>; therefore they should zip nicely to form triples of
        % [trialstarts, fixptonidx, trialends] with one triple per trial.
        [~, extratrialstarts, extrafixpts] = dg_zip(trialstarts, fixptonidx);
        [~, extrafixpts2, extratrialends] = dg_zip(fixptonidx, trialends);
        % When there is an interruption in recording, we will typically get
        % a bogus "trial" that consists of the end of one trial and the
        % beginning of another.  This is easy to detect if it results in
        % there being no FixPtOn.  If not, see "If there is a recording".
        trials2del = [];
        if ~isempty(extratrialstarts) || ~isempty(extratrialends)
            % This <if> construct must set <evts2del> in all cases.
            if isequal(extratrialstarts, extratrialends)
                % Handle the simple case where they are all from trials
                % that contain no FixPtOn, in which case they will be
                % equal.
                trials2del = extratrialstarts;
                evts2del = [];
            elseif isempty(extratrialstarts) && isequal(extratrialends, 1)
                % This could be either a missing trialstart at the
                % beginning of the session or a trial fragment.  We treat
                % it as a trial fragment for simplicity, and eliminate by
                % deleting the extra trial end event.  But we can't just
                % delete directly from <events> because that will
                % invalidate <trialstarts> and <trialends>, which we are
                % still using, so we put it in <evts2del>.
                evts2del = trialends(1);
                trialends(extratrialends) = [];
            elseif ( isempty(extrafixpts) || ...
                    all(extrafixpts) < length(fixptonidx) ) && ... 
                    ( isempty(extrafixpts2) || ...
                    all(extrafixpts2) < length(fixptonidx) ) ... 
                    && ~isempty(extratrialstarts) && all( ...
                    trialstarts(extratrialstarts) > fixptonidx(end) ) ...
                    && ~isempty(extratrialends) && all( ...
                    trialends(extratrialends) > fixptonidx(end) )
                % All three of these conditions are met:
                % 1. There are no extra <fixptonidx> at the end of
                % <fixptonidx>;
                % 2. the only extra <trialstarts> are after the last good
                % triple;
                % 3. the only extra <trialends> are after the last good
                % triple;
                % ...and that implies that there are just some extra ML
                % events hanging around at the end but the session is
                % otherwise fine, so the extra events are nothing to worry
                % about and can simply be deleted.
                evts2del = [
                    trialstarts(extratrialstarts)
                    trialends(extratrialends)
                    ];
            else
                error('lfp_otherPreprocess:extrastartends', ...
                    'This data error is too complicated to handle.');
            end
        else
            evts2del = [];
        end
        for trial2delidx = 1:length(trials2del)
            evts2del = [evts2del
                ( trialstarts(trials2del(trial2delidx)) ...
                : trialends(trials2del(trial2delidx)) )'
                ]; %#ok<AGROW> 
        end
        events(evts2del, :) = [];
        % That deletion invalidates these indices into <events>:
        clear fixptonidx trialstarts trialends defects event0s
        if isequal(lfp_SetupType, 'georgios_ApAvApAp')
            fixptonidx = find(ismember( events(:, 2), ...
                [FixPtOn FixPtOnApAp] ));
        else
            fixptonidx = find(events(:, 2) == FixPtOn);
        end
        trialstarts = find(events(:, 2) == MLTrialStart);
        trialends = find(events(:, 2) == MLTrialEnd);
        
        % If there are any missing Monkey Logic trialstarts or trialends,
        % then additional events must be fabricated to replace the missing
        % ones.  No attempt is made to replace missing FixPtOns; those
        % trials simply disappear without a trace.
        if ~isempty(extrafixpts) || ~isempty(extrafixpts2)
            warning('lfp_otherPreprocess:missing', ...
                'Missing required events; repairing.');
            newstartTS = zeros(0,1);
            newendTS = zeros(0,1);
            if ~isempty(extrafixpts)
                for extrafixidx = 1:length(extrafixpts)
                    % There is a trial start missing before each
                    % <extrafixpts>. If there is an IEI > 8 s within 5
                    % events of the "extra" FixPtOn, then insert the
                    % trial start at its end. Otherwise, simply start the
                    % trial right before the "extra" FixPtOn.
                    for evtidx = fixptonidx(extrafixpts(extrafixidx))
                        % <lookback> is needed to handle events that are
                        % within <5 events of the first event in the file:
                        lookback = min(evtidx - 1, 5);
                        IEIs = diff(events(evtidx + (-lookback:0), 1));
                        longidx = find(IEIs > 8, 1, 'last');
                        if isempty(longidx)
                            newstartTS(end+1, 1) = ...
                                events(evtidx, 1) - 10; %#ok<AGROW>
                        else
                            newstartTS(end+1, 1) = ...
                                events(evtidx - lookback + longidx - 1, 1) - 10; %#ok<AGROW>
                        end
                    end
                end
                warning('lfp_otherPreprocess:newstarts', ...
                    'Inserting remedial trial start(s) for trial(s) %s at %s', ...
                    dg_thing2str(extrafixpts), ...
                    dg_thing2str(newstartTS/1e6));
            end
            if ~isempty(extrafixpts2)
                for extrafixidx = 1:length(extrafixpts2)
                    % There is a trial end missing after each
                    % <extrafixpts2>. If there is an IEI > 8 s within 8
                    % events of the "extra" FixPtOn, then insert the the
                    % trial end at its beginning. Otherwise, simply
                    % terminate the trial right after the "extra" FixPtOn.
                    for evtidx = fixptonidx(extrafixpts2(extrafixidx))
                        IEIs = diff(events( evtidx + ...
                            (0 : min(size(events,1) - evtidx, 8)), 1 ));
                        longidx = find(IEIs > 8, 1);
                        if isempty(longidx)
                            newendTS(end+1, 1) = events(evtidx, 1) + 10;  %#ok<AGROW>
                        else
                            newendTS(end+1, 1) = ...
                                events(evtidx + longidx - 1, 1) + 10;  %#ok<AGROW>
                        end
                    end
                end
                warning('lfp_otherPreprocess:newstarts2', ...
                    'Inserting remedial trial end(s) for trial(s) %s at %s', ...
                    dg_thing2str(extrafixpts2), ...
                    dg_thing2str(newendTS/1e6));
            end
            events = sortrows([ events
                newstartTS repmat(MLTrialStart, size(newstartTS))
                newendTS repmat(MLTrialEnd, size(newendTS)) ]);
            trialstarts = find(events(:, 2) == MLTrialStart);
            trialends = find(events(:, 2) == MLTrialEnd);
        end
        % If there is a recording gap that occurs after FixPtOn, the only
        % way to know is by virtue of a ridiculously long trial.  If there
        % is no recording gap, the longest a trial can be is about 60
        % seconds from FixPtOn, which would be about 65 s including the
        % "baseline" period.  So here we check for both missing FixPtOn's
        % AND for trials that are longer than 90 s, and simply remove all
        % the events from MLTrialStart to MLTrialEnd.
        trialdurs = events(trialends, 1) - events(trialstarts, 1);
        trials2del = find(trialdurs > 90e6);
        evts2del = [];
        for trial2delidx = 1:length(trials2del)
            evts2del = [evts2del
                (trialstarts(trial2delidx) : trialends(trial2delidx))'
                ]; %#ok<AGROW> 
        end
        events(evts2del, :) = [];
        % That deletion invalidates these indices into <events>:
        clear fixptonidx trialstarts trialends
    otherwise 
        error('lfp_otherPreprocess:setup', 'Unknown value for lfp_SetupType: %s', lfp_SetupType);
end

