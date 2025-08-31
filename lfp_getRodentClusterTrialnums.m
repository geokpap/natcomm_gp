function trialnums = lfp_getRodentClusterTrialnums(trialdata)
%trialnums = lfp_getRodentClusterTrialnums(trialdata)
%   Applies the procedure formerly known as "brute force match" to divine
%   the correct correspondence between trials in a Rodent Cluster file and
%   trials in lfp_TrialIndex.  The general strategy is to find two events
%   that are separated by a variable time interval in each putative
%   matching pair of trials, and to accept the pair as a match only if the
%   time interval is the same to within 1.5e-4 sec (1.5x the resolution of
%   the Rodent Cluster format).  We explicitly ignore the trial numbers in
%   the <trialdata.header> field because they do not generally agree with
%   the trial numbers in lfp_TrialIndex.  In this implementation, for the
%   starting event we use event 11 (Gate) if it is present, 10 (Warning
%   Cue) otherwise.  For the ending event we use 14 (Turn Onset) if it is
%   present, or 23 (mid-T in Jianbin data) if present, or [31 38 20 21 22]
%   (Tone On or Tactile Cue On) otherwise.  All of these alternatives
%   "should" handle mouse data, Jianbin data, new (video tracker assisted)
%   data, and old (100% photobeam) data, in the face of known cases where
%   photobeam or video tracker events were missing.  <trialnums> is a list
%   of the trialnums that match the trials in <trialdata>, containing one
%   element for each trial in <trialdata> except for trailing empty trials.
%   If the first trial in <trialdata> does not match, the first element is
%   NaN.  There may be elements that are NaN interspersed in the middle as
%   well.  When difficulties are encountered matching trials between
%   lfp_lib and the T file, a graphical representation of the match matrix
%   is shown.

% MODIFICATION NOTES:
% The code that finds evtidx1 and evtidx2 in the lfp_lib events must match
% the code that finds startEvents and endEvents in the Rodent Cluster file.
% Specifically, note that the list of event IDs [31 38 20 21 22] is NOT a
% list of alternative end events, it is just a list of alternative Stimulus
% On events that constitutes the last choice for which event to use.  Also,
% the logic of which events to accept must be identical as regards
% duplicate events, missing events, incorrectly timestamped events.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

trialnums = []; % later, this is made to be a row vector
numtrials = size(lfp_TrialIndex, 1);
numrodenttrials = length(trialdata);

% eventses has all events for each trial, one trial per row, as represented
% in Rodent Cluster file:
eventses = cell(size(trialdata));
[eventses{:}] = deal(trialdata.events);
eventses = cell2mat(eventses');

% Check for empty trials:
emptytrials = find(any(eventses,2) == 0);
if ~isempty(emptytrials)
    if isequal(emptytrials, (emptytrials(1):size(eventses,1))')
        % simply delete trailing empty trials:
        eventses(emptytrials,:) = [];
        numrodenttrials = size(eventses,1);
    else
        error('lfp_getRodentClusterTrialnums:baddata', ...
            'The Rodent Cluster data include empty trials that are not at the end of the session.' );
    end
end

% Do feasibility check:
old_selectedtrials = lfp_SelectedTrials;
lfp_selectByRule('HasEvent([11 10])&&HasEvent([14 23 31 38 20 21 22])');
if ~all(lfp_SelectedTrials)
    warning('lfp_getRodentClusterTrialnums:evt1', ...
        'Cannot find start and end events in lfp_Events for all trials.' );
end
lfp_SelectedTrials = old_selectedtrials;
if ~all((eventses(:,11)| eventses(:,10)) & ...
        (eventses(:,14) | eventses(:,23) | eventses(:,31) | eventses(:,38) ))
    warning('lfp_getRodentClusterTrialnums:evt2', ...
        'Cannot find start and end events in Rodent Cluster file for all trials.' );
end

% Calculate intervals from lfp_Events (i.e. from Events.Nev / Events.dat)
% SEE HEADER COMMENTS!
lfpstartidx = zeros(numtrials,1);
lfpendidx = zeros(numtrials,1);
for trial = 1:numtrials
    dumpevts = false;
    % find start event
    evtidx1 = find(ismember(lfp_Events( ...
          lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), 2 ),11  )) ...
        + lfp_TrialIndex(trial,1) - 1;
    if  isempty(evtidx1)
         evtidx1 = find(ismember(lfp_Events( ...
          lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), 2 ),10)) ...
        + lfp_TrialIndex(trial,1) - 1;
    end
    if length(evtidx1) > 1
        % Eliminate multiple triggers; the first trig is probably more
        % accurate, but we use the last for consistency with Delphi
        % analysis programs.
        warning('lfp_getRodentClusterTrialnums:multitrig1', ...
            'Multiple triggers of starting event in lfp_lib trial %d', trial );
        evtidx1 = evtidx1(end);
        dumpevts = true;
    end
    
    % find end event; if it results in a negative interval, proceed as if
    % the event did not exist.
    evtidx2 = find(ismember(lfp_Events( ...
        lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), 2 ), 14 )) ...
        + lfp_TrialIndex(trial,1) - 1;
    if ~isempty(evtidx2) && ~isempty(evtidx1) && evtidx2(end) < evtidx1
        warning('lfp_getRodentClusterTrialnums:negint1', ...
            'Negative interval in lfp_lib trial %d; skipping', trial );
        evtidx2 = [];
        dumpevts = true;
    end
    if  isempty(evtidx2)
        evtidx2 = find(ismember(lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), 2 ), 23 )) ...
            + lfp_TrialIndex(trial,1) - 1;
    end
    if ~isempty(evtidx2) && ~isempty(evtidx1) && evtidx2(end) < evtidx1
        warning('lfp_getRodentClusterTrialnums:negint2', ...
            'Negative interval in lfp_lib trial %d; skipping', trial );
        evtidx2 = [];
        dumpevts = true;
    end
    if  isempty(evtidx2)
        evtidx2 = find(ismember(lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), 2 ),[31 38 20 21 22] )) ...
            + lfp_TrialIndex(trial,1) - 1;
    end
    if ~isempty(evtidx2) && ~isempty(evtidx1) && evtidx2(end) < evtidx1
        warning('lfp_getRodentClusterTrialnums:negint3', ...
            'Negative interval in lfp_lib trial %d; skipping', trial );
        evtidx2 = [];
        dumpevts = true;
    end
    if length(evtidx2) > 1
        warning('lfp_getRodentClusterTrialnums:multitrig2', ...
            'Multiple triggers of ending event in lfp_lib trial %d', trial );
        evtidx2 = evtidx2(end);
        dumpevts = true;
    end
    
    % compute interval
    if isempty(evtidx1) || isempty(evtidx2)
        interval(trial,1) = NaN;
        warning('lfp_getRodentClusterTrialnums:interval', ...
            'Could not find comparison interval in lfp_lib trial %d', trial );
        dumpevts = true;
    else
        lfpstartidx(trial) = evtidx1;
        lfpendidx(trial) = evtidx2;
        ts1 = lfp_Events(evtidx1, 1);
        ts2 = lfp_Events(evtidx2, 1);
        interval(trial,1) = ts2 - ts1;
    end
    
    % show events details
    if dumpevts
        detailstr = 'TS           ID';
        for evtidx = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2)
            detailstr = sprintf('%s\n%11.6f %3.0f', ...
                detailstr, lfp_Events(evtidx,1), lfp_Events(evtidx,2));
        end
        warning('lfp_getRodentClusterTrialnums:evts', ...
            'The events for lfp_lib trial %d were:\n%s', trial, detailstr );
    end
end

% Calculate intervals <rodentinterval> from Rodent Cluster file (find that
% special interval, and reject negative intervals as if event were
% missing).  I am assuming that the Delphi programs that generate the
% T-files and VT events files simply process each event as they arise,
% meaning that for duplicate events we end up with the timestamp of the
% LAST copy. SEE HEADER COMMENTS!
startEvents = eventses(:,11);
rodentid1 = 11*ones(size(startEvents));
startEvents(startEvents == 0) = eventses(startEvents == 0,10);
rodentid1(startEvents == 0) = 10;
endEvents = eventses(:,14);
rodentid2 = 14*ones(size(endEvents));
rodentinterval = 1e-4 * (endEvents - startEvents);
for endevtid = [23 31 38 20 21 22]
    reject = endEvents == 0 | rodentinterval < 0;
    endEvents(reject) = eventses(reject,endevtid);
    rodentid2(reject) = endevtid;
    rodentinterval = 1e-4 * (endEvents - startEvents);
end
for trial = 1:length(startEvents)
    if startEvents(trial) == 0
        rodentid1(trial) = 0;
    end
    if endEvents(trial) == 0
        rodentid2(trial) = 0;
    end
end

% Create pairings (numrodenttrials is number of trials in the T file while
% numtrials is the number of trials in the LFP trial index gotten from the
% Events.nev or events.dat).  <matchmatrix> is in rodenttrials X trialnums
% form.
diffmatrix = repmat(rodentinterval, 1, size(interval,1)) ...
    - repmat(interval', size(rodentinterval,1), 1);
matchmatrix = abs(diffmatrix) <= 1.5e-4;
[k, runstart, runsum] = dg_findDiagRuns(matchmatrix);
trialnums = NaN(1, numrodenttrials);
showmatchmatrix = false;
prevrunsum = 1;
for row = 1:length(runsum)
    if runsum(row) < 1
        break
    end
    rodentstarttrial = runstart(row) - min(0, k(row));
    lfpstarttrial = runstart(row) + max(0, k(row));
    thisrunlen = runsum(row);
    previouslyassigned = ~isnan(trialnums( ...
        rodentstarttrial : (rodentstarttrial + thisrunlen - 1) ));
    if any(previouslyassigned)
        % If it's just the first or last in the run, use the others in the
        % run; otherwise, ignore the run.
        showmatchmatrix = true;
        previouslyassignedidx = find(previouslyassigned);
        warning('lfp_getRodentClusterTrialnums:diagruns3', ...
            'Some of rodent cluster trials %s were already matched in a longer run', ...
            dg_canonicalSeries(previouslyassignedidx + rodentstarttrial - 1) );
        if sum(previouslyassigned) == 1
            if previouslyassignedidx == 1
                warning('lfp_getRodentClusterTrialnums:diagruns3a', ...
                    'Skipping first match in run' );
                rodentstarttrial = rodentstarttrial + 1;
                thisrunlen = thisrunlen - 1;
            elseif previouslyassignedidx == thisrunlen
                warning('lfp_getRodentClusterTrialnums:diagruns3b', ...
                    'Skipping last match in run' );
                thisrunlen = thisrunlen - 1;
            else
                warning('lfp_getRodentClusterTrialnums:diagruns3c', ...
                    'Skipping run: trialnums %d:%d, rodent cluster trials %d:%d', ...
                    k(row)+1, k(row)+thisrunlen, ...
                    rodentstarttrial, rodentstarttrial + thisrunlen - 1 );
                continue
            end
        end
    end
    trialnums(rodentstarttrial : (rodentstarttrial + thisrunlen - 1)) = ...
        lfpstarttrial : (lfpstarttrial + thisrunlen - 1);
    if thisrunlen < prevrunsum/3
        showmatchmatrix = true;
        warning('lfp_getRodentClusterTrialnums:diagruns1', ...
            'Abrupt drop in length of match starting trialnum %d, rodent cluster trial %d', ...
            trialnums(rodentstarttrial), rodentstarttrial );
    end
    if thisrunlen < size(lfp_TrialIndex,1)/10
        showmatchmatrix = true;
        warning('lfp_getRodentClusterTrialnums:diagruns2', ...
            'Run matches less than 10%% of lfp_lib trials starting trialnum %d, rodent cluster trial %d', ...
            trialnums(rodentstarttrial), rodentstarttrial );
    end
end

% Check for multiple rodent cluster trials matched to a single trialnum
matchedrodenttrials = ~isnan(trialnums);
if sum(matchedrodenttrials) ~= ...
        length(unique(trialnums(matchedrodenttrials)))
    matchcounts = zeros(size(lfp_SelectedTrials));
    for trialnum = trialnums(matchedrodenttrials)
        matchcounts(trialnum) = matchcounts(trialnum) + 1;
    end
    detailstr = '';
    for trialnum = find(matchcounts>1)
        detailstr = sprintf('%s\ntrialnum %d matched rodent trials %s', ...
            detailstr, trialnum, dg_canonicalSeries(find(trialnums == trialnum)) );
    end
    warning('lfp_getRodentClusterTrialnums:multimatch', ...
        'Some trialnums matched more than one rodent cluster trial:%s', ...
        detailstr );
    showmatchmatrix = true;
end

if any(isnan(trialnums))
    showmatchmatrix = true;
    detailstr = '';
    badRCtrials = find(isnan(trialnums));
    % At this point, individual trials that do match uniquely but that are
    % not part of a diagonal run are still marked NaN, so they should be
    % matched if possible.
    matchedBad = false(size(badRCtrials));
    for badRCidx = 1:length(badRCtrials)
        RCtrial = badRCtrials(badRCidx);
        matchidx = find(matchmatrix(RCtrial,:));
        if length(matchidx) == 1
            trialnums(RCtrial) = matchidx;
            matchedBad(badRCidx) = true;
        end
    end
    % Ambiguously matched RC trials can now potentially be disambiguated by
    % choosing the match that fits with the surrounding matches:
    for badRCidx = 1:length(badRCtrials)
        RCtrial = badRCtrials(badRCidx);
        matchidx = find(matchmatrix(RCtrial,:));
        if length(matchidx) > 1
            % all comparisons to NaN are false except for ~=, so this formula
            % automatically returns false when the previous or following trial
            % is unmatched:
            isqualified = matchidx > trialnums(RCtrial-1) ...
                & matchidx < trialnums(RCtrial+1);
            if sum(isqualified) == 1
                trialnums(RCtrial) = matchidx(isqualified);
                matchedBad(badRCidx) = true;
            end
        end
    end
    if any(~matchedBad)
        % Gather detail info on any remaining unmatched RC trials
        for badRCidx = 1:length(badRCtrials)
            RCtrial = badRCtrials(badRCidx);
            if isnan(trialnums(RCtrial))
                if rodentid1(RCtrial)
                    detailstr = sprintf('%s\ntrial %d comparison interval start evt ID=%d t=%.6f sec', ...
                        detailstr, RCtrial, rodentid1(RCtrial), ...
                        eventses(RCtrial, rodentid1(RCtrial))*1e-6 );
                else
                    detailstr = sprintf('%s\ntrial %d no comparison interval start evt found', ...
                        detailstr, RCtrial );
                end
                if rodentid2(RCtrial)
                    detailstr = sprintf('%s\n\tcomparison interval end evt ID=%d t=%.6f sec', ...
                        detailstr, rodentid2(RCtrial), ...
                        eventses(RCtrial, rodentid2(RCtrial))*1e-6 );
                else
                    detailstr = sprintf('%s\n\tno comparison interval end evt found', ...
                        detailstr );
                end
            end
        end
        lfp_log(sprintf( ...
            'Rodent Cluster Trial(s) %s could not be matched.%s', ...
            dg_canonicalSeries(badRCtrials(~matchedBad)), detailstr ));
        warning('lfp_getRodentClusterTrialnums:badRCtrials', ...
            'Rodent Cluster Trial(s) %s could not be matched.', ...
            dg_canonicalSeries(badRCtrials(~matchedBad)) );
    end
end

badtrials = setdiff(1:size(lfp_TrialIndex,1), trialnums);
if ~isempty(badtrials)
    message = sprintf(...
        'These lfp_lib trialnums have no match in the Rodent Cluster file:\n%s', ...
        dg_canonicalSeries(badtrials) );
    message = sprintf('%s\nAdding to lfp_BadTrials, whose previous value was %s', ...
        message, dg_canonicalSeries(lfp_BadTrials) );
    warning('lfp_getRodentClusterTrialnums:badtrials', message);
    for k = 1:length(badtrials)
        trial = badtrials(k);
        if lfpstartidx(trial)
            message = sprintf('%s\ntrial %d comparison interval start evt ID=%d t=%.6f sec', ...
                message, trial, lfp_Events(lfpstartidx(trial),2), ...
                lfp_Events(lfpstartidx(trial),1) );
        else
            message = sprintf('%s\ntrial %d no comparison interval start evt found', ...
                message, trial );
        end
        if lfpendidx(trial)
            message = sprintf('%s\ntrial %d comparison interval end evt ID=%d t=%.6f sec', ...
                message, trial, lfp_Events(lfpendidx(trial),2), ...
                lfp_Events(lfpendidx(trial),1) );
        else
            message = sprintf('%s\ntrial %d no comparison interval end evt found', ...
                message, trial );
        end
    end
    lfp_log(message);
    lfp_BadTrials = union(lfp_BadTrials, badtrials);
end

if showmatchmatrix
    alreadyshown = false;
    titlestr = sprintf('matchmatrix for %s', lfp_SessionNames{end});
    hF = findobj(0, 'Type', 'figure');
    for k = 1:length(hF)
        hA = findobj(hF(k), 'Type', 'axes');
        for j = 1:length(hA)
            oldtitlestr = get(get(hA(j), 'Title'), 'String');
            if isequal(oldtitlestr, titlestr)
                alreadyshown = true;
                break
            end
        end
    end
    if ~alreadyshown
        figure;
        imagesc(matchmatrix);
        caxis([-1 2]);
        grid on;
        title(titlestr, 'Interpreter', 'none');
        ylabel('Rodent Cluster trial', 'Interpreter', 'none');
        xlabel('lfp_lib trialnum', 'Interpreter', 'none');
    end
end

