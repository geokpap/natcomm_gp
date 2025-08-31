function [trialindex, events] = lfp_handleNewTrialTarget(trialindex, ...
    events, numTP)
% The name says it all.  Call it with the [lfp_NominalTrialStart
% lfp_NominalTrialEnd] values as specified in the setup, and it will change
% them to [lfp_NewTrialStart lfp_NewTrialEnd].

%$Rev: 351 $
%$Date: 2015-05-28 12:32:51 -0400 (Thu, 28 May 2015) $
%$Author: dgibson $

global lfp_NewTrialTarget lfp_NewTrialStart lfp_NewTrialEnd ...
    lfp_NominalTrialStart lfp_NominalTrialEnd ...
    lfp_SelectedTrials lfp_SelectionRule

insertingNewTrialStarts = false;

% Insert new trial start and end events if setup is configured for that
% and it has not already been done:
if any(ismember(events(:,2), lfp_NewTrialStart))
    warning('lfp_read2:newevents1', ...
        ['This events file already contains new trial starts;' ...
        ' no new events\nwill be inserted.'] );
else
    warning('lfp_read2:newevents2', ...
        'Inserting new trial starts/ends.' );
    insertingNewTrialStarts = true;
    targTS = [];
    % first "old trial" is from start of data through first old trial
    % end:
    trialevtidx = 1:trialindex(1,2);
    targTS = [targTS lfp_findTargTS(trialevtidx, 1)];
    % subsequent "old trials" are from (not including) one old trial
    % end through the next:
    for oldtrial = 2:size(trialindex,1)
        trialevtidx = ...
            trialindex(oldtrial-1,2) + 1 : trialindex(oldtrial,2);
        targTS = [targTS lfp_findTargTS(trialevtidx, oldtrial)]; %#ok<AGROW>
    end
    newevents = zeros(2*size(targTS,1), 2);
    % One pair of new events for each target event:
    newevents(1:length(targTS),1) = targTS + lfp_NewTrialTarget{2};
    newevents(1:length(targTS),2) = lfp_NewTrialEnd(1);
    newevents((length(targTS)+1):(2*length(targTS)),1) = ...
        targTS + lfp_NewTrialTarget{2} + 1e-6;
    newevents((length(targTS)+1):(2*length(targTS)),2) = ...
        lfp_NewTrialStart(1);
    % The first newevent is an undesired lfp_NewTrialEnd; we also must
    % add an extra lfp_NewTrialEnd to end the last trial:
    events = [ events
        newevents(2:end,:)
        [events(end,1) + 1e-6, lfp_NewTrialEnd(1)] ];
    events = sortrows(events);
    clear newevents;
end
% update lfp_NominalTrialStart and lfp_NominalTrialEnd, and create new
% Trial Index:
lfp_NominalTrialStart = lfp_NewTrialStart;
lfp_NominalTrialEnd = lfp_NewTrialEnd;
if exist('taskStopID', 'var') && any(events(:,2)==taskStopID)
    s1 = warning('off', 'lfp_createTrialIndex:extraStarts');
    s2 = warning('off', 'lfp_createTrialIndex:extraEnds');
end
% This is the final value of trialindex:
trialindex = lfp_createTrialIndex(events, numTP);
if exist('taskStopID', 'var') && any(events(:,2)==taskStopID)
    warning(s1.state, 'lfp_createTrialIndex:extraStarts');
    warning(s2.state, 'lfp_createTrialIndex:extraEnds');
end
% if StopTask is defined and any trial contains one, then the first
% StopTask in that trial should be used as the trial end event:
if exist('StopTask', 'var')
    oldtrialselection = lfp_SelectedTrials;
    oldselectionrule = lfp_SelectionRule;
    lfp_selectByRule(sprintf('HasEvent(%d)', StopTask));
    trials = lfp_enabledTrials(1:size(trialindex,1));
    for trial = trials
        stoptaskidx = find(events(...
            trialindex(trial,1):trialindex(trial,2), 2 ) ...
            == StopTask );
        trialindex(trial,2) = ...
            stoptaskidx(1) + trialindex(trial,1) - 1;
    end
    lfp_SelectedTrials = oldtrialselection;
    lfp_SelectionRule = oldselectionrule;
end

% Note that the first ID in lfp_NominalTrialEnd is used as a replacement
% for taskStopID, so if lfp_NominalTrialEnd includes taskStopID then
% taskStopID should NOT be first!
if insertingNewTrialStarts && exist('taskStopID', 'var') ...
        && any(events(:,2)==taskStopID)
    taskstops = find(events(:,2)==taskStopID);
    if taskStopID==255
        % remove searcheventcount 255s from <taskstops>.  The -30, -5, and
        % +5 added to curr_idx were empirically determined by Theresa to
        % work on lots of sessions, but it could still be logically
        % possible to miss some events or raise some false alarms.
        evt255check = false(length(taskstops),1);
        if isempty(taskstops)
            if ~silentflag
                fprintf('%s All Clear (no taskstops)\n', sess_name);
            end
        else
            for i = 1:length(taskstops)
                curr_idx = taskstops(i);
                % previous count
                exist_countb4 = ismember(254, ...
                    events(max(1,curr_idx-30):curr_idx,2));
                % previous exit
                exist_32b4 = ismember(32, ...
                    events(max(1,curr_idx-30):curr_idx,2));
                exist_33to57b4 = ...
                    sum(ismember(33:57, ...
                    events(max(1,curr_idx-5):curr_idx,2))) > 0;
                exist_32after = ismember(32, events(curr_idx : ...
                    min(size(events,1),curr_idx+5),2 ));
                if exist_countb4 && exist_32b4 && exist_33to57b4 ...
                        && exist_32after
                    evt255check(i) = 1;
                    lfp_log(sprintf( ...
                        'searcheventcount=255 at %.6f s', ...
                        events(taskstops(i), 1) ));
                else
                    if taskstops(i) < size(events,1) ...
                            && events(taskstops(i)+1,1) ...
                            < events(taskstops(i),1)+1
                        % Idiot check - actually, Theresa says that all of
                        % the cases where this is triggered were
                        % missed searcheventcount 255s, so they should be
                        % treated as such.
                        evt255check(i) = 1;
                        warning('lfp_read2:idiot', ...
                            'Treating taskstop at %.6f s as searcheventcount=255', ...
                            events(taskstops(i), 1) );
                        lfp_log(sprintf( ...
                            'Treating taskstop at %.6f s as searcheventcount=255; exist_countb4=%d exist_32b4=%d exist_33to57b4=%d exist_32after=%d', ...
                            events(taskstops(i), 1), ...
                            exist_countb4, exist_32b4, exist_33to57b4, ...
                            exist_32after ));
                    else
                        lfp_log(sprintf( ...
                            'taskStopID at %.6f s; exist_countb4=%d exist_32b4=%d exist_33to57b4=%d exist_32after=%d', ...
                            events(taskstops(i), 1), ...
                            exist_countb4, exist_32b4, exist_33to57b4, ...
                            exist_32after ));
                        eventdisp = '';
                        for k = -40:10
                            if curr_idx + k > 0
                                eventdisp = sprintf('%s%20.6f%8.0f\n', ...
                                    eventdisp, events(curr_idx + k, 1), ...
                                    events(curr_idx + k, 2) );
                            end
                        end
                        lfp_log(sprintf( ...
                            'Events -40:10 re: taskStopID at %.6f s:\n%s', ...
                            events(taskstops(i), 1), eventdisp ));
                    end
                end
            end
        end
        msg = sprintf( ...
            '%s: %d of %d taskstops due to searcheventcount', ...
            sessionname, sum(evt255check), ...
            length(taskstops));
        lfp_log(msg);
        warning('lfp_read2:taskstops', '%s', msg);
        taskstops(evt255check) = [];
    end
    % Process the actual taskstops:
    events2delete = [];
    for k = 1:length(taskstops)
        evtix = taskstops(k);
        while evtix <= size(events,1) && ...
                ~ismember(events(evtix, 2), lfp_NominalTrialEnd)
            evtix = evtix + 1;
        end
        if evtix <= size(events,1)
            events2delete(end+1) = evtix; %#ok<AGROW>
        end
        events(taskstops(k), 2) = lfp_NominalTrialEnd(1);
    end
    events(events2delete,:) = [];
    trialindex = lfp_createTrialIndex(events);
end
