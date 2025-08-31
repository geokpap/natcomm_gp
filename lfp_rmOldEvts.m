function lfp_rmOldEvts(evtIDs, name)
%lfp_rmOldEvts(evtIDs)
%lfp_rmOldEvts(evtIDs, name)
% Removes old events that have been superseded by newly added events in any
% trials that contain the new event.  Pays no attention to selection state
% of trials.
%INPUTS
% evtIDs: a two-element vector of event IDs where the first element is the
%   old ID that has been replaced in some trials by a new event with ID
%   given by evtIDs(2).
% name: a name to which to save the modified .evtsav file, not including
%   the .evtsav extension.  If empty or not given, then no file is saved.

%$Rev: 115 $
%$Date: 2010-03-31 15:38:29 -0400 (Wed, 31 Mar 2010) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 2
    name = '';
end

if evtIDs(1) == lfp_NominalTrialStart
    warning('lfp_rmOldEvts:trialstart', ...
        'You are removing lfp_NominalTrialStart events');
end
if evtIDs(1) == lfp_NominalTrialEnd
    warning('lfp_rmOldEvts:trialend', ...
        'You are removing lfp_NominalTrialEnd events');
end

evts2delete = zeros(0,1);
for trial = 1:length(lfp_SelectedTrials)
    trialevents = ...
        lfp_Events(lfp_TrialIndex(trial, 1):lfp_TrialIndex(trial,2), :);
    if ismember(evtIDs(2), trialevents(:,2))
        oldevtidx = find(ismember(trialevents(:,2), evtIDs(1)));
        if ~isempty(oldevtidx)
            evts2delete = [ evts2delete
                oldevtidx + lfp_TrialIndex(trial, 1) - 1 ];
        end
    end
end
lfp_Events(evts2delete,:) = [];
if ~isempty(name)
    lfp_save('preset', name, 'evt');
end
lfp_createEvents(@(a) a, 0, []);
