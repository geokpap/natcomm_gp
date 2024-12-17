function result = lfp_trialHasEvent(trial, event)
%LFP_TRIALHASEVENT helps select trials based on specified conditions.
%result = lfp_trialHasEvent(trial, event)
% Returns true if trial number <trial> has at least one event with event ID
% <event>.  "Has an event" means that the event occurs after the
% lfp_NominalTrialStart event and before the lfp_NominalTrialEnd event for
% that trial.  <event> may be a list, in which case <true> is returned if
% the trial has any of the events in the list; thus lfp_AlignmentRef is a
% legitimate value for the second argument.

%$Rev: 389 $
%$Date: 2017-06-12 15:07:03 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

global lfp_TrialIndex lfp_Events;

startindex = lfp_TrialIndex(trial, 1);
endindex = lfp_TrialIndex(trial, 2);
foundevents = find(ismember(lfp_Events(startindex:endindex, 2), event), 1);
if isempty(foundevents)
    result = false;
else
    result = true;
end