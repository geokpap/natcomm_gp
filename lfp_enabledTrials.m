function enabledtrials = lfp_enabledTrials(trials)
%LFP_ENABLEDTRIALS combines lfp_SelectedTrials and lfp_BadTrials
%enabledtrials = lfp_enabledTrials
%enabledtrials = lfp_enabledTrials(trials)
%  Accepts a list of trials and returns the same list minus any trials that
%  are not enabled for general display and processing purposes (e.g.
%  computing averages).  If called without arguments, returns all enabled
%  trials as a row vector.
% Examples (both calls return the same value):
%   trials = lfp_enabledTrials(1:length(lfp_SelectedTrials))
%   trials = lfp_enabledTrials

%  For a trial to be enabled, two conditions must hold:  the value of
%  lfp_SelectedTrials must be true, and the trial must not appear in
%  lfp_BadTrials.

%$Rev: 399 $
%$Date: 2019-06-20 13:34:13 -0400 (Thu, 20 Jun 2019) $
%$Author: dgibson $

global lfp_SelectedTrials lfp_BadTrials

if nargin == 0
    trials = 1:length(lfp_SelectedTrials);
end
enabledtrials = trials(logical(lfp_SelectedTrials(trials)));
enabledtrials = enabledtrials(~ismember(enabledtrials, lfp_BadTrials));