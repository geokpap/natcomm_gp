function trialID = lfp_getTrialID(trialnum)
%LFP_GETTRIALID constructs the unique Trial ID from session and trial
%number.
%trialID = lfp_getTrialID(trialnum)
%  <trialnum> must be less than or equal to length(lfp_OrigTrialNums).

%$Rev: 219 $
%$Date: 2011-04-20 17:14:48 -0400 (Wed, 20 Apr 2011) $
%$Author: dgibson $

global lfp_SessionNames lfp_SessionOffsets lfp_Events lfp_TrialIndex ...
    lfp_OrigTrialNums lfp_SelectedTrials;

if length(lfp_OrigTrialNums) ~= length(lfp_SelectedTrials)
    trialID = 'no trial ID';
    return
end
sessionnums = find(lfp_SessionOffsets < ...
    lfp_Events(lfp_TrialIndex(trialnum, 1), 1));
trialID = ...
    [lfp_SessionNames{sessionnums(end)} '-' num2str(lfp_OrigTrialNums(trialnum))];

