function lfp_selectByTrialID(trialIDs)
%lfp_selectByTrialID(trialIDs)
% Selects trials based on whether they are present in a list of Unique
% Trial IDs.  <trialIDs> is a cell string array containing Unique Trial
% IDs, or when there is only a single session loaded it may also contain
% string representations of only the original trial numbers
% (lfp_OrigTrialNums) without the session names.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

trialnums = [];
for k = 1:numel(trialIDs)
    trialnums(end+1) = lfp_getTrialNum(trialIDs{k});
end
lfp_SelectedTrials(:) = false;
lfp_SelectedTrials(trialnums) = true;
