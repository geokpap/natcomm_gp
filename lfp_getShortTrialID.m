function trialID = lfp_getShortTrialID(trialnum)
%LFP_GETSHORTTRIALID constructs a shorthand Trial ID from session number
%and trial number.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_SessionNames lfp_SessionOffsets lfp_Events lfp_TrialIndex ...
    lfp_OrigTrialNums;

sessionnums = find(lfp_SessionOffsets < ...
    lfp_Events(lfp_TrialIndex(trialnum, 1), 1));
trialID = ...
    [num2str(sessionnums(end)) '-' num2str(lfp_OrigTrialNums(trialnum))];

