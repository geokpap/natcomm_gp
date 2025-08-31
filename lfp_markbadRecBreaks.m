function lfp_markbadRecBreaks
% Searches for trials containing recording breaks, and marks them bad by
% adding them to lfp_BadTrials.

%$Rev: 381 $
%$Date: 2016-07-28 18:17:52 -0400 (Thu, 28 Jul 2016) $
%$Author: dgibson $

global lfp_TrialRec lfp_RecSegments lfp_BadTrials

for trial = 1:size(lfp_TrialRec ,1)
    if find(lfp_RecSegments(:,1) == lfp_TrialRec(trial, 1)) ...
            ~= find(lfp_RecSegments(:,2) == lfp_TrialRec(trial, 2))
        lfp_BadTrials(end+1) = trial;
    end
end
