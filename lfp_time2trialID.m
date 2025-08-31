function [trialID, trialnum] = lfp_time2trialID(time)
% Returns the Unique Trial ID and the trial number of the trial containing
% the absolute timestamp <time>.  This is done solely on the basis of the
% trial start events listed in lfp_TrialIndex, so <time> is considered to
% belong to the last trial that started at or before <time>.  If <time> is
% before the start of the first trial, then <trialID> = '' and trialnum =
% 0.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

trialsbefore = find(lfp_Events(lfp_TrialIndex(:,1), 1) <= time);
if isempty(trialsbefore)
    trialnum = 0;
    trialID = '';
else
    trialnum = trialsbefore(end);
    trialID = lfp_getTrialID(trialnum);
end