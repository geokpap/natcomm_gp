function [trialstartindex, trialendindex] = ...
    lfp_monkeyCreateTrialIndex(events)
%[trialstartindex, trialendindex] = lfp_monkeyCreateTrialIndex(events)
% <trialstartindex>, <trialendindex> are row vectors.
% Does not change any global values.

%$Rev: 318 $
%$Date: 2014-02-08 22:12:19 -0500 (Sat, 08 Feb 2014) $
%$Author: dgibson $

global lfp_NominalTrialStart lfp_NominalTrialEnd

lfp_getEvtIDs;

% This is specifically based on lfp_NominalTrialStart = BDOutputStop and
% lfp_NominalTrialEnd = BDOutputStart.  However, in RAW data, there is
% usually no BDOutputStop before the first trial, so we use the first event
% as the start of the first trial.  Also, there is usually a BDOutputStop
% after the last trial which we must ignore, but that will be done
% automatically in lfp_createTrialIndex when we dg_zip.
%
% The "exceptional" cases would be:
% 1. the very rare case where recording started just in time to catch the
%   BDOutputStop from the preceding (unrecorded) trial.  In this case, the
%   first BDOutputStop IS the first event in the trial, so there is no need
%   to insert another trialstartindex entry.
% 2. the common case where the events come from a perfectly cleaned-up
%   EVTSAV fragment file that begins with a lfp_NominalTrialStart and ends
%   with a BDOutputStop.  The terminating BDOutputStop becomes an extra
%   trialstart which is removed automatically in lfp_createTrialIndex when
%   we dg_zip as for RAW data.  This includes the case where the fragment
%   includes the end of the session.

% Modification of 4/6/2006
% To make matters even messier, we now have the possibility that Joey's
% 'monkey' data is set up like Theresa's as a result of using
% lfp_NewTrialTarget, but we can't just turn Joey from a 'monkey' into a
% 'theresa' because that will screw up other things.  So we now have to
% explicitly test to see if lfp_NominalTrialStart = BDOutputStop.


if isequal(lfp_NominalTrialStart, BDOutputStop)
    % treat like a 'monkey':
    bdoutputstops = find(events(:,2) == BDOutputStop);
    if bdoutputstops(1) > 1
        % The "usual" RAW data case:
        trialstartindex = [ 1 bdoutputstops' ];
    else
        % Case 1 or 2:
        trialstartindex = bdoutputstops';
    end
    trialendindex = find(ismember(events(:,2), lfp_NominalTrialEnd))';
else
    % treat like a 'theresa':
    trialstartindex = find(ismember(events(:,2), lfp_NominalTrialStart))';
    trialendindex = find(ismember(events(:,2), lfp_NominalTrialEnd))';
end
