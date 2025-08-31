function lfp_fakeEvents

% Fake up whole session as one trial, timestamps in seconds.
% Must be called after lfp_TimeStamps has final values in seconds.
% This is rodent-specific.

%$Rev: 363 $
%$Date: 2015-08-27 19:50:58 -0400 (Thu, 27 Aug 2015) $
%$Author: dgibson $

lfp_declareGlobals;

lfp_Events(1,:) = [ lfp_TimeStamps(1) + 1e-6, lfp_NominalTrialStart ];
lfp_Events(2,:) = [ lfp_TimeStamps(end) + ...
        (lfp_SamplesPerFrame-1)*lfp_SamplePeriod - 10e-6, lfp_NominalTrialEnd ];
lfp_Events(3,:) = [ lfp_TimeStamps(end) + ...
        (lfp_SamplesPerFrame-1)*lfp_SamplePeriod, lfp_GoodTrialEvent(1) ];
