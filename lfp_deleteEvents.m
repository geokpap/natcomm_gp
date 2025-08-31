function lfp_deleteEvents(eventIDs)
%lfp_deleteEvents(eventIDs)
%   Deletes all events from memory whose event IDs are listed in <eventIDs>.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

lfp_Events(find(ismember(lfp_Events(:,2), eventIDs)), :) = [];
disp('Recalculating lfp_TrialIndex...');
lfp_TrialIndex = lfp_createTrialIndex(lfp_Events);
if ~isempty(lfp_SamplePeriod)
    [startSampleIndex, endSampleIndex] = ...
        lfp_addTIcols3n4(lfp_TrialIndex(:,1), lfp_TrialIndex(:,2), lfp_Events);
    lfp_TrialIndex = [lfp_TrialIndex startSampleIndex endSampleIndex];
end
lfp_log(sprintf('Deleted events %s', dg_thing2str(eventIDs)));
