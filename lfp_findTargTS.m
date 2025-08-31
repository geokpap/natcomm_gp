function ts = lfp_findTargTS(evtidx, oldtrial)
% Helper for trial re-marking in lfp_handleNewTrialTarget; returns target
% timestamp, or [] if target not found.

%$Rev: 318 $
%$Date: 2014-02-08 22:12:19 -0500 (Sat, 08 Feb 2014) $
%$Author: dgibson $

global lfp_Events lfp_NewTrialTarget
ts = [];
targidx = find(ismember(lfp_Events(evtidx,2), lfp_NewTrialTarget{1}));
if isempty(targidx)
    warning('lfp_findTargTS:noTarg', ...
        'Old trial #%d contains no lfp_NewTrialTarget event.', ...
        oldtrial );
else
    ts = lfp_Events(...
        evtidx(1) + targidx(1) - 1, 1 );
end
