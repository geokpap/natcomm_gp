function targID = lfp_xy2targID(x, y, targets)
%lfp_xy2targID converts cartesian eye coordinates to a target ID.
%targID = xy2targID(x, y, targets)
%   <x>, <y> are as in the <fixations> table returned by lfp_EyeTabulation.
%   <targets> is as supplied to lfp_TargSeq.  If (x,y) is not in any
%   target, returns 0.  It is assumed that the <targets> are strictly
%   non-overlapping and therefore there can be no ambiguities.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

targID = 0;
for targidx = 1:size(targets,1)
    if sqrt((x - targets(targidx, 1))^2 + (y - targets(targidx, 2))^2) ...
            < targets(targidx, 3)
        targID = targidx;
        break
    end
end
