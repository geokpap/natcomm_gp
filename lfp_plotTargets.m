function lfp_plotTargets(targets)
%lfp_plotTargets(targets) plots lfp_TargSeq-style targets.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

figure;
for idx = 1:size(targets,1)
    r = targets(idx,3);
    if r > 0
        rectangle( ...
            'Position', ...
            [targets(idx,1)-r, targets(idx,2)-r, 2*r, 2*r], ...
            'Curvature', [1,1]);
        text(targets(idx,1), targets(idx,2), num2str(idx), 'HorizontalAlignment','center');
    end
end
daspect([1,1,1]);
