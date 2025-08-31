function hF = dg_plotBadVT(ts, v, badtimes)
%hF = dg_plotBadVT(ts, v, badtimes)
% Plot the value <v> as a function of time <ts> against a background that
% is colored red in all intervals specified by <badtimes>.  <ts> and
% <badtimes> can be in any units as long as they are consistent with each
% other.

%$Rev: 74 $
%$Date: 2010-08-20 17:17:42 -0400 (Fri, 20 Aug 2010) $
%$Author: dgibson $

hF = figure;
hA = axes('Parent', hF);
set(hA, 'NextPlot', 'add');
maxv = max(v);
for k = 1:size(badtimes,1)
    rectangle('Parent', hA, 'Position', ...
        [badtimes(k,1), 0, badtimes(k,2) - badtimes(k,1), maxv], ...
        'FaceColor', 'r', 'EdgeColor', 'none');
end
hL = plot(hA, ts, v);
xlabel('Time');
ylabel('Value');

