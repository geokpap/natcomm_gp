function [hF, r, p, c] = dg_scatterplot(x, y, titlestr)
% <x> and <Y> must be column vectors

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

c = polyfit(x, y, 1);
[r, p] = corr(x, y);
titlestr = sprintf('%s\nr = %.2g  p = %1.1g', titlestr, r, p);
hF = figure;
plot(x, y, '.');
hA = get(hF, 'CurrentAxes');
hold on;
plot([min(x) max(x)], [c(1)*min(x) + c(2) c(1)*max(x) + c(2)]);
title(hA, titlestr);
