function dg_plotPosAndAngle(x, y, theta)

%$Rev: 66 $
%$Date: 2010-08-04 18:00:51 -0400 (Wed, 04 Aug 2010) $
%$Author: dgibson $

linelength = 10;
x2 = x + linelength * cos(theta);
y2 = y - linelength * sin(theta);
hF = figure;
plot(x, y, 'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 18);
hold on;
plot([x'; x2'], [y'; y2'], 'Color', [0 .5 0]);
axis ij;

