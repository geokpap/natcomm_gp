function h = dg_plottick(X, y, height, color, width)
%DG_PLOTTICK plots tick marks
%dg_plottick(X, Y, height, color, width)
% Plots tick marks of height <height> at the specified X(k), y coordinates.
% <X> must be a vector.  <y> must be scalar.  <color> can be an RGB triple
% or a Matlab predefined color name.  <h> is a column vector of handles to
% line graphics objects, one handle per tick.  <width> specifies the
% LineWidth used, is optional, and defaults to 0.5.

%$Rev: 265 $
%$Date: 2019-08-27 16:50:49 -0400 (Tue, 27 Aug 2019) $
%$Author: dgibson $

if nargin < 5
    width = 0.5;
end
y1 = y - height/2;
y2 = y + height/2;
h = zeros(length(X),1);
for k = 1:length(X)
    h(k) = line([X(k) X(k)], [y1 y2], 'Color', color, 'LineWidth', width);
end