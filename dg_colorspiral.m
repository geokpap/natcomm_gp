function cmap = dg_colorspiral(N, hstep, slim, vlim)
% Produces a colormap of N colors that traverses a spiral (sort of) in the
% HSV color space.  Hue varies by successive increments of <hstep>, but is
% interpreted mod 1.
% Saturation alternates between slim(1) and slim(2) with each wrap around
% the hue circle. Value ramps from svlim(1) to svlim(2).

% To get hue values to repeat every kth wrap, use n*hstep - k = hstep, i.e.
% hstep = k/(n-1) where n is any integer s.t. (n-1) has no factors in
% common with k.

% Practically speaking, it's hard to get more than about 100
% distinguishable colors or 64 readily distinguishable colors.  However,
% this is a significant improvement over jet(64), whose colors are NOT
% readily distinguishable.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

myhues = (0:(N-1)) * hstep;
mysats = slim(mod(floor(myhues),2) + 1);
myvals = vlim(1) + (0:(N-1)) * (vlim(2) - vlim(1)) / (N-1);

cmap = zeros(N, 3);
for cidx = 1:N
    cmap(cidx,:) = hsv2rgb( ...
        [mod(myhues(cidx), 1) mysats(cidx) myvals(cidx)] );
end
