function cmap = dg_colorscale(c1, c2, brightlim, n)
%cmap = dg_colorscale(c1, c2, brightlim, n)
%cmap = dg_colorscale(c1, c2, brightlim)
%cmap = dg_colorscale(c1, c2)
%cmap = dg_colorscale
% Returns a colormap with <n> elements that grade smoothly from <c1> to
% <c2> (which are both RGB triples) and that are also scaled with a gray
% scale that ranges from <brightlim(1)> to <brightlim(2)>.
%DEFAULTS
% n: 64
% brightlim: [0.5 1]
% c1: [0 0 1]
% c2: [1 0 0]

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 4
    n = 64;
end
if nargin < 3
    brightlim = [0.5 1];
end
if nargin < 1
    c1 = [0 0 1];
    c2 = [1 0 0];
end

brightness = repmat(linspace(brightlim(1), brightlim(2), n)', 1, 3);
up = repmat((0:(n-1))'/(n-1), 1, 3);
down = repmat(((n-1):-1:0)'/(n-1), 1, 3);
cmap = (down .* repmat(c1, n, 1) + up .* repmat(c2, n, 1)) .* brightness;



