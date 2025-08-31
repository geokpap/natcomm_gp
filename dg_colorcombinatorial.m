function cmap = dg_colorcombinatorial(hues, sats, vals)
% Produces a colormap containing all combinations of <hues>, <sats>, and
% <vals>.  All three args must be row vectors.

% This produces 70 mostly-readily-distinguishable colors:
% ([.7 1.25 1.8 3.2 4 5 6]/6, [1 .3], [1 .8 .6 .4 .25])

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

cmap = zeros(numel(hues) * numel(sats) * numel(vals), 3);
cidx = 1;
for s = sats
    for v = vals
        for h = hues
            cmap(cidx,:) = hsv2rgb(h, s, v);
            cidx = cidx + 1;
        end
    end
end