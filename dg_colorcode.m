function cmap = dg_colorcode
% Produces a colormap containing many readily distinguishable colors.

% This produces 70 mostly-readily-distinguishable colors:

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

hues = {
    [0 .09 .13 .18 .23 .43 .53 .6 .72 .8 .88]   % 1:11
    [0 1 2 3 4 5 6]/7                           % 12:18
    [0 .1 .2 .3 .5 .58 .72 .8 .9]               % 19:27
    [0 1 1.6 3 4 5]/6                           % 28:33
    [0 .1 .2 .3 .5 .58 .72 .8 .9]                 % 34:42
    [0 1 2 3 4]/5
    [0 .1 .2 .3 .5 .6 .8 .9]
    [0 1 2 3]/4
    [0 .65]
    [.2 .5 .8]
    };
sats = [1 .3];
vals = [1 .8 .6 .44 .3];

cidx = 1;
for v = 1:5
    for s = 1:2
        myhues = hues{2*(v-1) + s};
        for h = myhues
            cmap(cidx,:) = hsv2rgb(h, sats(s), vals(v));
            cidx = cidx + 1;
        end
    end
end