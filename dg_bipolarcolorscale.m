function cmap = dg_bipolarcolorscale(clim, n, numyellows, numreds)
%cmap = dg_bipolarcolorscale(clim, n, numyellows, numreds)
%cmap = dg_bipolarcolorscale(clim, n, numyellows)
%cmap = dg_bipolarcolorscale(clim, n)
%cmap = dg_bipolarcolorscale(clim)
%cmap = dg_bipolarcolorscale
% Returns a colormap with <n> elements to represent the range of values
% <clim> in imagesc such that zero maps to black, increasingly negative
% values map to increasing brightnesses of blue, and positive values map to
% an incandescent-style color scale, i.e. black to red to orange to yellow
% to white.  <numyellows> is the number of non-saturated yellows up to and
% including white.
%DEFAULTS
% n: 64
% clim: [-1 1]
% numyellows: one fourth of the number of hot colors

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 2
    n = 64;
end
if nargin < 1
    clim = [-1 1];
end

if numel(clim) ~= 2 || clim(1) >= clim(2)
    error('dg_bipolarcolorscale:clim', ...
        '<clim> must be a 2-element vector with clim(2) > clim(1).' );
end

range = clim(2) - clim(1);
binrange = range/n;
zerobin = ceil(-clim(1)/binrange);
if zerobin > 0
    zerocolor = [0 0 0];
    numhots = n - zerobin;
%     hotbins = (zerobin+1) : n ;
else
    zerocolor = [];
    numhots = n;
%     hotbins = 1 : n ;
end
numblues = max(zerobin-1, 0);
if numblues > 0
    blues = [ zeros(numblues, 2) linspace(1, 1/numblues+0.2, numblues)' ];
else
    blues = [];
end

if nargin < 4
    % number of brightnesses from black to red:
    numreds = round(numhots/4);
end
if nargin < 3
    % number of saturations from yellow to white:
    numyellows = round(numhots/4);
end
numoranges = numhots - numreds - numyellows;
if numyellows > 0
    yellows = ...
        [ ones(numyellows, 2) linspace(1/numyellows+0.2, 1, numyellows)' ];
else
    yellows = [];
end
hots = [ linspace(1/numreds, 1, numreds)' zeros(numreds, 2)
    ones(numoranges, 1) linspace(1/numoranges+0.2, 1, numoranges)' zeros(numoranges, 1)
    yellows ];
cmap = [ blues
    zerocolor
    hots ];
