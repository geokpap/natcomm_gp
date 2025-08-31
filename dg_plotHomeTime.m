function isathome = dg_plotHomeTime( ...
    homebase, xy, name, Ts, windowlen, radius, startingpoint, hA, ...
    labelval)
% Helper function for dg_homeTime; can also be used to re-analyze results
% of dg_homeTime with different parameters without re-analyzing data to
% find homebase.  <isathome> is a logical column vector that is true for
% points where the animal was within <radius> of home base.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 7
    startingpoint = 1;
end
if nargin < 8
    hA = [];
end
if nargin < 9
    labelval = 1;
end

% Compute home base times
distance = sqrt( (xy(:,1) - homebase(1)).^2 ...
    + (xy(:,2) - homebase(2)).^2 );
isathome = distance - radius <= 1e-6;

% Smooth & plot
hometime = conv(ones(1,windowlen), double(isathome)) / windowlen;
hometime = hometime(windowlen:(end-windowlen+1));
if isempty(hA)
    figure;
    hA = axes;
    titlestr = sprintf('%s, Ts=%.3f s, win=%d samples, r=%.2f', ...
        name, Ts, windowlen, radius);
else
    titlestr = name;
end
plot(hA, ((windowlen:length(isathome)) + startingpoint - 1) * Ts, hometime);
ylim(hA, [-.05 1.05]);
switch labelval
    case 1
        title(hA, titlestr);
        ylabel(hA, 'Fractional Home Time');
        xlabel(hA, 'Elapsed Time, s');
    case 0
    otherwise
        ylabel(hA, labelval);
end
end


