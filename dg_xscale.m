function yscaled = dg_xscale(x,y,m,varargin)
% Scales the function y=f(x) along the x axis by a factor of m.  Scale
% factors greater than 1 result in the plot of y becoming broader.
% <yscaled> is the value of the scaled rendition of y at all of the same x
% values as were in the original x.
%INPUTS
% x: vector of x values
% y: vector of y values, same length as x
% m: scalar scale factor
%OUTPUT
% yscaled: vector of scaled y values, same size & shape as y
%OPTIONS
% 'propagate' - when m < 1, the right-hand values are filled in with NaNs
%   by default; this option propagates the last non-NaN value.  When m > 1,
%   if the first x value is nonzero, there will be initial NaNs; this
%   option replaces them by linearly extrapolating the first two points
%   back to zero.

%$Rev: 88 $
%$Date: 2010-10-22 16:47:43 -0400 (Fri, 22 Oct 2010) $
%$Author: dgibson $

propagateflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'propagate'
            propagateflag = true;
        otherwise
            error('dg_xscale:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

ysize = size(y);
xscaled = x / m;
% In the case where m > 1 and the first data point is not zero, instead of
% propagating on the left, we extrapolate back to x=0 in order to prevent
% interp1 from returning NaNs at the initial point(s):
if propagateflag && m > 1 && x(1) ~= 0
    slope = (y(2) - y(1)) / (x(2) - x(1));
    y0 = y(1) - x(1) * slope;
    y = [ y0 reshape(y,1,[]) ];
    x = [ 0 reshape(x,1,[]) ];
end
yscaled = interp1(x, y, xscaled);
if propagateflag
    nans = isnan(yscaled);
    nanidx = find(nans);
    if m <= 1
        if ~isempty(nanidx)
            if nanidx(1) < 2
                error('dg_xscale:badnan', ...
                    'There is a nan in the first data point');
            else
                yscaled(nans) = yscaled(nanidx(1)-1);
            end
        end
    end
end
yscaled = reshape(yscaled, ysize);

