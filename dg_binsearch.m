function index = dg_binsearch(v, c, h, h2)
%index = dg_binsearch(v, c)
%   In the sorted increasing vector of numbers v, returns the index of the
%   first value that is greater than c.  If there is no such number in v,
%   returns length(v)+1.  If v is empty, returns v.
%index = dg_binsearch(v, c, h)
%index = dg_binsearch(v, c, h, h2)
%   h is a hint to use as the index for the first comparison instead of the
%   center value.  The purpose here is that if the desired element is
%   estimated to be much closer to one bound than the other, h can be used
%   to prune the search interval without actually throwing away the data on
%   the far side of h; for this purpose, h should be a little further
%   from the closer bound than the estimated position of the target element.
%   If h2 is given, it is handled exactly like h.  This makes it possible
%   to immediately bracket a small subrange that probably contains the
%   target value.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if isempty(v)
    index = v;
    return
end

if c > v(end)
    index = length(v) + 1;
    return
end
if v(1) > c
    index = 1;
    return
end

lowerbound = 1;
upperbound = length(v);
if nargin > 2
    if v(h) > c
        upperbound = h;
    else
        lowerbound = h;
    end
end
if nargin > 3
    if v(h2) > c
        upperbound = h2;
    else
        lowerbound = h2;
    end
end
while upperbound > lowerbound + 1
    center = round((upperbound + lowerbound)/2);
    if v(center) > c
        upperbound = center;
    else
        lowerbound = center;
    end
end
index = upperbound;
