function [result, C] = dg_analyzeGram(S, t, f, thresh, varargin)
%result = dg_analyzeGram(S, t, f, thresh)
% Analyzes bursts in the spectrogram or hhtogram <S>, whose time and
% frequency bin centers are given by <t> and <f>.
%INPUTS
% S: 2-D numeric array
% t: numeric vector of length size(S,2)
% f: numeric vector of length size(S,1)
% thresh: numeric scalar.  A burst is considered to be any region of <S>
%   whose values are at least <thresh>.  The string value 'auto' can also
%   be specified, in which case <thresh> is set to the maximum value in <S>
%   minus 3.
%OUTPUTS
% result: a numeric array containing one row for each burst.  Col 1 = start
%   time; col 2 = end time; col 3 = min freq; col 4 = max freq; col 5 = 0
%   if the burst boundary is completely contained within the gram, = 1 if
%   it touches the edge of the gram and may thus be incomplete.
% C: the set of isocontour lines at <thresh> that define the bursts, as
%   returned by Matlab func 'contourc'.
%OPTIONS
% 'findmax' - adds a sixth column to <result> containing the maximum value
%   within the burst.   This can be time consuming.

%$Rev: 61 $
%$Date: 2010-07-20 21:06:30 -0400 (Tue, 20 Jul 2010) $
%$Author: dgibson $

findmaxflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'findmax'
            findmaxflag = true;
        otherwise
            error('funcname:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

if isequal(thresh, 'auto')
    maxval = max(S(:));
    thresh = maxval - 3;
end

result = zeros(0, 5);

C = contourc(t, f, S, [thresh thresh]);
% For finding elements of S that are inside the burst boundary
X = repmat(reshape(t, 1, []), size(S,1), 1);
Y = repmat(reshape(f, [], 1), 1, size(S,2));

idx = 1;
while idx <= size(C,2)
    tvals = C(1, idx+1:idx+C(2,idx));
    fvals = C(2, idx+1:idx+C(2,idx));
    touches = any(tvals == t(1)) || any(tvals == t(end)) || ...
        any(fvals == f(1)) || any(fvals == f(end));
    result(end+1, 1:5) = ...
        [min(tvals) max(tvals) min(fvals) max(fvals) touches];
    if findmaxflag
        isinburst = inpolygon(X, Y, tvals, fvals);
        result(end, 6) = max(S(isinburst));
    end
    idx = idx + C(2,idx) + 1;
end