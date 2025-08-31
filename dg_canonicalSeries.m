function str = dg_canonicalSeries(v)
%DG_CANONICALSERIES produces a compact string representation of a vector
%of integers.
%str = dg_canonicalSeries(v)

%  Collapses any consecutive integer subsequence of v into the equivalent
%  Matlab colon expression; non-consecutive elements are separated by
%  spaces.  The entire string is surrounded by square brackets.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if isempty(v) 
    str = '[]';
    return
end
if length(v)<2
    str = mat2str(v);
    return
end
prevlast = 0;
str = '[';
[first last] = dg_findconsecutive(v, 1);
while ~isempty(last)
    if first > 1
        str = [str sprintf(' %d', v(prevlast+1 : first-1))];
    end
    str = [ str sprintf(' %d:%d', v(first), v(last)) ];
    prevlast = last;
    [first last] = dg_findconsecutive(v, prevlast+1);
end
if prevlast < length(v)
    str = [str sprintf(' %d', v(prevlast+1 : end))];
end
str = [str ' ]'];


