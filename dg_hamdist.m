function dist = dg_hamdist(a, b, varargin)
%dist = dg_hamdist(a, b) computes Hamming distance between two nonnegative
%numbers.  If <b> is an array, then computes the Hamming distance between a
%and each element of <b>.
%OPTIONS
% 'min' - returns only the minimum distance between <a> and all the
%   elements of <b>.

%$Rev: 47 $
%$Date: 2010-03-11 17:36:30 -0500 (Thu, 11 Mar 2010) $
%$Author: dgibson $

minflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'min'
            minflag = true;
        otherwise
            error('dg_hamdist:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

dist = NaN(size(b));

for e = 1:numel(b)
    diffbits = bitxor(a,b(e));
    dist(e) = 0;
    for k = 0:63
        if bitand(diffbits, 1)
            dist(e) = dist(e) + 1;
        end
        diffbits = bitshift(diffbits, -1);
        if diffbits == 0
            break
        end
    end
end

if minflag
    dist = min(dist(:));
end

