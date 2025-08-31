function counts = dg_symbolcounts(seq)
%counts = dg_symbolcounts(seq)
% Returns a column vector <counts> of the number of instances of each value
% in the vector <seq>.  <counts> is indexed by the values in <seq> plus
% one, so <seq> may only contain values greater than 0 and less than some
% maximum number.  That number is theoretically the maximum integer value
% available on the machine, but practically other bad things tend to happen
% (like running out of memory), so the upper limit is set substantially
% lower, and its value is included in the error message if it is violated.
% length(counts) = max(seq) + 1.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if ~isnumeric(seq) || ~isequal(fix(seq), seq)
    error('dg_symbolcounts:badtype', ...
        '<seq> must contain integer values.' );
end
[str,maxsize] = computer;
maxval = round((maxsize + 1)/128) - 1;
if any(seq > maxval)
    error('dg_symbolcounts:badval1', ...
        '<seq> contains value(s) that are greater than %d.', maxval );
end
if any(seq < 0)
    error('dg_symbolcounts:badval2', ...
        '<seq> contains value(s) less than zero.' );
end

seq = reshape(seq, [], 1) + 1;
counts = zeros(max(seq), 1);
for k = unique(seq)'
    counts(k) = sum(seq == k);
end
