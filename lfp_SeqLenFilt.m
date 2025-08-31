function s2 = lfp_SeqLenFilt(s1, n1, n2)
%s2 = lfp_SeqLenFilt(s1, v1, v2)
%   <s1> is a set of sequences as returned by lfp_TargSeq.  Returns a
%   copy of <s1> containing only those sequences whose length (i.e. the
%   number of targets in the sequence) is between n1 and n2 inclusive.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

goodlength = zeros(size(s1));
for k = 1:length(s1)
    goodlength(k) = (length(s1{k}) >= n1) && (length(s1{k}) <= n2);
end
s2 = s1(find(goodlength));
