function s1 = lfp_SeqDelete(seq, set)
%s1 = lfp_SeqDelete(seq, set)

%   <seq> is a set of sequences as returned by lfp_TargSeq.  Searches for
%   all instances of each value in <set> and deletes them from all the
%   sequences. <seq> may contain empty sequences.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

s1 = seq;

for k = 1:length(s1)
    if ~isempty(s1{k})
        hits = find(ismember(s1{k}(:,2), set));
        s1{k}(hits,:) = [];
    end
end
