function s2 = lfp_seqFirstOnly(s1)
%LFP_SEQFIRSTONLY reduces a sequence to just the first occurrence of each
%target within the sequence.
%s2 = lfp_seqFirstOnly(s1)

%   <s1> is a set of sequences as returned by lfp_TargSeq.  <s2> is
%   returned in the same format.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

s2 = cell(size(s1));
for k = 1:length(s1)
    if ~isempty(s1{k})
        elements = unique(s1{k}(:,2))';
        keepers = [];
        for el = elements
            hits = find(s1{k}(:,2) == el);
            keepers(end+1) = hits(1);
        end
        s2{k} = s1{k}(sort(keepers),:);
    end
end
