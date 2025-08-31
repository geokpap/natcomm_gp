function s2 = lfp_TargOnly(s1)
%s2 = lfp_TargOnly(s1)
%   <s1> is a set of sequences as returned by lfp_TargSeq.  Returns a
%   cell column vector containing the same data but with timestamps
%   removed, so sequences contain targets only.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

s2 = cell(numel(s1),1);
for k=1:length(s1)
    if isempty(s1{k})
        s2{k} = [];
    else
        s2{k} = s1{k}(:,2);
    end
end