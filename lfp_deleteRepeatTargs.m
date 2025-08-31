function s2 = lfp_deleteRepeatTargs(s1)
%LFP_DELETEREPEATTARGS removes repeat entries in a target sequence

%   <s1> and <s2> are in the format returned by lfp_TargSeq.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

s2 = s1;
for trialidx = 1:length(s1)
    if ~isempty(s1{trialidx})
        repeats = find(s1{trialidx}(2:end,2) == s1{trialidx}(1:end-1,2));
        s2{trialidx}(repeats+1,:) = [];
    end
end