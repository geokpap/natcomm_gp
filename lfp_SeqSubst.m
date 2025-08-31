function s2 = lfp_SeqSubst(s1, v1, v2)
%s2 = lfp_SeqSubst(s1, v1, v2)

%   <s1> is a set of sequences as returned by lfp_TargSeq.  Searches for
%   instances of the value v1 as a target ID, and replaces them with v2.
%   If v2 is [], then the instances of v1 are deleted from the sequence
%   entirely.  <s1> may contain empty sequences.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

if ~isequal(size(v1), [1 1]) || ~(isempty(v2) || isequal(size(v2), [1 1]))
    error('lfp_SeqSubst:badarg', ...
        '<v1> and <v2> must be scalars' );
end

s2 = s1;

for k = 1:length(s2)
    if ~isempty(s2{k})
        hits = find(s2{k}(:,2) == v1);
        if isempty(v2)
            s2{k}(hits,:) = [];
        else
            s2{k}(hits,:) = v2;
        end
    end
end
