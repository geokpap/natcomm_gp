function p = dg_firstordersignif(delimiters, seq)
%p = dg_firstordersignif(delimiters, seq)
% Computes the probability of the null hypothesis that all symbols -
% ignoring delimiters - in <seq> occur equally frequently.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

N = dg_symbolcounts(seq)';
if ~isempty(delimiters)
    N(delimiters+1) = [];
end
p = dg_chi2test2(N, true);