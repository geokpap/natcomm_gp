function h = dg_grubbstest(x, alpha)
%INPUTS
% x: a set of numerical observations.
% alpha: statistical significance level.
%OUTPUT
% h: true if the null hypothesis that there are no outliers can be rejected
%   at significance level <alpha>.
%NOTES
% From formulas in Wikipedia entry "Grubbs' test for outliers".

%$Rev: 255 $
%$Date: 2016-12-19 16:48:11 -0500 (Mon, 19 Dec 2016) $
%$Author: dgibson $

G = max(abs(x - mean(x))) / std(x);
N = length(x);
tcrit = -tinv(alpha/(2*N), N-2);
h = G > (N-1) * sqrt(tcrit^2/(N-2+tcrit^2)) / sqrt(N);
end

