function y = dg_twopiece(x, b, m, r, k)
% Two-segment piecewise linear function.
% "All I want is two piece, and sheet on the bed."
%INPUTS
% b: intercept of left line
% m: slope of left line
% r: slope of right line
% k: breakpoint
%NOTES
% To calculate the intercept c of the right line, note that
% the two pieces intersect at the break point, i.e.
%  b + m * k = c + r * k
% and so
%  c = b + (m-r) * k

%$Rev:  $
%$Date:  $
%$Author: dgibson $

y = zeros(size(x));
c = b + (m-r) * k;
isL = x < k;
y(isL) = b + m * x(isL);
y(~isL) = c + r * x(~isL);
end
