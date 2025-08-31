function [hF, hA] = dg_plotglm(stats, xlimits, hA)
% Plots decision boundary with confidence limits from output of 
%   glmfit(X, y, 'binomial', 'link', 'logit')'
% where <X> has two columns (predictors) and <y> is logical.
%INPUTS
% stats: as returned by 'glmfit'.
% xlimits: the minimum and maximum x values over which to plot the curves.
% hA: axes handle into which to plot; may be omitted, in which case a new
%   figure containing a new axes object is created.
%OUTPUTS
% hF: figure handle that owns <hA>.
% hA: the axes containing the new curves.
%NOTES
% The CLs are constructed on the basis of the parameters' values being
% normally distributed around means <stats.beta> with standard deviations
% <stats.se>.

if nargin < 3
    hF = figure;
    hA = axes('Parent', hF);
end

% the decision boundary line (where p = 0.5):
% "The coefficient of the constant term is the first element of b."
b = stats.beta;
plot(hA, xlimits, -(b(1) + b(2) * xlimits) / b(3));
set(hA, 'NextPlot', 'add');
se = stats.se;
numobs = 2000;
numpts = 100;
coeffs = repmat(b, 1, numobs);
coeffs = coeffs + randn(size(coeffs)) .* repmat(se, 1, numobs);
xvals = linspace(xlimits(1), xlimits(2), numpts);
yvals = NaN(numobs, numpts);
for obsidx = 1:numobs
    yvals(obsidx, :) = -(coeffs(1, obsidx) + ...
        coeffs(2, obsidx) .* xvals) ./ coeffs(3, obsidx);
end
CLs = prctile(yvals, [2.5 97.5]);
plot(hA, xvals, CLs');
