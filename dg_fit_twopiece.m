function [B, xisect, meanerror] = dg_fit_twopiece(x, y, xisect)
% Fit two-segment piecewise linear function to data.
% "All I want is two piece, and sheet on the bed."
%INPUTS
% x: array of x coordinates to fit with two lines; must be in strictly
%   ascending order.
% y: array of y coordinates corresponding to <x>.
% xisect: initial guess as to where the two lines intersect.
%OUTPUTS
% B: fit parameters.  <B(1,:)> is the y-intercepts; <B(2,:)> is the slopes.
%   <B(:,1)> is the left-hand segment fitted to values less than
%   <xisect>; <B(:,2)> is the right-hand segment fitted to x values
%   greater thatn <xisect>.
% xisect: x coordinate of the intersection of the fitted segments.
% meanerror: mean error per point.
%NOTES
%   The strategy is to begin with the left half of the left portion of
% <y> and the right half of the right portion of <y>, fit each of those,
% and calculate their intersection point, yielding an updated value for
% <xisect>.  The process then repeats iteratively until either <xisect>
% doesn't change by at least the smallest difference between successive x
% values, or the summed absolute error of the two fits doesn't get smaller,
% or <xisect> hits the start or end of <y>.  On each iteration, the portion
% of data excluded from the fit is cut in half.  The summed absolute error
% is always computed over the entirety of the two fitted segments.
%   The 'while' loop is designed to be executed at least once to
% start with.

%$Rev: 277 $
%$Date: 2021-08-20 16:25:43 -0400 (Fri, 20 Aug 2021) $
%$Author: dgibson $

if numel(x) ~= numel(y)
    error('dg_fit_twopiece:numel', ...
        '<x> and <y> must have the same number of elements.');
end
if xisect > max(x) || xisect < min(x)
    error('dg_fit_twopiece:xisect', ...
        '<xisect> is out of bounds.');
end
x = x(:); % column vector
y = y(:); % column vector
if any(diff(x) <= 0)
    error('dg_fit_twopiece:x', ...
        '<x> is not monotonically increasing.');
end
tolx = min(diff(x));
newxisect = Inf; % to be calculated after fitting the 2 segments
lastL = 2; % to be calculated after fitting the 2 segments
firstR = 1; % to be calculated after fitting the 2 segments
abserr = [];
prevtoterr = [];
toterr = Inf; % just to keep mlint happy
enddata = numel(y);
B = [];
newB = [];
iternum = 0;
while abs(newxisect - xisect) > tolx && ...
        (isempty(prevtoterr) || (toterr < prevtoterr)) ...
        && ~(lastL < 2 || firstR > enddata - 1)
    iternum = iternum + 1;
    endL = round(lastL*(1 - 1/2^iternum));
    startR = enddata - round((enddata - firstR)*(1 - 1/2^iternum));
    B = newB;
    % compute slopes and offsets for both linear fits.
    % <B(1,:)> is the y-intercepts, <B(2,:)> is the slopes:
    newB(:, 1) = [ones(endL, 1) x(1:endL)] \ y(1:endL);
    newB(:, 2) = [ones(enddata - startR + 1, 1) x(startR:enddata)] ...
        \ y(startR:enddata);
    prevtoterr = toterr;
    if ~isinf(newxisect)
        xisect = newxisect;
        abserr = newabserr;
    end
    % Find x-coordinate <newxisect> of new intersection point:
    newxisect = round((newB(1,1) - (newB(1,2))) / ...
        (newB(2,2) - newB(2,1)));
    lastL = find(x < newxisect, 1, 'last'); % last <x> value before <newxisect>
    if isempty(lastL)
        lastL = 2;
        newxisect = x(1);
    end
    firstR = find(x > newxisect, 1); % first <x> value after <newxisect>
    if isempty(firstR)
        firstR = enddata - 1;
        newxisect = x(end);
    end
    % Find new total error:
    newabserr = [
        abs( y(1:lastL) - ...
        (newB(1,1) + newB(2,1)*(1:lastL)') )
        abs( y(firstR:enddata) - ...
        (newB(1,2) + newB(2,2)*(firstR:enddata)') )
        ];
    toterr = sum(newabserr);
end
if isempty(B)
    % The loop has only executed once, in which case the initial value of
    % <xisect> is as good as it gets.  We therefore use the fit
    % parameters that resulted from that value.
    B = newB;
    abserr = newabserr;
end
meanerror = mean(abserr);

