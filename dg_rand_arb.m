function nums = dg_rand_arb(binedges, counts, N)
% Generates a pseudo-random set of numbers that fit an arbitrary
% numerically specified distribution.
%INPUTS
% binedges: vector of bin edges of a histogram specifying the desired
%   distribution.
% counts: vector of a histogram specifying the distribution.  The rule for
%   handling the bin edges is the same as in Matlab 'histc', and the length
%   of <counts> is one less than the length of <binedges> (again, for
%   consistency with 'histc', where the last element of the counts is
%   normally a useless value).
% N: number of random values to create.
%OUTPUTS
% nums: column vector of randomly generated ages, in time steps.
%NOTES
% Follows "inversion method" of
% http://web.mit.edu/urban_or_book/www/book/chapter7/7.1.3.html, which is
% to map uniformly distributed random numbers to target values using the
% inverse of the cumulative distribution function (CDF) of the desired
% probability distribution.  This appears to be the same as the Matlab
% 'icdf' method, except that we don't mess with toolbox objects.

if ~(isvector(binedges) && isvector(counts))
    error('dg_rand_arb:vector', ...
        '<binedges> and <counts> must both be vectors.');
end
if numel(counts) ~= numel(binedges) - 1
    error('dg_rand_arb:counts', ...
        '<counts> must contain one less element than <binedges>.');
end

% <mycdf> will be the desired CDF, properly normalized to run from 0 to 1.
% By definition, the first value is zero and the last value is 1. In order
% to be invertible, <mycdf> cannot contain any repeated values. That means
% that <counts> must not contain any zero values, so we need to delete both
% the zero values and the bins that correspond to them.
iszerobin = [ counts(:)==0 ];
counts(iszerobin) = [];
binedges([false; iszerobin]) = [];

% We start with a value of zero so that each point in <mycdf> corresponds
% to one value in <binedges>.
mycdf = [ 0; cumsum(counts(:)) / sum(counts(:)) ];


% Now we interpret results of 'rand' as values in <mycdf>, and linearly
% interpolate the point within the bin that corresponds to the inverse of
% the <mycdf> values.
r = rand(1, N);
nums = interp1(mycdf, binedges, r, 'linear');

