function [mantissa, exp] = dg_findScale(data, varargin)
%[mantissa, exp] = dg_FindScale(data)
%INPUTS
% data:  a vector of numbers.
%OUTPUTS
% Finds the smallest number of the form [ 1 2 5 ] * 1.0e+N that is
% greater than or equal to max(abs(data)).  The integer from the list is
% returned as <mantissa>, and the value of N is returned as <exp>.  If
% max(abs(data)) is 0, returns [1 0].
%OPTIONS
% 'levels', N - returns column vectors in decreasing order of
%   mantissa*10^exp.  If N = 0, then enough levels are returned to span the
%   entire range of data.  Otherwise, at most N levels are returned.

%$Rev: 105 $
%$Date: 2011-04-29 16:25:17 -0400 (Fri, 29 Apr 2011) $
%$Author: dgibson $

numlevels = 1;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'levels'
            argnum = argnum + 1;
            numlevels = varargin{argnum};
        otherwise
            error('dg_FindScale:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
if ~isnumeric(numlevels) || ~isequal(size(numlevels), [1 1])
    error('dg_findScale:levels', ...
        'Number of levels must be a numeric scalar.');
end

datasize = size(data);
if sum(datasize > 1) > 1
    error('dg_FindScale:BadArg', 'The argument must be a vector');
end
maxdata = max(abs(data));
if maxdata == 0
    mantissa = 1;
    exp = 0;
    return;
end
exp = floor(log10(maxdata));
done = false;
mantissavals = [ 1 2 5 ];
for midx = 1:length(mantissavals)
    if mantissavals(midx)*10^exp >= maxdata
        done = true;
        break
    end
end
if ~done
    % maxdata is between 5*10^exp and 10^(exp+1):
    exp = exp + 1;
    midx = 1;
end
mantissa = mantissavals(midx);

mindata = min(abs(data));
while mantissa(end)*10^exp(end) > mindata && ( ...
        numlevels == 0 || length(exp) < numlevels )
    if midx == 1
        midx = length(mantissavals);
        exp(end+1, 1) = exp(end) - 1; %#ok<*AGROW>
    else
        midx = midx - 1;
        exp(end+1, 1) = exp(end);
    end
    mantissa(end+1, 1) = mantissavals(midx);
end
    
