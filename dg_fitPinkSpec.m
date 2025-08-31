function [pinkBL, a, b] = dg_fitPinkSpec(BL, varargin)
%pinkBL = dg_fitPinkSpec(BL)
% Creates a new baseline spectrum structure by fitting a smooth pink noise
% spectrum to <BL>.  There are two free parameters: slope and amplitude.
% Fitting is done by the Curve Fitting Toolbox using p = a*f^b, where f is
% frequency, p is power, a is fitted amplitude and b is fitted slope.  The
% first two points are ignored in doing the fit and are simply extrapolated
% linearly back to x=0.
%INPUTS
% BL: as returned by lfp_BLspectrum; requires the following fields:
%   BL.sum: a column vector containing the sum of all linear spectral power
%       observations; the mean spectrum computed as BL.sum/BL.N is exactly
%       as computed by lfp_mtspectrum(..., 'avg').
%   BL.N: the number of observations (windows times tapers) that produced
%       the sums.
%   BL.f: the frequency represented by each point in BL.sum and BL.sumsqrs.
%       It is assumed that BL.f(1) = 0.
%OUTPUTS
% pinkBL: a partial baseline spectrum structure, with fields
%   'f' - same as BL.f
%   'N' - set to 1.
%   'sum' - the fitted pink spectrum (same length as <pinkBL.f>).  Values
%       are NaN at frequencies greater than <freqlim(2)> (see OPTIONS) and
%       are linearly extrapolated from the first two points inside
%       <freqlim> at frequencies below <freqlim(1)>.
% a, b: the parameters of the fitted curve p = a*f^b.
%OPTIONS
% 'fixb', fixb - fixes the value of <b> in the fitted p = a*f^b, where
%   b = -<fixb>.  In this case the fitting is done analytically because I
%   can't figure out how to get the curve fitting toolbox to work correctly
%   when I'm using an arbitrary function instead of one of the built-in
%   functions, and <cfun> = [];
% 'freqlim', freqlim - <freqlim> is a two-element numeric vector specifying
%   the lower limit (f >= freqlim(1)) and upper limit (f < freqlim(2)) of
%   the frequency range to fit.  All f < freqlim(1) are extrapolated
%   linearly.  Default range of points to fit is BL.f(3:end).
% 'logfit' - tranforms to log-log scale before doing the fit so that
%   instead of minimizing squared absolute error, minimizes square error
%   relative to value.  In this case, <cfun> is a first order polynomial
%   function fitted to the natural logarithms of the frequencies and
%   spectrum in <BL>.

%$Rev: 251 $
%$Date: 2016-09-09 15:52:44 -0400 (Fri, 09 Sep 2016) $
%$Author: dgibson $

if any(BL.sum < 0)
    error('dg_fitPinkSpec:negBL', ...
        '<BL.sum> contains negative values.  It must contain linear power, not log power.');
end

argnum = 1;
fixb = [];
freqlim = [];
weighting = 'linear';
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'fixb'
                argnum = argnum + 1;
                fixb = varargin{argnum};
            case 'freqlim'
                argnum = argnum + 1;
                freqlim = varargin{argnum};
            case 'logfit'
                weighting = 'log';
            otherwise
                error('dg_fitPinkSpec:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('dg_fitPinkSpec:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if isempty(freqlim)
    idxrange = 3:length(BL.f);
else
    idxrange = find(BL.f >= freqlim(1) & BL.f < freqlim(2));
end
if isempty(idxrange)
    error('dg_fitPinkSpec:badfreqlim', ...
        'There are no frequency points in the %s range.', ...
        dg_thing2str(freqlim));
end
if idxrange(1) == 1
    % Under no circumstances can we attempt to fit a power law to zero Hz!
    idxrange(1) = [];
end
f = reshape(BL.f(idxrange), [], 1);
p = reshape(BL.sum(idxrange) / BL.N, [], 1);
pinkBL.f = BL.f;
pinkBL.N = 1;
pinkBL.sum = NaN(size(BL.sum));
switch weighting
    case 'linear'
        % Fit Model: p = a * BL.f(idxrange) .^ b
        if isempty(fixb)
            [cfun,gof,output] = fit( f, p, 'Power1' ); %#ok<ASGLU>
            a = cfun.a;
            b = cfun.b;
        else
            cfun = [];
            b = -fixb;
            % "Trust me" (I did the math): the a that minimizes the squared
            % error is:
            a = sum((f.^b) .* p) / sum(f.^(2*b));
            output.exitflag = 1;
        end
    case 'log'
        % Transform coordinates: y = log(p); x = log(f).
        % Transform Model: log(p) = log(a*f^b) = log(a) + b*log(f)
        % Substitute x and y in transformed model: y = log(a) + b*x.
        % Or in Matlab-speak, where we are fitting the line p1*x + p2:
        % b = p1, log(a) = p2 and so a = exp(p2).
        x = log(f);
        y = log(p);
        % consequently, a = exp(c).
        if isempty(fixb)
            [cfun,gof,output] = fit( x, y, 'poly1' ); %#ok<ASGLU>
            a = exp(cfun.p2);
            b = cfun.p1;
        else
            b = -fixb;
            % The least squares fit of a constant to a set of data is
            % simply the mean of the data.  So p2 = mean(y - b*x), and
            % therefore a = exp(mean(y - b*x)).
            a = exp(mean(y - b*x));
            output.exitflag = 1;
        end
end
pinkBL.sum(idxrange(1):end) = a * BL.f(idxrange(1):end) .^ b;
if output.exitflag <= 0
    error('dg_fitPinkSpec:fit', ...
        'Curve fitting failed.');
end
% Extrapolate the low frequency point(s) linearly as y = mx + c:
m = diff(pinkBL.sum(idxrange(1:2))) / diff(BL.f(idxrange(1:2)));
c = pinkBL.sum(idxrange(1)) - m * BL.f(idxrange(1));
pinkBL.sum(1:idxrange(1)-1) = m * BL.f(1:idxrange(1)-1) + c;


