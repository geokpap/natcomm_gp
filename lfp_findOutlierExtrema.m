function TS = lfp_findOutlierExtrema(filenum, varargin)
%TS = lfp_findOutlierExtrema(filenum)
% For use with lfp_createEvents.  Marks extrema (defined by
% dg_findFlattops) that are in the tails of the distributions of extrema.
% Statistics are computed separately for maxima and minima.  Can use a
% Gaussian model (based on the mean and standard deviation; this is the
% default method) or the actual empirical distribution.  The Gaussian model
% is best to use when there might not be outliers that should be excluded,
% whereas the actual distirbution method will always exclude the specified
% proportion of points.
%INPUTS
% filenum: CSC channel to analyze.
%OUTPUTS
% TS: timestamps of first sample of each (potentially flat-topped) peak.
%OPTIONS
% 'pval', P - <P> is a one-sided p-value for defining the outliers based on
%   the actual distribution of values.
% 'sigmas', N - <N> is the number of standard deviations from the mean for
%   defining outliers; default = 5.  This corresponds to p < 3e-7, which in
%   a 2 or 3 hour session sampled at 1kHz is expected to produce about 3
%   false alarms from normally distributed samples.  <N> = 6 would
%   eliminate false alarms in 99% of sessions, but may also miss serious
%   artifacts.

%$Rev: 299 $
%$Date: 2013-04-24 16:04:30 -0400 (Wed, 24 Apr 2013) $
%$Author: dgibson $

global lfp_Samples

pval = 0;
nsigmas = 5;

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'nsigmas'
            argnum = argnum + 1;
            nsigmas = varargin{argnum};
        case 'pval'
            argnum = argnum + 1;
            pval = varargin{argnum};
        otherwise
            error('lfp_findOutlierExtrema:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

vmaxidx = dg_findFlattops(lfp_Samples{filenum}(:));
maxvals = lfp_Samples{filenum}(vmaxidx);
vminidx = dg_findFlattops(-lfp_Samples{filenum}(:));
minvals = lfp_Samples{filenum}(vminidx);
if pval
    maxthresh = prctile(maxvals, 100*(1-pval));
    minthresh = prctile(minvals, 100*pval);
else
    meanmax = mean(maxvals);
    meanmin = mean(minvals);
    sdmax = std(maxvals);
    sdmin = std(minvals);
    maxthresh = meanmax + nsigmas*sdmax;
    minthresh = meanmin - nsigmas*sdmin;
end
allextrema = [reshape(vmaxidx, [], 1)
    reshape(vminidx, [], 1)];
isbadextremum = lfp_Samples{filenum}(allextrema) > maxthresh ...
    | lfp_Samples{filenum}(allextrema) < minthresh;
TS = lfp_index2time(allextrema(isbadextremum));

