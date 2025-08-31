function [result, headers, binedges] = ...
    lfp_assocHist(xfnum, yfnum, varargin)
% <result> is a column vector.  <xfnum> is a filenum whose contents are
%   used as a basis for binning the data.  <yfnum> is a filenum whose
%   contents are used as the values of the binned data.  If <yfnum> is
%   zero, then <result> contains the number of points in each bin.  That
%   is, <result> is what would be returned by hist(lfp_Samples{xfnum}).  If
%   <yfnum> is not zero, then <result> contains the mean value of
%   lfp_Samples{yfnum} for all the samples that are in each bin of
%   lfp_Samples{xfnum}.  For historical reasons, <result> is zero for any
%   bins that contain no data points, except in the case where every bin
%   contains no data and every element of <result> is NaN.
% <headers> is cell string row vector containing a column header for each
%   column of <result>.
% <binedges> is a list of edge values for the bins that are used as in
%   Matlab's histc function (edges(k)<= x(i) < edges(k+1)).  Bin centers
%   can be computed as: 
%       binctrs = (binedges(1:end-1) + binedges(2:end))/2
%OPTIONS
% 'binwidth', b - specifies the width of the bins in the same units as the
%   data contained in <xfnum>.  Default is (edgeN - edge1)/nbins. 
% 'count' - the number of non-NaN-valued data points in the bin
% 'edge1', e - sets the value of binedges(1) equal to <e>.  Default is
%   min(lfp_Samples{xfnum}).
% 'edgeN', e - sets the value of binedges(end) equal to <e>.  Default is
%   max(lfp_Samples{xfnum}).
% 'median' - returns an additional column containing medians.
% 'nbins', n - specifies the number of bins (rows) in <result>.  Default is
%   100.
% 'sem' - returns an additional column containing standard errors of the
%   means.
% 'std' - returns an additional column containing standard deviations.
% 'trialonly' - limits the analysis to the time periods belonging to each
%   enabled trial (see lfp_enabledTrials), including the endpoints listed
%   in lfp_TrialIndex columns 3 and 4. 
%NOTES
% If more than one of the additional result column options is given,
%   then the columns are returned in the order in which the options were
%   given.  These options have no effect if <yfnum> is 0.
% It is not possible to specify all of 'binwidth', 'edge1', 'edgeN',
%   'nbins' because they are constrained by the equation:
%       binwidth = (edgeN - edge1)/nbins
%   If they are all specified, then 'nbins' is ignored.
% None of the statistics returned in <result> is nan-tolerant aside from
%   'count', so even just one data point with the value NaN will produce
%   a Nan-valued bin.

%$Rev: 53 $
%$Date: 2009-03-20 18:16:34 -0400 (Fri, 20 Mar 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if yfnum
    headers{1,1} = 'mean';
else
    headers{1,1} = 'count';
end

binwidth = [];
countflag = false;
edge1 = [];
edgeN = [];
medianflag = false;
nbins = [];
result = zeros(0, 1);
semflag = false;
stdflag = false;
trialonlyflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'binwidth'
            argnum = argnum + 1;
            binwidth = varargin{argnum};
        case 'count'
            countflag = true;
            countcol = size(result,2) + 1;
            result = zeros(0, countcol);
            headers{countcol} = 'count';
        case 'edge1'
            argnum = argnum + 1;
            edge1 = varargin{argnum};
        case 'edgeN'
            argnum = argnum + 1;
            edgeN = varargin{argnum};
        case 'median'
            medianflag = true;
            mediancol = size(result,2) + 1;
            result = zeros(0, mediancol);
            headers{mediancol} = 'median';
        case 'nbins'
            argnum = argnum + 1;
            nbins = varargin{argnum};
        case 'sem'
            semflag = true;
            semcol = size(result,2) + 1;
            result = zeros(0, semcol);
            headers{semcol} = 'sem';
        case 'std'
            stdflag = true;
            stdcol = size(result,2) + 1;
            result = zeros(0, stdcol);
            headers{stdcol} = 'std';
        case 'trialonly'
            trialonlyflag = true;
        otherwise
            error('lfp_assocHist:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~ismember(xfnum, lfp_ActiveFilenums)
    error('lfp_assocHist:xfilenum', ...
        'There is no filenum %d loaded', xfnum );
end

if yfnum && ~ismember(xfnum, lfp_ActiveFilenums)
    error('lfp_assocHist:yfilenum', ...
        'There is no filenum %d loaded', yfnum );
end

if trialonlyflag
    samplebounds = lfp_TrialIndex( ...
        lfp_enabledTrials(find(lfp_SelectedTrials)), 3:4 );
else
    samplebounds = [1 numel(lfp_Samples{xfnum})];
end
include = false(size(lfp_Samples{xfnum}));
for row = 1:size(samplebounds,1)
    include(samplebounds(row,1):samplebounds(row,2)) = true;
end
if yfnum
    data = [ reshape(lfp_Samples{xfnum}(include), [], 1) ...
        reshape(lfp_Samples{yfnum}(include), [], 1) ];
else
    data = reshape(lfp_Samples{xfnum}(include), [], 1);
end
data = sortrows(data);

% If other constraints are sufficient to determine nbins, calculate its
% value; otherwise use its default value.  Either way, nbins is sure to
% have a value after this.
if ~(isempty(binwidth) || isempty(edge1) || isempty(edgeN))
    nbins = round((edgeN - edge1)/binwidth);
else
    if isempty(nbins)
        nbins = 100;
    end
end

if isempty(binwidth)
    if ~(isempty(edge1) || isempty(edgeN))
        binwidth = (edgeN - edge1)/nbins;
    end
end

if isempty(edge1)
    if ~(isempty(binwidth) || isempty(edgeN))
        edge1 = edgeN - nbins * binwidth;
    else
        edge1 = min(data(:,1));
    end
end

if isempty(edgeN)
    if ~(isempty(binwidth) || isempty(edge1))
        edgeN = edge1 + nbins * binwidth;
    else
        edgeN = max(data(:,1));
    end
end

if isempty(binwidth)
    % edge1 and edgeN MUST have values at this point:
    binwidth = (edgeN - edge1)/nbins;
end

binedges = linspace(edge1, edgeN, nbins+1);

binend = find(data(:,1) >= binedges(1));
if isempty(binend)
    result = NaN(nbins, size(result,2));
else
    % zero instead of NaN for backwards compatibility w/ existing code
    result = zeros(nbins, size(result,2)); 
    binend = binend(1);
    for binnum = 1:nbins
        binstart = binend + 1;
        if binstart > size(data,1)
            break
        end
        binend = binend + 1;
        while binend < size(data,1) && data(binend,1) < binedges(binnum+1)
            binend = binend + 1;
        end
        if data(binend,1) >= binedges(binnum)
            binend = binend - 1;
        end
        if yfnum
            result(binnum, 1) = mean(data(binstart:binend, 2));
            if medianflag
                result(binnum, mediancol) = median(data(binstart:binend, 2));
            end
            if countflag
                result(binnum, countcol) = sum(~isnan( ...
                    data(binstart:binend, 2) ));
            end
            if stdflag
                result(binnum, stdcol) = std(data(binstart:binend, 2));
            end
            if semflag
                result(binnum, semcol) = ...
                    std(data(binstart:binend, 2)) ...
                    / sqrt(binend - binstart + 1);
            end
        else
            result(binnum, 1) = binend - binstart + 1;
        end
    end
end

    
        
        
