function [hL, binvals, binctrs] = dg_linehisto(data, edges, hA, varargin)
%hL = dg_linehisto(data, edges, hA)
% Constructs histograms and plots them as lines on a common axes.
%INPUTS
% data:  a cell array of numerical data arrays ("datasets")
% edges:  Optional.  Edges of bins, as for Matlab's "histc" function.  If
%   empty or not given, 100 bins with first edge at minimum value in <data>
%   and last edge at (maximum value in <data>) * (1+eps).  The last bin of
%   "histc", containing only data points that are exactly equal to
%   edges(end), is not used.
% hA:  Optional.  axes into which to plot.  If empty or not given, a new
%   figure window containing one axes is created.
%OUTPUT
% hL:  a vector of handles to line objects, one for each element of <data>,
%   in the same order as the elements of <data>.
% binvals:  the count in each bin, in bins X datasets format
% binctrs:  the center point of each bin
%OPTIONS
% Any options not listed below are passed on to the Matlab "plot" function.
% 'norm' - normalizes binvals so that each column of binvals sums to 1
%   (i.e. makes it into a probability distribution)

%$Rev: 151 $
%$Date: 2012-06-18 17:06:32 -0400 (Mon, 18 Jun 2012) $
%$Author: dgibson $

if nargin < 3 || isempty(hA)
    hF = figure;
    hA = axes('Parent',hF);
end
if nargin < 2
    edges = [];
end

opts2delete = [];
normflag = false;

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'norm'
                normflag = true;
                opts2delete(end+1) = argnum; %#ok<*AGROW>
        end
    end
    argnum = argnum + 1;
end
varargin(opts2delete) = [];

if isempty(edges)
    maxdata = -Inf;
    mindata = Inf;
    for k = 1:numel(data)
        newmax = max(data{k}(:));
        if newmax > maxdata
            maxdata = newmax;
        end
        newmin = min(data{k}(:));
        if newmin < mindata
            mindata = newmin;
        end
    end
    edges = linspace(mindata, maxdata*(1+eps), 101);
end

counts = NaN(length(edges), numel(data));
for k = 1:numel(data)
    counts(:,k) = histc(data{k}(:), edges);
    if normflag
        counts(:,k) = counts(:,k) / sum(counts(1:end-1,k));
    end
    
end
binvals = counts(1:end-1,:);
binctrs = (edges(1:end-1) + edges(2:end)) / 2;
hL = plot(hA, binctrs, binvals, varargin{:});

