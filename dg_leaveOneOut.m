function [hitrate, pred, grpidx, plevel] = ...
    dg_leaveOneOut(data, group, shuffle)
%[hitrate, pred, grpidx, plevel] = dg_leaveOneOut(data, group, shuffle)
% "Leave-one-out" validation of linear classifier (aka "decoder") as
% described in:
%   Quiroga RQ and Panzeri S, 2009. Extracting information from neuronal
%   populations: information theory and decoding approaches.  Nature Rev.
%   Neuroscience 10: 173-185.
%INPUTS
% data: observations X variables, as for <sample> and <training> in
%   Matlab's 'classify' function.
% group: data categories, as for Matlab 'classify' function.
% shuffle: optional argument.  If 0 or missing, no shuffle.  If 1, shuffle
%   all in one batch.  If >1, shuffle in batches of <shuffle>.  <group> is
%   the variable that gets shuffled (or not).  The last batch cannot be
%   guaranteed to be of length <shuffle>, but it is guaranteed to be
%   between <shuffle>/2 and 1.5*<shuffle> inclusive.
% plevel: the analytic p level for the number of correct predictions,
%   relative to the null hypothesis of random guessing.
%OUTPUTS
% hitrate: fraction of classes that are correct.
% pred: the value from <group> predicted by the classifier from each row
%   of <data>.  Suitable for use as input to dg_confusionMatrix.
% grpidx: column vector containing the index that specified the shuffling.
%   This is simply the series of counting numbers 1:length(group) when no
%   shuffling is done, or some randomized order of them when shuffling is
%   done. If all(diff(grpidx)==1), then no shuffling was done.
%NOTES
% If the decoder works at all, it is likely to work at an astonishingly
% small p-level of significance, like 1e-100.  Decoders that return results
% that are only significant at plevel = .01 or so are performing badly, so
% even though they are probably not just randomly picking groups to assign
% <data> rows to, they might also not be assigning the groups on the basis
% of meaningful physiological differences.

%$Rev: 187 $
%$Date: 2014-01-13 16:54:05 -0500 (Mon, 13 Jan 2014) $
%$Author: dgibson $

if nargin < 3 
    shuffle = 0;
end

if size(group,2) ~=1
    error('dg_leaveOneOut:group', ...
        '<group> must be a column vector.');
end

switch shuffle
    case 0
        grpidx = (1:length(group))';
    case 1
        [~, grpidx] = sort(rand(size(group)));
    otherwise
        numbatches = floor(length(group) / shuffle);
        remainder = length(group) - numbatches * shuffle;
        if remainder < length(group) / 2
            numbatches = numbatches + 1;
        end
        grpidx = NaN(size(group));
        for batchnum = 1:numbatches
            batchidx = (batchnum - 1) * shuffle + 1 ...
                : min(batchnum * shuffle, length(group));
            [~, grpidx(batchidx)] = sort(rand(size(batchidx)));
            grpidx(batchidx) = grpidx(batchidx) + (batchnum - 1) * shuffle;
        end
end

if iscell(group)
    pred = cell(size(group));
elseif isnumeric(group)
    pred = zeros(size(group));
else
    pred = repmat(group(1), size(group));
end

numhits = 0;
for k = 1:length(group)
    pred(k) = classify(data(k, :), data([1:k-1 k+1:end], :), ...
        group(grpidx([1:k-1 k+1:end])), 'linear');
    if isequal(pred(k), group(grpidx(k)))
        numhits = numhits + 1;
    end
end
hitrate = numhits / length(group);
plevel = sum(binopdf(numhits:length(group), ...
    length(group), 1/length(unique(group)) ));

