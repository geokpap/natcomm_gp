function result = dg_decisionMatrix(trialparams, values, varargin)
% Computes a decision matrix result as in Amemori, Amemori, and Graybiel
% 2015 J Neurosci.
%INPUTS
% trialparams: one row per trial; col 1: reward magnitude, col 2: aversion
%   magnitude.
% values: one row per trial, one column; any number associated with each
%   trial.
%OUTPUTS
% result: 8 x 8 matrix containing <values> binned by <trialparams> and
%   smoothed as in Amemori et al., with reward across columns and aversion
%   down rows.
%OPTIONS
% 'range', range - <range> is a two-element vector containing the minimum
%   and maximum values to use for binning <trialparams>.  Both columns of
%   <trialparams> are binned on the same scale.  If the first element of
%   <range> is not zero, this may cause unexpected results.
%NOTES
% The default range of <trialparams> is 0:199 as in the Georgios data.  The
% bin index is calculated as idx = floor(trialparams(:,k) / binwidth) + 1,
% which means that a value that falls exactly on a bin edge (e.g. 25 for
% the default values) goes into the higher bin (e.g. bin 2).  I.e. binning
% is done based on binedge(n) <= value & value < binedge(n+1).
%   The binning is done by converting directly to bin indices, based on the
% assumption that <trialparams> is integer-valued and has a minimum value
% of 0. This makes this function rather application-specific (hence its
% name).

%$Rev: 305 $
%$Date: 2023-09-08 14:19:01 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

if size(values, 1) ~= size(trialparams, 1)
    error('dg_decisionMatrix:rows', ...
        '<trialparams> and <values> must have the same number of rows.');
end
minrwd = 0;
maxrwd = 199;
minave = 0;
maxave = 199;
numbins = 8;

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'range' % does not get deleted
            argnum = argnum + 1;
            minrwd = varargin{argnum}(1);
            maxrwd = varargin{argnum}(2);
            minave = varargin{argnum}(1);
            maxave = varargin{argnum}(2);
        otherwise
            error('funcname:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

% <matrix> is initialized to NaN, and the <values> element for each trial
% is entered into the appropriate bin in matrix(aveidx, rwdidx, trialidx).
% When all data have been entered, we take the nan-tolerant mean over
% trials.
matrix = NaN(8, 8, size(values, 1));
rwdbinwidth = (maxrwd + 1 - minrwd) / numbins;
avebinwidth = (maxave + 1 - minave) / numbins;
rwdidx = floor(trialparams(:,1) / rwdbinwidth) + 1;
aveidx = floor(trialparams(:,2) / avebinwidth) + 1;
if any([rwdidx; aveidx] > numbins)
    error('dg_decisionMatrix:brainfog', ...
        'This can''t happen.');
end
for trialidx = 1:length(values)
    matrix(aveidx(trialidx), rwdidx(trialidx), trialidx) = ...
        values(trialidx);
end
result = mean(matrix, 3, 'omitnan');
