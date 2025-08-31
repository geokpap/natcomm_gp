function [results, filenames] = dg_iterateOverFiles(sessiondir, funch, ...
    matchExpr)
%[results, filenames] = dg_iterateOverFiles(sessiondir, funch, matchExpr)
% Possibly useful to those who write the occasional bulk process.
%INPUTS
% sessiondir: directory containing files to process
% funch: function handle to function that accepts a single absolute
%   pathname as an argument and returns a single value (use a struct if you
%   need to return more than just one array).
% matchExpr: regexp that defines which files to use, or a cell array of
%   strings that contain literal filenames relative to <sessiondir>.
%OUTPUT
% results: a cell column vector containing the result returned by <funch>
%   after processing each file.
% filenames: a cell string array of the same size as <results> specifying
%   which file produced each result.

%$Rev: 218 $
%$Date: 2015-04-15 13:48:13 -0400 (Wed, 15 Apr 2015) $
%$Author: dgibson $

results = {0,1}; %#ok<NASGU>
if iscell(matchExpr)
    results = cell(numel(matchExpr), 1);
    filenames = reshape(matchExpr, [], 1);
else
    files = dir(sessiondir);
    filenames = reshape({files(~cell2mat({files.isdir})).name}, [], 1);
    fnidx = find(~cellfun(@isempty, regexpi(filenames, matchExpr)));
    if isempty(fnidx)
        error('dg_iterateOverCSCs:nomatch', ...
            'There is no match in %s for ''%s''', ...
            sessiondir, matchExpr);
    end
    filenames = filenames(fnidx);
    results = cell(numel(filenames), 1);
end
for fidx = 1:length(results)
    results{fidx} = feval(funch, fullfile(sessiondir, filenames{fidx}));
end

