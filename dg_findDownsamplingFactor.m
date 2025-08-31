function N = dg_findDownsamplingFactor(sessiondir, requiredtrodes, varargin)
%N = dg_findDownsamplingFactor(sessiondir, requiredtrodes)
% If there is a 'dg_downsampleIfNeeded.mat' file in <sessiondir>, the value
% of N is read from it and <requiredtrodes> is ignored.  If there is no
% 'dg_downsampleIfNeeded.mat' file, dg_findDownsamplingFactor finds the
% maximum downsampling factor <N> for which there is a full set of files in
% <sessiondir> for each electrode in <requiredtrodes>.  It is assumed that
% the downsampling factor is coded into the CSC filenames to
% case-insensitively match the regexp '^csc(\d+)_down(\d+).mat$'.
%INPUTS
% sessiondir: absolute or relative path to directory containing downsampled
%   CSC files.
% requiredtrodes: numeric array of electrode numbers that must be
%   represented among a set of downsampled files.
%OUTPUTS
% N: the maximum common downsampling factor, or NaN if there is no such
%   downsampling factor.
%OPTIONS
% 'force' - do the full analysis even if there is a
%   'dg_downsampleIfNeeded.mat' file.

forceflag = false;

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'force'
            forceflag = true;
        otherwise
            error('dg_findDownsamplingFactor:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

N = NaN;
matfilepath = fullfile(sessiondir, 'dg_downsampleIfNeeded.mat');
if ~forceflag && exist(matfilepath, 'file')
    s = load(matfilepath, '-mat');
    N = s.N;
    fprintf('dg_downsampleIfNeeded.mat: N = %d\n', N);
else
    % Find all the different downsampling factors, count the number of
    % trodes we have for each one, and determine the Chosen Downsampling
    % Factor <N> (if there is one).  Note that this has become an
    % unnecessarily complicated method now that I don't attempt to "finish
    % up" partially downsampled sessions, but I left the old code in case I
    % decide to go back.
    srcfiles = dir(fullfile(sessiondir, 'csc*_down*.MAT'));
    srcfiles = [srcfiles
        dir(fullfile(sessiondir, 'csc*_down*.mat'))];
    srcfilenames = {srcfiles.name};
    toks = regexpi(srcfilenames, '^csc(\d+)_down(\d+).MAT$', 'tokens');
    trodenums = NaN(length(toks), 2);
    for fnidx = 1:length(srcfilenames)
        if ~isempty(toks{fnidx})
            trodenums(fnidx,:) = cellfun(@str2double, toks{fnidx}{1});
        end
    end
    % delete the NaN rows which represent filenames that didn't match the
    % regexp (e.g. files that are already local-avg-ref-subtracted):
    trodenums(isnan(trodenums(:,1)), :) = [];
    Nvals = unique(trodenums(:,2));
    Nvals(isnan(Nvals)) = [];
    numtrodes = zeros(size(Nvals));
    isgoodNval = false(size(Nvals));
    for Nvalidx = 1:length(Nvals)
        numtrodes(Nvalidx) = sum(ismember(trodenums(:,1), requiredtrodes) ...
            & trodenums(:,2)==Nvals(Nvalidx));
        isgoodNval(Nvalidx) = all(ismember( requiredtrodes, ...
            trodenums(trodenums(:,2)==Nvals(Nvalidx), 1) ));
    end
    if sum(isgoodNval) > 1
        N = max(Nvals(isgoodNval));
    elseif sum(isgoodNval) < 1
        warning('dg_findDownsamplingFactor:reqmissing', ...
            'In dir %s:\nno N has all of electrodes %s', ...
            sessiondir, dg_canonicalSeries(requiredtrodes));
    else
        N = Nvals(isgoodNval);
    end
    fprintf('dg_findDownsamplingFactor: N = %d\n', N);
end
