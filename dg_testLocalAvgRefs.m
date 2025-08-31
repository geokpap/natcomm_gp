function [results, filenames, refnames] = dg_testLocalAvgRefs( ...
    localavgrefs, sessiondir, prefix, testmode)
% For each local average reference group defined in <localavgrefs>,
% computes the cross-covariance of each raw file that went into the local
% average reference with its own local average reference.
%INPUTS
% localavgrefs: as for dg_makeLocalAvgRefs.
% sessiondir:  as for dg_makeLocalAvgRefs.
% prefix:  as for dg_makeLocalAvgRefs.
% testmode: if <true>, then computes cross covariance against all *other*
%   local average files in the session, and does not use the correctly
%   matched local average file.  Test mode will *not* work correctly on
%   values of <localavgrefs> that specify expressions for matching to file
%   names rather than literal lists of filenames.
%OUTPUTS
% results: column vector of 0-lag xcov results for every pair of files
%   tested.
% filenames: names of the raw CSC files (column vector).
% refnames: names of the local average reference files (column vector).
%NOTES
% There will usually be many more rows of output when <testmode> is <true>.
%   A major drawback of this method of assessing whether a given channel
% "belongs" in a group is that the channel does contribute to the average,
% so the xcov will never hit zero, and indeed the lower limit could be
% quite high when there are very few channels involved (I think it might be
% 1/sqrt(N) where N is the number of channels).
%   See also dg_xcovMatrix, which takes longer to run but produces more
% useful results.

%$Rev: 211 $
%$Date: 2015-03-12 21:08:04 -0400 (Thu, 12 Mar 2015) $
%$Author: dgibson $

results = [];
filenames = {};
refnames = {};

% Iterate over local average refs:
for refidx = 1:length(localavgrefs)
    CSCfilespec = localavgrefs{refidx};
    if isempty(CSCfilespec)
        continue
    end
    if ischar(CSCfilespec) && testmode
        error('dg_testLocalAvgRefs:testmode', ...
            'Test mode does not work on filename-matching expressions.');
    elseif testmode
        % construct list of all files that go into the other local average
        % groups:
        CSCfilespec = {};
        for nonrefidx = 1:length(localavgrefs)
            if nonrefidx ~= refidx && ...
                    ~isempty(localavgrefs{nonrefidx})
                CSCfilespec = [ CSCfilespec
                    reshape(localavgrefs{nonrefidx}, [], 1) ]; %#ok<AGROW>
            end
        end
    end
    refname = sprintf('%s%d.mat', prefix, refidx);
    localavgfile = fullfile(sessiondir, refname);
    [res, fnames] = testOneLocalAvg(sessiondir, CSCfilespec, ...
        localavgfile);
    results = [ results
        res ]; %#ok<AGROW>
    filenames = [ filenames
        fnames ]; %#ok<AGROW>
    refnames = [ refnames
        repmat({refname}, size(fnames)) ]; %#ok<AGROW>
end

end

function [results, filenames] = testOneLocalAvg(sessiondir, CSCfilespec, ...
    localavgfile)
% Computes the xcovs for every file in <CSCfilespec> against
% <localavgfile>.  <results> contains one value for each file in
% <CSCfilespec>.
if ischar(CSCfilespec)
    filelist = dir(fullfile(sessiondir, CSCfilespec));
    filenames = reshape({filelist.name}, [], 1);
elseif iscell(CSCfilespec)
    filenames = reshape(CSCfilespec, [], 1);
else
    error('dg_testLocalAvgRefs:CSCfilespec', ...
        'Each element of <localavgrefs> must be a cell array or a string.');
end
results = NaN(size(filenames));
[~, ~, ext] = fileparts(localavgfile);
switch ext
    case {'.mat' '.MAT'}
        x = load(localavgfile);
        localavg = x.dg_Nlx2Mat_Samples;
    case {'.ncs' '.NCS' '.dat' '.DAT'}
        [~, localavg] = dg_readCSC(localavgfile);
    otherwise
        error('dg_testLocalAvgRefs:ext', ...
            'Unknown filename extension: %s', ext);
end
for fidx = 1:length(filenames)
    onechanfile = fullfile(sessiondir, filenames{fidx});
    [~, ~, ext] = fileparts(onechanfile);
    switch ext
        case {'.mat' '.MAT'}
            x = load(onechanfile);
            onechan = x.dg_Nlx2Mat_Samples;
        case {'.ncs' '.NCS' '.dat' '.DAT'}
            [~, onechan] = dg_readCSC(onechanfile);
        otherwise
            error('dg_testLocalAvgRefs:ext2', ...
                'Unknown filename extension: %s', ext);
    end
    results(fidx) = xcov(localavg(:), onechan(:), 0, 'coeff');
end

end


