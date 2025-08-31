function suffix = dg_findRealSessionFiles(sessiondir)
% Uses file time stamps to identify the set of files that constitute the
% "real" session data, which I take to be the files that are separated by
% the longest time from the preceding batch of files, provided that the
% longest temporal separation is more than one hour.  If it isn't, then the
% earliest batch of files is used.
%INPUT
% sessiondir: absolute or relative pathname to a Neuralynx session
%  directory.
%OUTPUT
% suffix: string to be inserted in filenames like "SE63<suffix>.nse",
%   "CSC5<suffix>.ncs", "Events<suffix>.nev".  Presumably this also applies
%   to "VT*.mp4", "VT*.smi", and "VT*.nvt", but they all seem to get saved
%   on a different schedule during the session, so it's hard to tell.

if ~exist(sessiondir, 'dir')
    error('dg_findRealSessionFiles:sessiondir', ...
        'No such directory: %s', sessiondir);
end
nlxfiles = dir(fullfile(sessiondir, '*.n*'));
if length(nlxfiles) < 3
    error('dg_findRealSessionFiles:nonsession', ...
        '%s is not a Neuralynx session directory.', sessiondir);
end

% Get listing of non-empty ordinary (non-directory) files:
files = dir(sessiondir);
files([files.isdir]) = [];
files([files.bytes] <= 16384) = [];
% <ts> is a datenum, which is in units of days; <onehour> in days:
onehour = 1/24;
ts = [files.datenum];
% The session time will be within 6 hours either way of the median <ts>.
% Histogram all the timestamps in 15-minute bins:
medts = median(ts);
binsize = onehour / 4;
% number of bins on each side of <medts> should span 6 hours, which is half
% a day:
nbins = 0.5 / binsize;
binedges = medts + (-nbins:nbins) * binsize;
counts = histc(ts, binedges); %#ok<HISTC>
% Eliminate small counts:
counts(counts < max(counts)/2) = 0;
% Eliminate counts that represent side-bins of a larger peak:
for binnum = 1:nbins
    if counts(binnum) > 0
        if counts(binnum+1) > counts(binnum)
            counts(binnum) = 0;
        elseif counts(binnum-1) > counts(binnum)
            counts(binnum) = 0;
        end
    end
end
% We now have one non-zero count for each set of files, indicating the
% middle of the timestamp range for that set of files.  Find the timestamp
% of the "real" session files, <realts>:
sessionts = binedges(counts > 0);
if length(sessionts) == 1
    realts = sessionts(1);
else
    sessionIEI = diff(sessionts);
    [~, sortidx] = sort(sessionIEI);
    if sessionIEI(sortidx(end)) > onehour
        % <sortidx> indexes sessionIEI, and so <sortidx+1> indexes the set
        % of session files that were saved at the end of
        % <sessionIEI(sortidx)>:
        realts = sessionts(sortidx(end)+1);
    else
        realts = sessionts(1);
    end
end
% Look at the filenames that were saved between binedges(realsessidx) and
% binedges(realsessidx+1) to divine what suffix, if any, to use:
filenames = {files(ts >= realts & ts < realts + binsize).name};
toks = regexp(filenames, '^.*(_\d+)\..*$', 'tokens');
foundmatch = cellfun(@(a) ~isempty(a), toks);
suffices = cell(size(filenames));
suffices(foundmatch) = cellfun(@(a) a{1}, toks(foundmatch));
suffices(~foundmatch) = {''};
suffix = unique(suffices);
if length(suffix) > 1
    % There are four known filenames that have a different schedule of
    % getting saved ('CheetahLogFile.txt', 'DataProcessingErrors.nde',
    % 'VT2*.mp4', 'VT2*.smi'), so we can safely ignore any suffices that
    % occur fewer than 5 times.  Some sessions also have other distractor
    % files (e.g. 'VT1*.*'), so we also ignore suffices that occur in less
    % than 10% of the filenames.
    numtimes = zeros(size(suffix));
    for sufidx = 1:length(suffix)
        numtimes(sufidx) = sum(ismember(suffices, suffix(sufidx)));
    end
    isbadsuffix = numtimes < 5 | numtimes < max(numtimes)/9;
    suffix(isbadsuffix) = [];
    if length(suffix) > 1
        % We still have a problem.
        error('dg_findRealSessionFiles:suffix1', ...
            'More than one suffix found: %s', dg_thing2str(suffix));
    elseif isempty(suffix)
                error('dg_findRealSessionFiles:suffix2', ...
            'No suffix left after removing trivial cases.');
    end
end
suffix = suffix{1};

