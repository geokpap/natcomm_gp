function dg_downsampleIfNeeded(srcroots, destroot, sessionID, ...
    allowedmissing)
%dg_downsampleIfNeeded(srcroots, destroot, sessionID)
% First, finds <srcroot>, which is a member of <srcroots> and contains
% <sessionID>. Compares the list of .ncs files in <srcroot> with .mat files
% in <destroot>, and if any aside from the ones listed in <allowedmissing>
% are missing, (re)runs dg_downsampleAndConvert( ... 'reconcile', 'all').
% If not, reads lists of timestamps from the existing downsampled files,
% and if there are mismatched timestamps, (re)runs dg_downsampleAndConvert(
% ... 'reconcile', 'all').  Before (re)running dg_downsampleAndConvert,
% moves existing downsampled session directory to subdirectory 'old' of
% <destroot>.
%   When comparing file lists, the downsampling factor is extracted from
% the filenames of the form csc<CH>_down<N>.MAT, where <CH> is the channel
% number and <N> is the downsampling factor.  If there is more than one
% downsampling factor, then there must be at least half of the .ncs files
% all sharing the SAME downsampling factor in order to escape (re)running
% dg_downsampleAndConvert.  If there is more than one downsampling factor
% that has more than half of the .ncs files downsampled, then the largest
% one is chosen. If no downsampling factor can be identified, the default
% value of 32 is used.
%INPUTS
% srcroots: cell array of directories containing raw session directories.
% destroot: directory containing downsampled session directories.
% sessionID: exactly as listed in <srcroot>.
% allowedmissing: a vector of integers representing CSC files that can be
%   missing in the downsampled output directory and still consider that
%   directory to be complete.
%OUTPUTS
% Downsampled files in directory fullfile(destroot, sessionID); plus the
% file 'dg_downsampleIfNeeded.mat', containing the downsampling factor <N>
% for the complete set of downsampled files, and the reconciled timestamp
% list <recTS>.  If <recTS> is empty, that indicates that the whole session
% was downsampled de novo using dg_downsampleAndConvert.

%$Rev: 206 $
%$Date: 2014-10-14 13:19:55 -0400 (Tue, 14 Oct 2014) $
%$Author: dgibson $

srcroot = '';
for k = 1:length(srcroots)
    if exist(fullfile(srcroots{k}, sessionID), 'dir')
        srcroot = srcroots{k};
    end
end
if isempty(srcroot)
    return
end
srcfiles = dir(fullfile(srcroot, sessionID, '*.ncs'));
srcfilenames = {srcfiles.name};
toks = regexpi(srcfilenames, '^csc(\d+).ncs$', 'tokens');
srcnums = NaN(length(toks), 1);
for fnidx = 1:length(srcfilenames)
    if ~isempty(toks{fnidx})
        srcnums(fnidx) = str2double(toks{fnidx}{1});
    end
end
trodenums = srcnums(~isnan(srcnums));
requiredtrodes = setdiff(trodenums, allowedmissing);

destfiles = dir(fullfile(destroot, sessionID, 'csc*_down*.MAT'));
destfiles = [destfiles
    dir(fullfile(destroot, sessionID, 'csc*_down*.mat'))];
destfilenames = {destfiles.name};
toks = regexpi(destfilenames, '^csc(\d+)_down(\d+).MAT$', 'tokens');
destnums = NaN(length(toks), 2);
for fnidx = 1:length(destfilenames)
    if ~isempty(toks{fnidx})
        destnums(fnidx,:) = cellfun(@str2double, toks{fnidx}{1});
    end
end

% Find all the different downsampling factors, count the number of
% trodes we have for each one, and determine the Chosen Downsampling Factor
% <N> (if there is one).  Note that this has become an unnecessarily
% complicated method now that I don't attempt to "finish up" partially
% downsampled sessions, but I left the old code in case I decide to go
% back.
Nvals = unique(destnums(:,2));
Nvals(isnan(Nvals)) = [];
numtrodes = zeros(size(Nvals));
for Nvalidx = 1:length(Nvals)
    numtrodes(Nvalidx) = sum(ismember(destnums(:,1), requiredtrodes) ...
        & destnums(:,2)==Nvals(Nvalidx));
end
isgoodNval = numtrodes == length(requiredtrodes);
if sum(isgoodNval) > 1
    N = max(Nvals(isgoodNval));
elseif sum(isgoodNval) < 1
    N = NaN;
else
    N = Nvals(isgoodNval);
end
fprintf('dg_downsampleIfNeeded: N = %d\n', N);

% If we might skip the full dg_downsampleAndConvert, check for timestamp
% conflicts:
recTS = [];
foundproblem = false;
if ~isnan(N)
    for fnidx = reshape(find( ismember(destnums(:,1), requiredtrodes) ...
            & destnums(:,2)==N ), 1, [])
        try
            S = load(fullfile(destroot, sessionID, destfilenames{fnidx}), ...
                '-mat', 'dg_Nlx2Mat_Timestamps');
        catch e
            logmsg = sprintf( 'Load failed for %s\n', ...
                fullfile(destroot, sessionID, destfilenames{fnidx}) );
            logmsg = sprintf('%s\n%s\n%s', ...
                logmsg, e.identifier, e.message);
            for stackframe = 1:length(e.stack)
                logmsg = sprintf('%s\n%s\nline %d', ...
                    logmsg, e.stack(stackframe).file, ...
                    e.stack(stackframe).line);
            end
            disp(logmsg);
            foundproblem = true;
            break
        end
        if isempty(recTS)
            recTS = S.dg_Nlx2Mat_Timestamps;
        else
            foundproblem = ~isequal(S.dg_Nlx2Mat_Timestamps, recTS);
            if foundproblem
                break
            end
        end
    end
end
fprintf('dg_downsampleIfNeeded: foundproblem = %d\n', foundproblem);

if isnan(N) || foundproblem
    if exist(fullfile(destroot, sessionID), 'dir')
        olddir = fullfile(destroot, 'old');
        if ~exist(olddir, 'dir')
            if exist(olddir, 'file')
                error('dg_downsampleIfNeeded:olddir', ...
                    'There is a file named %s, cannot create backup directory', ...
                    olddir );
            end
            mkdir(olddir);
        end
        if exist(fullfile(olddir, sessionID), 'dir')
            % ...then a previous mv got stopped in the middle, and we must
            % move the remaining files one by one.
            files = dir(fullfile(destroot, sessionID));
            filenames = {files.name};
            filenames(ismember(filenames, {'.', '..'})) = [];
            for fidx = 1:length(filenames)
                movefile(fullfile(destroot, sessionID, filenames{fidx}), ...
                    fullfile(olddir, sessionID, filenames{fidx}));
            end
            rmdir(fullfile(destroot, sessionID));
        else
            movefile(fullfile(destroot, sessionID), ...
                fullfile(olddir, sessionID));
        end
    end
    N = 32;
    fprintf('dg_downsampleIfNeeded: running dg_downsampleAndConvert\n');
    dg_downsampleAndConvert( fullfile(srcroot, sessionID), ...
        fullfile(destroot, sessionID), N, 'reconcile', 'all', ...
        'verbose', 'suffix', 'overwrite', ...
        'allowedmissing', allowedmissing );
end
save(fullfile(destroot, sessionID, 'dg_downsampleIfNeeded.mat'), ...
    'N', 'recTS');

