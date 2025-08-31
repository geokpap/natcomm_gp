function fracdup = dg_findDupRoots_calcfracdup(Ridx)
% Returns the fraction of files that are duplicated (.fracdup) in an
% arbitrary subdirectory specified by <Ridx>, an index into <R>.
% If R(Ridx) still has an empty .fracdup field, calculates .fracdup of
% subdirs as needed. Inconveniently, "duplicated" means with respect to the
% contents of the actual filesystem, not just the files that got flagged
% and entered into <R>.  Specifically, it means the fraction of files in
% the specified directory in the actual filesystem that are listed in <R>
% as being duplicates.  Subdirectories contribute their own .fracdup score
% to the final number, i.e. R(Ridx).fracdup = ( # of duplicate files +
% sum(duplicated directory .fracdups) ) / (# of files and dirs in actual
% filesystem).  Thus a directory tree can be marked as containing files
% that did not get duplicated, even when those files are deeply buried many
% subdirectories down.
%   This function has The Power to call itself recursively, which it uses.
%   Directories that have already been deleted are assigned a .dupfrac
% value of 0.  Thus they contribute nothing to the .dupfrac score being
% computed.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

global R

fracdup = NaN; %#ok<NASGU>

% The trivial case:
if ~isempty(R(Ridx).fracdup)
    fracdup = R(Ridx).fracdup;
    return
end

% Actually calculate something:
files = dir(R(Ridx).path);
if isempty(files)
    actualnum = 0;
else
    % The "- 2" is for the obligatory "hidden" directories '.' and '..':
    actualnum = length(files) - 2;
end
if isempty(R(Ridx).subdiridx)
    % No recursion required, just compute the fraction of duplicated
    % ordinary files.
    if actualnum == 0
        R(Ridx).fracdup = 0;
    else
        R(Ridx).fracdup = length(R(Ridx).files) / actualnum;
    end
else
    % Recurse (foil again :-P)
    subdirfracdups = NaN(size(R(Ridx).subdirs));
    for sdidx = 1 : length(R(Ridx).subdirs)
        subdirfracdups(sdidx) = ...
            dg_findDupRoots_calcfracdup(R(Ridx).subdiridx(sdidx));
    end
    if actualnum == 0
        R(Ridx).fracdup = 0;
    else
        R(Ridx).fracdup = (length(R(Ridx).files) + sum(subdirfracdups)) ...
            / actualnum;
    end
end

fracdup = R(Ridx).fracdup;

