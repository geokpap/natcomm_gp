function dg_findDupRoots(infile)
%INPUTS
% infile: path to a Pouzznerian duplicated file listing.
%NOTES
% This function is driven by a monstrous input file, each line of which
% gives information for a duplicated file, which is identified by its
% absolute path.  For reasons of modularity, the top-level loop simply
% calls a function with the current and previous lines of the file as
% arguments.  For ease of testing and debugging, that function resides in a
% separate file of its own.  To avoid duplicating a huge data structure on
% every iteration, that structure is a global variable.
% Strategy:
%     1.	Iterate line-by-line through the Pouzznerian dups file
%     <infile>.
%     2.	construct a series of parallel directory trees.  Each node has
%     subdirectories, and parallel directories (nodes) that are duplicates.
%     Subdirs all go in a list that recursively uses the same data
%     structure. Duplicates are pointed to as indices to preceding
%     records.  Each record contains only one such pointer, so cases of
%     multiple duplication may require the traversal of multiple records to
%     get the entire list of duplicates.
%     3.	recursively calculate the fraction of files duplicated in each
%     directory's sub-tree.
%     4.	Recursively find the root dir of each duplicated tree:  compare
%     the fraction of files duplicated in each directory to the fraction
%     duplicated in each of its copies.  If all copies are 100% duplicated,
%     then they are all equivalent, but I prefer the shortest pathname;
%     mark the longer ones for deletion.  If there is one copy that is not
%     100% duplicated, it must be the original; mark the others for
%     deletion.  If there is more than one copy that is not 100%
%     duplicated, panic.
%     5.	Produce a report for each nominally original root dir giving
%     its pathname and the duplicates' pathnames that were marked for
%     deletion.
% This function is specifically written for a Unix environment, since the
% monster file in question is a list of files on chunky.  I've tried to
% make it as OS-independent (e.g. using <filesep> instead of the literal
% character '/'), but there is no guarantee that it will work correctly
% e.g. under Windows.

%$Rev: 263 $
%$Date: 2018-07-30 17:53:28 -0400 (Mon, 30 Jul 2018) $
%$Author: dgibson $

% <R> is a tree structure representing the directory tree of duplicated
% files.  It is implemented as a flat array, where R(1) is the root
% directory, and each node is another struct containing the same fields as
% <R(1)>.  It is recursive in
% the sense that any element may contain indices that refer to other
% elements.  Each element of <R> represents a directory.  Ordinary files
% are listed in the 'files' field of the directory that contains them.
% 'subdiridx' is not logically necessary, but reduces the amount of
% recurses foiled again.
% The fields are:
%   path: the string representation of the path to the current directory.
%   subdirs: all subdirectories of the current directory, represented
%       as absolute pathname strings.  Not sorted. Column vector.
%   subdiridx: same information as in 'subdirs' field, but presented as
%       indices into <R>.
%   files: filenames of ordinary (non-directory) files in the current
%       directory. Column vector.
%   dupidx: into <R> for the (duplicate of the current) file that is listed
%       on the preceding line of <infile>.  Scalar.
%   fracdup: to be calculated later; here it just serves as a memorandum.
%   delete: to be calculated later; here it just serves as a memorandum.

global R

% In the beginning, there is nothing (except an empty record denoting the
% absolute root of all pathnames).  R(1) is by definition the root dir.
R = [];
R.path = filesep;
R.depth = 1;
R.subdirs = {};
R.subdiridx = [];
R.files = {};
R.dupidx = [];
R.fracdup = [];
R.delete = [];

% Steps 1 and 2:
fid = fopen(infile);
if fid == -1
    error('dg_findDupRoots:file', ...
        'Could not open %s.', infile);
end
linenum = 1;
tic;
try
    prevline = '';
    thisline = fgetl(fid);
    while ~isequal(thisline, -1)
        if mod(linenum, 1000) == 0
            fprintf('%d seconds: read line #%d\n', round(toc), linenum);
        end
        dg_findDupRoots_processline(thisline, prevline, linenum);
        prevline = thisline;
        thisline = fgetl(fid);
        linenum = linenum + 1;
    end
catch e
    fclose(fid);
    fprintf('While processing line %d\n', linenum);
    rethrow(e);
end

% Step 3:
dg_findDupRoots_calcfracdup(1);
fprintf('%d seconds: finished dg_findDupRoots_calcfracdup\n', round(toc));

% Step 4:

% Step 5:

disp('Still more code to write...');

    