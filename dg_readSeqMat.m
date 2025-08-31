function [seqs, positions, params, counts] = ...
    dg_readSeqMat(filename, seqvarname, CPvals, TPvals)
%DG_READSEQMAT reads Joey's sequence .MAT files.
%[seqs, positions, params, counts, TPvals] = dg_readSeqMat(filename, ...
%   seqvarname, CPvals, TPvals)
% NaNs in <positions> or <params> denote aggregated values.

% Joey's sequence .MAT files are structured as follows.
% TP:
%   The rows of <TP> contain every possible combination of the 6 relevant
%   non-position Trial Params, each of which is a binary parameter.  The
%   parameters are arranged across the columns as:
%   1. cue colorset (4-RGB, 10-PAO)
%   2. # of joystick movements required ("0" is 3 movements, "1" is 1
%       movement)
%   3. presentation type ("0" is simultaneously presented cues, "1" is
%       sequentially cued)
%   4. H1 duration ("0" is short, "1" is long)
%   5. H2 duration (ditto)
%   6. H3 duration (ditto)
% CP:
%   Cue Permutations, one per row.
% H1i, H1r, and/or some other variable name specified by <seqvarname>:
%   The rows correspond to rows of <TP>.  The columns correspond to rows of
%   <CP>.  These are sparsely populated, as most rows of TP are not
%   represented in the data set. Each element of H1i is a cell containing
%   all the spatial sequence strings in the data set, and H1r contains
%   rank-order/color sequences. 
% ID:
%   Cell array, contains the sets of unique trial IDs for the non-empty
%   cells in H1i, H1r.  Not used in this application.
%
% INPUTS
% <TPvals>, <CPvals> - these function similarly, but TPvals relates to
%   values in the <TP> array, and <CPvals> relates to values in the <CP>
%   array.  They control the selection and aggregation of trials based on
%   trial parameters and cue permutations.  This could be done in the
%   database query that produced the original sequences.  As a convenience
%   when reading from .MAT files, the same functionality is provided here.
%   <TPvals> is a row vector that has the same number of columns as the
%   <TP> variable in the .MAT file, and <CPvals> is a row vector that has
%   the same number of columns as the <CP> variable.  In each column, if
%   the value is not NaN, and not -1, then it specifies a value that must
%   be matched exactly by the value in <TP> in order for the corresponding
%   collection of sequences to be included.  If the value is NaN, then the
%   values in the corresponding column of <TP> are ignored and all the data
%   in the .MAT file are used in aggregate.  If the value is -1, then all
%   the data are used, but each value is analyzed separately.  E.g. suppose
%   the .MAT file contains trials with cues at positions 5, 7, and 8.  If
%   CPvals = [5 NaN NaN], then only trials whose first target is at
%   position 5 will be included in the sequence counts, the counts will be
%   aggregated over all second and third targets, and so every row of
%   <positions> will be [5 NaN NaN]; whereas if CPvals = [5 -1 -1], then
%   (as before) only trials whose first target is target ID 5 will be
%   included in the sequence counts, but separate tallies will be made for
%   a given sequence's occurences in trials with target positions [5 7 8]
%   versus [5 8 7].  CPvals = [NaN 5 -1] says to include only trials that
%   have the second cue at position 5, ignore the position of the first
%   cue, and separately tally sequences for the third cue being at position
%   7 or 8 (note that for any .MAT file that contains only three cue
%   positions, this is equivalent to CPvals = [-1 5 -1]).  CPvals = [-1 NaN
%   NaN] says to include all trials in the file, and aggregate trials
%   regardless of the second and third cue positions, but compute separate
%   tallies for each first target position (so then <positions> would
%   contain rows with values [5 NaN NaN], [7 NaN NaN], and [8 NaN NaN]).
%   The scalar value -1 can be used as shorthand for all -1s.  The value []
%   is equivalent to all NaNs.
%
% OUTPUTS
% <seqs>, <positions>, <params>, <counts> all have the same number of rows.
% <seqs> is a cell string column vector containing variable-length string
%   data representing fixations sequences (could be position-coded or
%   order-coded). 
% <positions> is a numeric array with the same number of rows as <seqs>, 
%   and same number of columns as CP, containing the spatial positions of
%   the three targets, with NaNs inserted in columns where different values
%   have been aggregated. 
% <params> contains the rows of <TP> corresponding to each row of
%   <seqs>, with NaNs inserted in columns where different values have been
%   aggregated. 
% <counts> contains the number of instances within the
%   aggregate of the corresponding sequence from <seqs>.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if size(TPvals,1) > 1
    error('dg_readSeqMat:badTPvals2', ...
        '<TPvals> must be a row vector.')
end
if size(CPvals,1) > 1
    error('dg_readSeqMat:badCPagg2', ...
        '<CPvals> must be a row vector.')
end
    
% Load data; for coding convenience, copy contents of variable whose name
% is the value of <seqvarname> into <rawseqs>.
load(filename);
eval(sprintf('rawseqs = %s;', seqvarname));

if isequal(TPvals, -1)
    TPvals = zeros(1, size(TP,2)) - 1;
end
if isempty(TPvals)
    TPvals = nan(1, size(TP,2));
end
if size(TPvals,2) ~= size(TP,2)
    error('dg_readSeqMat:badTPvals', ...
        '<TPvals> must have the same number of columns as TP.')
end
if isequal(CPvals, -1)
    CPvals = zeros(1, size(CP,2)) - 1;
end
if isempty(CPvals)
    CPvals = nan(1, size(CP,2));
end
if size(CPvals,2) ~= size(CP,2)
    error('dg_readSeqMat:badCPvals', ...
        '<CPvals> must have the same number of columns as CP.')
end

% find and delete empty rows:
gotnodata = all(cell2mat(dg_mapfunc(@isempty, rawseqs)), 2);
rawseqs(gotnodata,:) = [];
TP(gotnodata,:) = [];

% filter as specified by TPvals
matchcols = find(~isnan(TPvals) & (TPvals ~= -1));
TPexclude = ~all( TP(:,matchcols) == repmat( ...
    TPvals(1,matchcols), size(TP,1), 1 ), 2 );
TP(TPexclude,:) = [];
rawseqs(TPexclude,:) = [];
% filter as specified by CPvals
matchcols = find(~isnan(CPvals) & (CPvals ~= -1));
CPexclude = ~all( CP(:,matchcols) == repmat( ...
    CPvals(1,matchcols), size(CP,1), 1 ), 2 );
CP(CPexclude,:) = [];
rawseqs(:,CPexclude) = [];

% aggregate as specified
if any(isnan(TPvals))
    aggcols = find(isnan(TPvals));
    [TP, aggseqs] = dg_aggregate(TP, rawseqs, aggcols);
end
if any(isnan(CPvals))
    aggcols = find(isnan(CPvals));
    [CP, seqstmp] = dg_aggregate(CP, aggseqs', aggcols);
    aggseqs = seqstmp';
end

%
% For each non-empty cell in the aggregated sequences table, compute the
% list of unique <aggseqs> and the corresponding counts, params & positions
%

seqs = {};
positions = [];
params = [];
counts = [];
for tpidx = 1:size(aggseqs,1)
    for cpidx = 1:size(aggseqs,2)
        uniqseqs = unique(aggseqs{tpidx, cpidx});
        for seqidx = 1:length(uniqseqs)
            counts(end+1,1) = ...
                sum(ismember(aggseqs{tpidx, cpidx}, uniqseqs(seqidx)));
            seqs{end+1,1} = uniqseqs{seqidx};
            params(end+1,:) = TP(tpidx,:);
            positions(end+1,:) = CP(cpidx,:);
        end
    end
end


