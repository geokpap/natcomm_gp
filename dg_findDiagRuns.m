function [k, runstart, runsum] = dg_findDiagRuns(M)
%[k, runstart, runsum] = dg_findDiagRuns(M)
%   For each diagonal k = -(size(M,1)-1) : (size(M,2)-1), computes the sum
%   of the elements of M over each contiguous run of nonzeros along that
%   diagonal, and returns the values of k, the index in the diagonal where
%   the run started <runstart>, and the sum over the run <runsum>, sorted
%   primarily in order of decreasing <runsum>, then in order of increasing
%   <k> and then in order of increasing <runstart>.
%   All result values are column vectors.  There will be as many instances
%   of each value of <k> as there are contiguous runs along that diagonal.
%   A "contiguous run" is by definition a minimum of two consecutive
%   elements.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

results = NaN(size(M,1) + size(M,2) - 1, 3);
% kvals = NaN(size(M,1) + size(M,2) - 1, 1);
% runstart = NaN(size(kvals));
% runsum = NaN(size(kvals));
row = 1;
for k = -(size(M,1)-1) : (size(M,2)-1)
    elements = diag(M, k);
    contig = elements(1:end-1) ~= 0 & elements(2:end) ~= 0;
    if ~any(contig)
        results(row,1) = k;
        results(row,2) = 0;
        results(row,3) = 0;
        row = row + 1;
    else
        runstartidx = find([ contig(1)
            ~contig(1:end-1) & contig(2:end) ]);
        % The first element cannot be the end of a run by definition.
        % if contig(end), then the last two elements were nonzero, which
        % means that the last run ends with length(elements).  If
        % ~contig(end), the last two elements could be [0 1], [1 0], [0 0];
        % in all three cases, the last run does NOT end with
        % length(elements).
        runendidx = find([ 0
            contig(1:end-1) & ~contig(2:end)
%             ~contig(end) && elements(end-1) ~= 0 && elements(end) == 0
            contig(end) ]);
        if length(runstartidx) ~= length(runendidx)
            error('dg_findDiagRuns:oops', 'oops');
        end
        for runidx = 1:length(runstartidx)
            results(row,1) = k;
            results(row,2) = runstartidx(runidx);
            results(row,3) = sum(elements( ...
                runstartidx(runidx) : runendidx(runidx) ));
            row = row + 1;
        end
    end
end
results = sortrows(results);
[runsum, sortidx] = sort(results(:,3), 'descend');
runstart = results(sortidx,2);
k = results(sortidx,1);

