function [minidx, mindist] = dg_seqCentroid(seqs, match, indel)
%DG_SEQCENTROID finds the centermost member of a set of sequences.

%   <seqs> is a cell vector of arrays of nonnegative integers; each cell
%   represents one sequence.  For each sequence, calculate the total edit
%   distance to all other sequences in <seqs>, and return the index to the sequence
%   with the smallest total distance together with the total distance.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

sumdist = zeros(1,length(seqs));

% Construct distance table in upper half-matrix excluding diagonal
dist = repmat(NaN, length(seqs), length(seqs));
for i1 = 1:length(seqs)
    for i2 = i1+1 : length(seqs)
        dist(i1,i2) = ...
            dg_edDist(seqs{i1}, seqs{i2}, match, indel);
    end
end

for seqidx = 1:length(seqs)
    others = 1:length(seqs);
    others(seqidx) = [];
    for seqidx2 = others
        if seqidx < seqidx2
            i1 = seqidx;
            i2 = seqidx2;
        else
            i1 = seqidx2;
            i2 = seqidx;
        end
        sumdist(seqidx) = sumdist(seqidx) + dist(i1,i2);
    end
end

[mindist, minidx] = min(sumdist);
