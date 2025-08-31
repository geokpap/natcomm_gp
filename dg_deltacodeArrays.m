function [idx, Bvals] = dg_deltacodeArrays(A, B, outfile)
%[idx, Bvals] = dg_deltacodeArrays(A, B, outfile)
% Delta-codes the differences between two arrays for more efficient
% storage.
%INPUTS
% A, B - two arrays of the same size
% outfile - relative or absolute path to which to save output.  This is
%   passed unmodified to the 'save' function, so if there is no filename
%   extension, it will receive the system default (e.g. '.mat') extension.
%   If empty, no output file is saved.
%OUTPUTS
% idx: logical index of elements that differ between <A> and <B>
% Bvals: B(idx)
%NOTES
% To reconstruct B, run the following code:
%   B = A;
%   load(outfile, '-mat');
%   B(idx) = Bvals;

%$Rev: 170 $
%$Date: 2013-03-01 17:50:20 -0500 (Fri, 01 Mar 2013) $
%$Author: dgibson $

if ~isequal(size(A), size(B))
    error('dg_deltacodeArrays:eq', ...
        '<A> and <B> must be the same size.');
end

idx = A ~= B;
Bvals = B(idx);
if ~isempty(outfile)
    save(outfile, 'idx', 'Bvals');
end
