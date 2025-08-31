function dg_filesplitTT(infile, outfiles, offset, nrec1)
%dg_filesplitTT(infile, outfiles)
%dg_filesplitTT(infile, outfiles, offset, nrec1)
% The inverse of dg_filecatTT, except of course the second file just gets
% another copy of the same header.
%INPUTS
%   infile: absolute or relative pathname of file to split.
%   outfile: a cell array containing two strings representing the absolute
%       or relative pathnames of the two output files.
%   offset: this value is subtracted from all the timestamps in the second
%       output file; must be specified in same units as the timestamps in
%       <infile> (typically microseconds).  Default value is 0.
%   nrec1: the last record number in the first output file.  The second
%       output file starts at record <nrec1> + 1.  Default value is half
%       the file.
%OUTPUTS
%   There are no return values.  All output is to the specified files.

%$Rev: 207 $
%$Date: 2014-10-16 19:07:56 -0400 (Thu, 16 Oct 2014) $
%$Author: dgibson $

if nargin < 3
    offset = 0;
end

[TimeStamps, ScNumbers, CellNumbers, Params, DataPoints, NlxHeader] = ...
    Nlx2MatTT(infile,1,1,1,1,1,1);

if nargin < 4
    nrec1 = round(length(TimeStamps)/2);
end

if ~ispc
    error('dg_filesplitTT:notPC', ...
        'Neuralynx format files can only be created on Windows machines.');
end
Mat2NlxTT(outfiles{1}, 0, 1, 1, nrec1, [1 1 1 1 1 1], ...
    TimeStamps(1:nrec1), ScNumbers(1:nrec1), CellNumbers(1:nrec1), ...
    Params(:, 1:nrec1), DataPoints(:, :, 1:nrec1), NlxHeader );

if ~ispc
    error('dg_filesplitTT:notPC', ...
        'Neuralynx format files can only be created on Windows machines.');
end
Mat2NlxTT(outfiles{2}, 0, 1, 1, size(DataPoints,3) - nrec1, [1 1 1 1 1 1], ...
    TimeStamps(nrec1+1:end) - offset, ScNumbers(nrec1+1:end), ...
    CellNumbers(nrec1+1:end), Params(:, nrec1+1:end), ...
    DataPoints(:, :, nrec1+1:end), NlxHeader );
