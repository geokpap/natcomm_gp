function dg_filesplitClust(infile, outfiles, nrec1) %#ok<INUSD>
%dg_filesplitClust(infile, outfiles)
%dg_filesplitClust(infile, outfiles, offset, nrec1)
% Similar to dg_filesplitTT, but operates on Multiple Cut Cluster .mat
% files, with or without waveform data.
%INPUTS
%   infile: absolute or relative pathname of file to split.
%   outfile: a cell array containing two strings representing the absolute
%       or relative pathnames of the two output files.
%   nrec1: the last record number in the first output file.  The second
%       output file starts at record <nrec1> + 1.  Default value is half
%       the file.
%OUTPUTS
%   There are no return values.  All output is to the specified files.

%$Rev: 158 $
%$Date: 2012-09-19 16:57:14 -0400 (Wed, 19 Sep 2012) $
%$Author: dgibson $

S = load(infile);
fields = fieldnames(S);

if nargin < 3
    nrec1 = round(size(S.(fields{1}),1)/2); %#ok<NASGU>
end

eval(sprintf('%s = S.(fields{1})(1:nrec1,:);', fields{1}));
save(outfiles{1}, fields{1});
eval(sprintf('%s = S.(fields{1})(nrec1+1:end,:);', fields{1}));
save(outfiles{2}, fields{1});
