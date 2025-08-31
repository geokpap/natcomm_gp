function [offset, nrec1] = dg_filecatTT(infiles, outfile)
%dg_filecatTT(infiles, outfile) concatenates the contents of two
% Neuralynx TT input files to a Neuralynx TT output file.  Stupidly uses
% the first file's header as the header for the output file, does not check
% for consistency of contents.  Adds a temporal offset to all timestamps in
% the second file sufficient so that the first timestamp in the second file
% is greater than the last timestamp in the first file (it is assumed that
% the timestamps are in ascending order).
% INPUTS
%   infiles: a cell array of strings, each of which points to a file to be
%       read and copied to the output.
%   outfile: a string that points to the destination file to which output
% OUTPUTS
%   offset: the temporal offset of the second file, expressed in the
%       units of the timestamps in the input files.
%   nrec1: the number of records in the first file.
% NOTES
% As of 2013-11-05, the following usually useless fields are no longer
%   copied to the output: ScNumbers, CellNumbers, Params.

%$Rev: 183 $
%$Date: 2013-11-05 14:15:48 -0500 (Tue, 05 Nov 2013) $
%$Author: dgibson $

for k = 1:2
    [spikeTS{k}, DataPoints{k}, NlxHeader{k}] = dg_readSpike(infiles{k}); %#ok<AGROW>
end

nrec1 = size(DataPoints{1},3);
maxtimestamp = spikeTS{1}(end);
timestamp_exp = fix(log10(maxtimestamp));
offset = (fix(maxtimestamp/(10^timestamp_exp)) + 2) * 10^timestamp_exp;
spikeTS{2} = spikeTS{2} + offset;

dg_writeSpike(outfile, [spikeTS{1} spikeTS{2}], ...
    cat(3, DataPoints{1}, DataPoints{2}), NlxHeader{1});
