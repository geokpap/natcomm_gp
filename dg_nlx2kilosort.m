function channelnums = dg_nlx2kilosort(sessiondir, outfile)
%channelnums = dg_nlx2kilosort(sessiondir, outfile)
% Creates a file suitable for submission to 'kilosort'
% (https://github.com/cortex-lab/KiloSort).
%INPUTS
% sessiondir: path to directory containing a collection of Neuralynx CSC
%   files with names in the form of 'csc*.ncs' where the '*' represents a
%   channel number.
% outfile: desired path to the file that is to be created.
%OUTPUTS
% Primary output is written to the file <outfile>.
% channelnums: the channel numbers extracted from each CSC file name.  Note
%   that there may be missing numbers in the set of files in <sessiondir>,
%   and also that any file whose length is less than the maximum length in
%   the session will be skipped over, and its number not included in
%   <channelnums>.  <channelnums> is thus a vector whose length is equal to
%   the number of rows in the array, <myDataMatrix>, written to the output
%   file.

%$Rev: 259 $
%$Date: 2017-08-02 14:37:57 -0400 (Wed, 02 Aug 2017) $
%$Author: dgibson $


files = dir(fullfile(sessiondir, 'csc*.ncs'));
filenames = {files.name};
toks = regexp(filenames, 'csc(\d+)\.ncs', 'tokens');
channelnums = cellfun(@(a) str2double(a{1}{1}), toks);
filesizes = [files.bytes];
isbad = filesizes < max(filesizes);
channelnums(isbad) = [];
channelnums = sort(channelnums);
for chidx = 1:length(channelnums)
    infilename = sprintf('csc%d.ncs', channelnums(chidx));
    [~, Samples] = dg_readCSC(fullfile(sessiondir, infilename));
    if chidx == 1
        myDataMatrix = zeros(length(channelnums), numel(Samples));
    end
    myDataMatrix(chidx, :) = reshape(Samples, 1, []);
end
fid = fopen(outfile, 'w');
if fid == -1
    error('nlx2kilosort:open', ...
        'Could not open %s for writing.', outfile);
end
fwrite(fid, myDataMatrix, 'int16');
fclose(fid);

