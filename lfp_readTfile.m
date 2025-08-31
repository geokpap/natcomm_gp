function [header, spikes] = lfp_readTfile(filepath)
% LFP_READTFILE reads an MClust/BubbleClust T-file containing timestamps in
% units of 100 us (MClust's out-of-the-box loading engines use 10 us units).
% [header, spikes] = lfp_readTfile(filepath)
%   The <header> output does not include the delimiters '%%BEGINHEADER',
%   '%%ENDHEADER', nor their terminating line breaks.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

fid = fopen(filepath, 'r', 'b');
line = dg_ReadLn(fid);
if ~strcmp(line, '%%BEGINHEADER')
    error('lfp_readTfile:badfile', '%s is not a T-file', filepath);
end
header = [];
linebreak = sprintf('\n');
line = dg_ReadLn(fid);
while ~strcmp(line, '%%ENDHEADER')
    header = [ header line linebreak ];
    line = dg_ReadLn(fid);
end
[spikes, count] = fread(fid, inf, 'uint32');
fprintf(1, 'Read %d spikes from %s\n', count, filepath);
spikes = spikes * 100;
fclose(fid);
