function line = dg_ReadLn(fid)
%DG_READLN reads one line of text from file ID <fid>.
% line = dg_ReadLn(fid)
%   reads one line of text from file ID <fid>, returns it as a string.
%   Essentially duplicates Matlab function fgetl.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

bufferIncrement = 131072; % 128K, enough to handle BBCLUST TFILE Version=1
line = char(zeros(1, bufferIncrement));
numchar = 0;
c = char(fread(fid, 1, 'uchar'));
linebreak = sprintf('\n');
while c ~= linebreak
    numchar = numchar + 1;
    if numchar > length(line)
        line = [ line char(zeros(1, bufferIncrement)) ];
    end
    line(numchar) = c;
    c = char(fread(fid, 1, 'uchar'));
end
line = line(1:numchar);