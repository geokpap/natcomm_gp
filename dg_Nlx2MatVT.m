function [ts x y a targ pts hdr] = dg_Nlx2MatVT(filename)
% Using <nrec> to denote the number of records:
% ts: row vector of uint64
% x: row vector of uint32
% y: row vector of uint32
% a: row vector of uint32
% targ: 50 x nrec array of uint32
% pts: 400 x nrec array of uint32
% hdr: 8k (8,192) row vector of uint16s that must be unpacked to convert to
%   text
%NOTES
%   The endian persuasion of the machine running this code is irrelevant,
% as the byte ordering in the Nlx file is determined by the machine that
% recorded the file, which is assumed to be Wintel (the only architecture
% supported by Neuralynx for Cheetah as of 27-May-2010).
%   This code ran 5x slower than Nlx2MatVT_v3 on a local test file
% (S36\acq10\VT1.Nvt), so it should only be used on non-Windows machines.

%$Rev: 55 $
%$Date: 2010-06-01 19:06:52 -0400 (Tue, 01 Jun 2010) $
%$Author: dgibson $

fid = fopen(filename);
if fid == -1
    error('dg_Nlx2MatVT:badfn', ...
        'Could not open file "%s"', filename);
end
A = fread(fid, Inf, '*uint16');
fclose(fid);
hdrsz = 8 * 2^10;   % in uint16s, = 16kB
recsz = 7 + 400*2 + 7 + 50*2; % in uint16s
xrow = 7 + 400*2 + 2; % position of extracted x coordinate in record
hdr = A(1:hdrsz);
numrecs = fix((length(A) - hdrsz) / recsz);
if numrecs * recsz ~= length(A) - hdrsz
    warning('dg_Nlx2MatVT:missdata', ...
        'File is missing data');
end
A = reshape(A(hdrsz+1:end), recsz, numrecs); % one rec per column
if ~all(A(1,:) == 2048)
    warning('dg_Nlx2MatVT:baddata', ...
        'File is corrupted');
end
ts = uint64(zeros(1, size(A,2)));
for k = 0:3
    ts = bitor(ts, bitshift(uint64(A(k+4,:)), k*16));
end
x = bitor(uint32(A(xrow,:)), bitshift(uint32(A(xrow+1,:)), 16));
y = bitor(uint32(A(xrow+2,:)), bitshift(uint32(A(xrow+3,:)), 16));
a = bitor(uint32(A(xrow+4,:)), bitshift(uint32(A(xrow+5,:)), 16));
targ = bitor(uint32(A(xrow+(6:2:105),:)), ...
    bitshift(uint32(A(xrow+(7:2:105),:)), 16));
pts = bitor(uint32(A(7+(1:2:800),:)), ...
    bitshift(uint32(A(7+(2:2:800),:)), 16));

% Just for the record, I attempted the following de-parallelization of the
% code in hopes that skipping computation of zero values would make it
% faster, but the net result was to slow it down instead:
% targ = zeros(50, numrecs);
% for rec = 1:numrecs
%     for targnum = 1:50
%         targ(targnum, rec) = bitor(uint32(A(xrow+5+2*targnum, rec)), ...
%             bitshift(uint32(A(xrow+6+2*targnum, rec)), 16));
%         if targ(targnum, rec) == 0
%             break
%         end
%     end
% end
% pts = zeros(400, numrecs);
% for rec = 1:numrecs
%     for ptnum = 1:400
%         pts(ptnum, rec) = bitor(uint32(A(6+2*ptnum, rec)), ...
%             bitshift(uint32(A(7+2*ptnum, rec)), 16));
%         if pts(ptnum, rec) == 0
%             break
%         end
%     end
% end


