function dg_mat2txt(pathname, A)
%DG_MAT2TXT saves a Matlab 2-D array as a text file.
%dg_mat2txt(pathname, A)
%  The columns within each row are comma-delimited with no whitespace.  Of
%  course, this means that you can NOT save character arrays that contain
%  commas!  (Does not work well with Notepad, use Wordpad or an Office
%  program.)

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

fid = fopen(pathname, 'w');
if fid < 0
    error('dg_mat2txt:nowrite', ...
        'Could not open file %s for writing.', pathname);
end
for row = 1:size(A,1)
    for col = 1:size(A,2)-1
        fprintf(fid, '%s,', mat2str(A(row,col)));
    end
    fprintf(fid, '%s\n', mat2str(A(row,end)));
end
fclose(fid);
