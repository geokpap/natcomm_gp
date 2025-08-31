function dirs = dg_getDirs(dirpath)
%DG_GETDIRS returns a list of directories
%dirs = dg_getDirs(dirpath)
%   Returns a cell vector of names of directories, not including '.' and
%   '..', that reside in directory <dirpath>.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

files = dir(dirpath);
dirs = {files(find(cat(1,files.isdir))).name};
dirs(find(ismember(dirs, '.'))) = [];
dirs(find(ismember(dirs, '..'))) = [];
