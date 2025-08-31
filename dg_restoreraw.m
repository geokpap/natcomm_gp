function dg_restoreraw(dirpath)
% Finds every file in <dirpath> whose name starts with 'raw' and renames it
% to the same name without the 'raw', replacing any existing file that
% has the new name.  Does NOT open any subdirectories in <dirpath>, i.e.
% operates only at files that reside immediately at the top level of
% <dirpath>.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

files = dir(dirpath);
for fileidx = 1:length(files)
    if ~files(fileidx).isdir
        if regexp(files(fileidx).name, '^raw')
            newname = files(fileidx).name(4:end);
            movefile(fullfile(dirpath, files(fileidx).name), ...
                fullfile(dirpath, newname), 'f');
        end
    end
end

            
