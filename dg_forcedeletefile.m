function dg_forcedeletefile(filepath)
%dg_deletefile(filepath)
% Delete file without querying user even if it is read-only.
%NOTES
% See also dg_deletefile.

%$Rev: 186 $
%$Date: 2013-12-17 19:38:58 -0500 (Tue, 17 Dec 2013) $
%$Author: dgibson $

if ispc
    system(sprintf('del /F %s', filepath));
elseif ismac || isunix
    system(sprintf('rm -f %s', filepath));
else
    error('dg_copyfile:arch', ...
        'Unrecognized computer platform');
end
