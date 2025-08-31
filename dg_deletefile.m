function dg_deletefile(filepath)
%dg_deletefile(filepath)
% This is a replacement for Matlab's 'delete', which might be perfectly
% fine, but I got superstitious after being burned by 'copyfile'.
%NOTES
% See also dg_forcedeletefile.

%$Rev: 186 $
%$Date: 2013-12-17 19:38:58 -0500 (Tue, 17 Dec 2013) $
%$Author: dgibson $

if ispc
    system(sprintf('del %s', filepath));
elseif ismac || isunix
    system(sprintf('rm %s', filepath));
else
    error('dg_copyfile:arch', ...
        'Unrecognized computer platform');
end
