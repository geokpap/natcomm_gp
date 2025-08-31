function dg_copyfile(src, dest)
%dg_copyfile(src, dest)
% This is a replacement for Matlab's 'copyfile', which can do bad things on
% a Mac.

%$Rev: 220 $
%$Date: 2015-05-28 12:27:28 -0400 (Thu, 28 May 2015) $
%$Author: dgibson $

if ispc
    [status, cmdout] = system(sprintf('copy "%s" "%s"', src, dest));
elseif ismac || isunix
    [status, cmdout] = system(sprintf('cp "%s" "%s"', src, dest));
else
    error('dg_copyfile:arch', ...
        'Unrecognized computer platform');
end
if status
    error('dg_copyfile:system', ...
        'Could not copy: %s', cmdout);
end

