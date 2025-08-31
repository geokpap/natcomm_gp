function tempname = dg_mktempname(basename, ext)
%tempname = dg_mktempname(basename, ext)
% Creates a temporary filename starting with <basename> that does not
% conflict with any existing file.  <basename> can be an absolute or
% relative path.  <ext> is always appended verbatim at the end of
% <tempname> and should include a '.' if you want there to be one.

%$Rev: 43 $
%$Date: 2009-12-08 17:28:24 -0500 (Tue, 08 Dec 2009) $
%$Author: dgibson $

tempfilenum = 1;
tempname = sprintf('%s_temp%d%s', basename, tempfilenum, ext);
while exist(tempname)
    tempfilenum = tempfilenum + 1;
    tempname = sprintf('%s_temp%d%s', basename, tempfilenum, ext);
end
