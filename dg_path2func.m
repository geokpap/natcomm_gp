function funchandle = dg_path2func(pathstr)
%funchandle = dg_path2func(pathstr)
%INPUT
% pathstr: absolute pathname to a specific version of a function
%OUTPUT
% funchandle: function handle that will run that specific version

%$Rev: 132 $
%$Date: 2011-10-12 12:37:39 -0400 (Wed, 12 Oct 2011) $
%$Author: dgibson $

if isempty(pathstr)
    funchandle = [];
    return
end
if ~exist(pathstr, 'file')
    error('dg_path2func:nosuchfile', ...
        'The file %s does not exist.', pathstr);
end
[p,f,e] = fileparts(pathstr);
if ~isequal(e, '.m')
    warning('dg_path2func:ext', ...
        'Unfamiliar function file extension: "%s"', e);
end
oldwd = cd(p);
funchandle = str2func(f);
cd(oldwd);
