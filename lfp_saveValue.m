function lfp_saveValue(varname, varargin)
%lfp_saveValue(varname)
%lfp_saveValue(..., 'session')
%lfp_saveValue(..., 'overwrite')
% Companion function to lfp_loadValue.
% Saves the variable whose name is <varname> in a <varname>.mat file in the
% current data directory, lfp_DataDir.  If <varname>.mat already
% exists, then an error is raised.
%
%INPUT
% varname: string containing the name of the variable.
%OPTIONS
% 'session' - saves the .mat file in the root session directory if
%   lfp_DataDir is a subdirectory (i.e. a fragment directory).
% 'overwrite' - raises a warning instead of an error and overwrites the
%   existing file if there is one.

%$Rev: 374 $
%$Date: 2016-03-10 19:42:56 -0500 (Thu, 10 Mar 2016) $
%$Author: dgibson $

lfp_declareGlobals;
overwriteflag = false;
sessionflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'overwrite'
            overwriteflag = true;
        case 'session'
            sessionflag = true;
        otherwise
            error('lfp_saveValue:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if sessionflag
    pathname = fullfile(dg_findSessionRoot(lfp_DataDir), [varname '.mat']);
else
    pathname = fullfile(lfp_DataDir, [varname '.mat']);
end

if exist(pathname, 'file')
    if overwriteflag
        warning('lfp_saveValue:overwrite', ...
            'Overwriting file %s', pathname );
        evalin('caller', sprintf('save(''%s'', ''%s'');', pathname, varname));
    else
        error('lfp_saveValue:overwrite', ...
            'The file %s already exists', pathname );
    end
else
    evalin('caller', sprintf('save(''%s'', ''%s'');', pathname, varname));
end
if sessionflag
    fragtype = 'session';
else
    fragtype = 'fragment';
end
lfp_log(sprintf('Saved %s value: %s', fragtype, varname));

