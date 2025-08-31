function status = lfp_loadValue(varname, varargin)
%status = lfp_loadValue(varname)
% Companion function to lfp_saveValue.
% In the workspace from which lfp_loadValue was called, loads the variable
% named <varname> from a .mat file of the same name in the current data
% directory (lfp_DataDir) and returns status 1.  If no such file exists,
% then it is loaded from the root session directory and returns status 2.
% If no such file exists there either, nothing happens, and status is -1.
%OPTIONS
% 'ic' - ignore case when searching for .mat file even if running under an
%   OS which is case sensitive.  If more than one file matches, the first
%   one in the array returned by the 'dir' function is used.

%$Rev: 311 $
%$Date: 2013-10-20 22:07:36 -0400 (Sun, 20 Oct 2013) $
%$Author: dgibson $

icflag = false;

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'ic'
            icflag = true;
        otherwise
            error('lfp_loadValue:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end
lfp_declareGlobals;
status = 1;
filename2try = [varname '.mat'];
filename = dg_fileExists(lfp_DataDir, filename2try, icflag);
if isempty(filename)
    filedir = dg_findSessionRoot(lfp_DataDir);
    filename = dg_fileExists(filedir, filename2try, icflag);
    if isempty(filename)
        status = -1;
    else
        if icflag
            filename = filename{1};
        end
        evalin('caller', sprintf( 'load(''%s'');', ...
            fullfile(filedir, filename) ));
        status = 2;
        lfp_log(sprintf('Loaded session value: %s', varname));
    end
else
    if icflag
        filename = filename{1};
    end
    evalin('caller', sprintf('load(''%s'', ''-mat'');', ...
        fullfile(lfp_DataDir, filename)));
    lfp_log(sprintf('Loaded fragment value: %s', varname));
end

