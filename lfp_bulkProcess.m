function lfp_bulkProcess(funcH, argvals, varargin)
%BULKPROCESS performs offline processing on a list of
%argument values (e.g. directories).
% bulkProcess(funcH, argvals)
%   <funcH> is a function handle to a function of one argument that defines
%   a processing protocol.  <argvals> is a cell array containing the list
%   of argument values that are to be processed; if only one is to be
%   processed, it must still be a cell array, not a scalar or string.

% OPTIONS
%   'debug' - simply executes the function without try...catch error
%   handling.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

debugflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'debug'
            debugflag = true;
        otherwise
            error('bulkProcess:badoption', ...
                ['The option "' varargin{argnum} ...
                    '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

for val_idx = 1:length(argvals)
    myval = argvals{val_idx};
    if debugflag
        feval(funcH, myval);
    else
        try
            feval(funcH, myval);
        catch
            lfp_declareGlobals;
            if isempty(lfp_LogFileName)
                logname = 'lfp_lib.log';
                lfp_LogFileName = which(logname);
                if isempty(lfp_LogFileName)
                    lfp_LogFileName = fullfile(pwd, logname);
                end
            end
            [msgstr, msgid] = lasterr;
            logmsg = sprintf('Error while processing %s\n%s\n%s', ...
                myval, msgid, msgstr );
            lfp_log(logmsg);
            disp(logmsg);
        end
    end
end
