function dg_batchProcess(funchandle, varargin)
%BATCHPROCESS executes any function wrapped in error handling with logging,
%and then exits Matlab.
% batchProcess(funchandle, varargin)
%   <funchandle> is an arbitrary function handle, and <varargin> is its
%   argument list, with any dg_batchProcess options prepended.
%   dg_batchProcess options will be removed from the head of <varargin>
%   until a value which is not a dg_batchProcess option is encountered.
%OPTIONS
% 'leaveMarkerFileOnExit' - creates the file "dg_batchProcess_done.txt"
%   immediately before exiting.
% 'verbose' - eyup.

%$Rev: 310 $
%$Date: 2025-02-07 17:59:39 -0500 (Fri, 07 Feb 2025) $
%$Author: dgibson $

disp(datestr(now, 0));
fprintf('Matlab version %s\npathdef: %s\n', version, which('pathdef'));
fprintf('Current working directory: %s\n', pwd);
fprintf('funchandle: %s\n', func2str(funchandle));
fprintf('varargin: %s\n', dg_thing2str(varargin));
foundfirstarg = false;
markerflag = false;
verboseflag = false;
while ~foundfirstarg && ~isempty(varargin)
    if isequal(varargin{1}, 'leaveMarkerFileOnExit')
        markerflag = true;
        varargin(1) = [];
    elseif isequal(varargin{1}, 'verbose')
        verboseflag = true;
        varargin(1) = [];
    else
        foundfirstarg = true;
    end
end
workdir = pwd;
try
    if verboseflag
        disp('Evaluating <funchandle>...');
    end
    feval(funchandle, varargin{:});
    if verboseflag
        disp('... <funchandle> done.');
    end
catch e
    cd(workdir);
    if verboseflag
        disp('Caught error.');
    end
    if isempty(varargin)
        logmsg = sprintf('Error while processing %s', ...
            func2str(funchandle));
    else
        logmsg = sprintf('Error while processing %s(%s)', ...
            func2str(funchandle), dg_thing2str(varargin{1}) );
    end
    logmsg = sprintf('%s: %s\n%s', ...
        logmsg, e.identifier, e.message);
    for stackframe = 1:length(e.stack)
        logmsg = sprintf('%s\n%s\nline %d', ...
            logmsg, e.stack(stackframe).file, e.stack(stackframe).line);
    end
    disp(logmsg);
end
cd(workdir);
if markerflag
    if verboseflag
        disp('Writing marker file...');
    end
    filename = 'dg_batchProcess_done.txt';
    [fid, msg] = fopen('dg_batchProcess_done.txt', 'w');
    if fid == -1
        if verboseflag
            fprintf('Failed to open marker file %s.\n', filename);
        end
        warning('dg_batchProcess:marker', msg);
    else
        if verboseflag
            fprintf('Succesfully opened marker file %s.\n', filename);
        end
        fprintf(fid, '%s\n', datestr(now, 0));
        ST = fclose(fid);
        if verboseflag
            fprintf('''fclose'' returned %d.\n', ST);
        end
    end
    if verboseflag
        disp('Done attempting to write marker file.');
    end
end
if verboseflag
    disp('Exiting Matlab.');
end
exit;
