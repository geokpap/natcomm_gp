function pid = dg_pid
%pid = dg_pid
% Returns the current process ID.
% If there is no executable mex file in the path, creates one in the user's
% MATLAB directory, which is assumed to be in the default location  as of
% Matlab version 7.10.0.499 (R2010a).
%NOTES
% A bit of a kludgerei, since in the Posix world 'getpid'  (see the .c
% files) is deprecated and who knows what goes on in the Windows world.
% But it works for now (23-Aug-2010).

%$Rev: 77 $
%$Date: 2010-08-24 18:45:35 -0400 (Tue, 24 Aug 2010) $
%$Author: dgibson $

if exist('dg_getpid') ~= 3
    if ismac
        destdir = fullfile(getenv('HOME'), 'Documents', 'MATLAB');
    elseif ispc
        % never did get this one to compile, though...
        homedir = getenv('HOME');
        if isempty(homedir)
            homedir = getenv('HOMEPATH');
        end
        destdir = fullfile(homedir, 'My Documents', 'MATLAB');
    elseif isunix
        destdir = fullfile(getenv('HOME'), 'matlab'); % as of R2009a
    else
        error('dg_pid:unkarch', ...
            'Hey, dude, like, what the hell IS that thing you''re using?');
    end
    [p, n] = fileparts(which('dg_pid.m'));
    if ispc
        str = computer;
        if isequal(str, 'PCWIN64')
            destext = '.mexw64';
        else
            destext = '.mexw32';
        end
        mex('-outdir', destdir, fullfile(p, 'dg_getpid_win.c'));
        movefile(fullfile(destdir, ['dg_getpid_win' destext]), ...
            fullfile(destdir, ['dg_getpid' destext]));
    else
        mex('-outdir', destdir, fullfile(p, 'dg_getpid.c'));
    end
    
    if isempty(which('dg_getpid'))
        error('dg_pid:nopidpath', ...
            'Please move the new dg_getpid executable from "%s" to a directory in your Matlab path.', ...
            destdir);
    end
end
pid = dg_getpid;
