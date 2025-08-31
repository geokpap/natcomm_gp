function [mazespec, sessiondatestr] = lfp_findMazespec(date)
%mazespec = lfp_findMazespec
%mazespec = lfp_findMazespec(date)
% Returns a <mazespec> with additional field 'mazespec.calibstr', which
% contains a human-readable ID string for the calibration used.  <date> is
% optional; if given, it overrides the process for reading the session date
% from file timestamps, and it must be a single date value suitable for
% passing as the sole argument to Matlab's datenum function.  If no
% suitable calibration is found, then the default calibration is used and
% 'mazespec.calibstr' = 'default'.  <mazespec> contains the following
% fields:
%   cmperpixel
%   X0
%   X1
%   Y0
%   Y1
%   calibstr
%   sessiondatestr
% where the 'sessiondatestr' field contains a human-readable representation
% of the date used as the date of the session.

%$Rev: 69 $
%$Date: 2009-05-26 19:12:30 -0400 (Tue, 26 May 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 1
    date = [];
end
sessiondatestr = 'none';
% The file 'lfp_tmazecalib.mat' contains a variable named <tmazecalib>
% whose value is a vector of structs having fields 'calibname' (a
% human-readable ID string), 'datenum' (on which the calibration was
% recorded), 'x', 'y' (vectors of x and y coords of the numbered points in
% "T maze markers for rat from Hu Dan.doc", formerly known as "Rat T-maze
% markers from Hu Dan.doc"; note that only points [1:10 13 20 21] are used
% by dg_tmazeCalibAnal, and of those, [9 10 13 20 21] are only used for
% plotting).  The calib structs may be in any any order in the file.
load('lfp_tmazecalib.mat');
% sort <tmazecalib> chronologically:
tmazedates = cell(size(tmazecalib));
[tmazedates{:}] = deal(tmazecalib.datenum);
tmazedates = cell2mat(tmazedates);
[B,IX] = sort(tmazedates);
tmazecalib = tmazecalib(IX);

if isempty(date)
    files2try = { 'events.nev'
        'events.dat'
        'lfp1.ncs'
        'lfp1.dat'
        'vt1.nvt'
        'vt1.dat' };
    filedatenum = [];
    for k = 1:length(files2try)
        fn = fullfile(lfp_EventsDataDirs{end}, files2try{k});
        if exist(fn, 'file')
            fstruct = dir(fn);
            filedatenum = fstruct.datenum;
            sessiondatestr = datestr(filedatenum, 1);
            disp(sprintf('lfp_findMazespec found file %s dated %s', ...
                fstruct.name, sessiondatestr ));
            break
        end
    end
else
    filedatenum = datenum(date);
    sessiondatestr = datestr(filedatenum, 1);
    disp(sprintf('Using date %s', ...
        sessiondatestr ));
end
if isempty(filedatenum)
    warning('lfp_findMazespec:filedatenum', ...
        'Could not find suitable file to date session; using latest calib');
    calibnum = length(tmazecalib);
else
    calibdates = cell(size(tmazecalib));
    [calibdates{:}] = deal(tmazecalib.datenum);
    calibdates = cell2mat(calibdates);
    [Y, M, D] = datevec(calibdates);
    calibdates = ...
        datenum(Y, M, D, zeros(size(Y)), zeros(size(Y)), ones(size(Y)));
    oldenough = find(calibdates < filedatenum);
    if isempty(oldenough)
        first_c02_session_datenum = datenum('2005-11-27');
        if filedatenum >= first_c02_session_datenum
            warning('lfp_findMazespec:tooOld', ...
                'File is new enough to use first calib');
            calibnum = 1;
        else
            warning('lfp_findMazespec:tooOld', ...
                'There is no calibration file from before the session date');
            calibnum = 0;
        end
    else
        % Use the newest calib that's older than the session
        calibnum = oldenough(end);
    end
end
if calibnum > 0
    msg = sprintf('lfp_findMazespec: using T-maze calibration "%s"', ...
        tmazecalib(calibnum).calibname);
    mazespec = dg_tmazeCalibAnal(tmazecalib(calibnum));
    mazespec.calibstr = datestr(tmazecalib(calibnum).datenum, 1);
else
    msg = 'Using default calibration';
    mazespec.X0 = 100;
    mazespec.X1 = 160;
    mazespec.Y0 = 215;
    mazespec.Y1 = 245;
    myX2 = 460;
    myY3 = 370;
    myY4 = 100;
    cmperpixelX = 119.5/(myX2 - mazespec.X0);
    cmperpixelY = 73.8/(myY3-myY4);
    mazespec.cmperpixel = (cmperpixelX + cmperpixelY) / 2;
    mazespec.calibstr = 'default';
end
lfp_log(msg);
disp(msg);