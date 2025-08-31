function dg_normalizeTS(sessiondir, outdir, varargin)
% Finds the earliest timestamp in any CSC, SE, or Events file in the
% session and subtracts it from all timestamps in all those files,
% effectively setting the start time of the session to zero.
%INPUTS
% sessiondir: directory containing Neuralynx session files.
% outdir: directory in which to create new 'fixed' session folder to save
%   normalized data.  If not specified or empty, then <outdir> is the
%   parent of <sessiondir>.  The fixed sessiondir has the same name as
%   <sessiondir>, but with 'fixed' appended (note that if <sessiondir> ends
%   with a directory delimiter, it will contain a new directory named
%   'fixed').
%OUTPUTS
% All outputs go to files in a folder with same session ID but the string
% 'fixed' appended.
%OPTIONS
% 'test' - only performs normalization if the last timestamp in ALL files
%   is greater than one week (Plexon OFS can handle one week, but not two
%   weeks).  Also, skips normalization if <outdir> already exists.
% 'express' - don't bother with the whole timestamp search for the first
%   timestamp of the session; just read the first CSC file's first TS and
%   set time zero one minute earlier.  Raises an error if any subsequent
%   timestamps are before that.

%$Rev: 288 $
%$Date: 2022-02-04 17:39:43 -0500 (Fri, 04 Feb 2022) $
%$Author: dgibson $

[parent, name] = fileparts(sessiondir);
if nargin < 2 || isempty(outdir)
    outdir = fullfile(parent, [name 'fixed']);
else
    outdir = fullfile(outdir, [name 'fixed']);
end

expressflag =false;
testflag = false;

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'express'
            expressflag = true;
        case 'test'
            testflag = true;
        otherwise
            error('dg_normalizeTS:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

suffix = dg_findRealSessionFiles(sessiondir);
files = dir(sessiondir);
% find earliest and latest timestamps
startTS = Inf;
endTS = 0;
numfilesread = 0;
if testflag
    if exist(outdir, 'dir')
        warning('dg_normalizeTS:exists', ...
            '%s already exists, skipping.', outdir);
        return
    end
    if exist(fullfile(sessiondir, 'fixed'), 'dir')
        warning('dg_normalizeTS:fixed', ...
            '%s already contains a ''fixed'' subdir, skipping.', sessiondir);
        return
    end
end
for k = 1:length(files)
    % If a file has gotten stepped on, it might terminate in the middle of
    % the header, so we use "<=":
    if files(k).isdir || files(k).bytes <= 16384
        continue
    end
    [~, ~, ext] = fileparts(files(k).name);
    ext = lower(ext);
    if ~ismember(ext, {'.nev' '.ncs' '.nse'})
        continue
    end
    if expressflag
        if numfilesread
            break
        end
        if isempty(regexpi( files(k).name, ...
                sprintf('^CSC\\d+%s\\.ncs$', suffix), 'once' ))
            continue
        end
    end
    fprintf( 'Reading timestamps from %s.\n', ...
        fullfile(sessiondir, files(k).name) );
    try
        switch lower(ext)
            case '.nev'
                TS = dg_readEvents(fullfile(sessiondir, files(k).name));
            case '.ncs'
                if expressflag
                    TS = dg_readCSC(fullfile(sessiondir, files(k).name), ...
                        'mode', 3, 1);
                else
                    TS = dg_readCSC(fullfile(sessiondir, files(k).name));
                end
            case '.nse'
                TS = dg_readSpike(fullfile(sessiondir, files(k).name));
        end
    catch
        warning('dg_normalizeTS:badfile', ...
            'Unreadable file: %s', fullfile(sessiondir, files(k).name));
        continue
    end
    numfilesread = numfilesread + 1;
    if testflag
        if TS(end) < 7*24*3600e6
            fprintf( '%s does not need normalizing, skipping.\n', ...
                sessiondir );
            return
        end
    end
    if TS(1) < startTS
        startTS = TS(1);
    end
    if TS(end) > endTS
        endTS = TS(end);
    end
end
if numfilesread
    fprintf('dg_normalizeTS found earliest timestamp.\n');
else
    fprintf('No Nlx files found.\n');
    return
end
if ~exist(outdir, 'dir')
    dg_mkdir_p(outdir);
end
if expressflag
    startTS = startTS - 60e6;
end
% set time scale to start at 0
for k = 1:length(files)
    % skip directories and empty files
    if files(k).isdir || files(k).bytes <= 16384
        continue
    end
    [~, name, ext] = fileparts(files(k).name);
    switch lower(ext)
        case '.nev'
            % skip files that aren't "real session files":
            if isempty(regexpi( name, sprintf('^Events%s$', suffix), 'once' ))
                continue
            end
            fprintf('%s Normalizing TS for %s\n', ...
                datestr(now), files(k).name);
            [TS, TTL, ES, Hdr, nlxjunk] = dg_readEvents( ...
                fullfile(sessiondir, files(k).name) );
            TS = TS - startTS;
            if TS(1) < 0
                error('dg_normalizeTS:badstartTS', ...
                    '<startTS> is not the earliest timestamp.');
            end
            dg_writeNlxEvents( fullfile(outdir, files(k).name), ...
                TS, TTL, ES, Hdr, 'nlxjunk', nlxjunk );
        case '.ncs'
            if isempty(regexpi( name, sprintf('^CSC\\d+%s$', suffix), 'once' ))
                continue
            end
            fprintf('%s Normalizing TS for %s\n', ...
                datestr(now), files(k).name);
            [TS, Samples, Hdr, nlxjunk] = dg_readCSC( ...
                fullfile(sessiondir, files(k).name) );
            TS = TS - startTS;
            if TS(1) < 0
                error('dg_normalizeTS:badstartTS', ...
                    '<startTS> is not the earliest timestamp.');
            end
            dg_writeCSC( fullfile(outdir, files(k).name), ...
                TS, Samples, Hdr, 'nlxjunk', nlxjunk );
        case '.nse'
            if isempty(regexpi( name, sprintf('^SE\\d+%s$', suffix), 'once' ))
                continue
            end
            fprintf('%s Normalizing TS for %s\n', ...
                datestr(now), files(k).name);
            [TS, Samples, Hdr, nlxjunk] = dg_readSpike( ...
                fullfile(sessiondir, files(k).name) );
            TS = TS - startTS;
            if TS(1) < 0
                error('dg_normalizeTS:badstartTS', ...
                    '<startTS> is not the earliest timestamp.');
            end
            dg_writeSpike( fullfile(outdir, files(k).name), ...
                TS, Samples, Hdr, 'nlxjunk', nlxjunk );
    end
end

