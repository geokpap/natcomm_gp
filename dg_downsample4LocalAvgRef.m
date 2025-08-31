function dg_downsample4LocalAvgRef(filenames, sessiondir, N, varargin)
%dg_downsample4LocalAvgRef(filenames, sessiondir, N)
% dg_downsample all the files listed in <filenames>, after checking them
% as a group for timestamp conflicts.  If there are conflicts or disordered
% timestamps and the 'reconcile' option is not used, an error is raised.
%INPUTS
% filenames: names of files to downsample, either as a cell string vector
%   or a cell vector of cell arrays of strings (i.e. in the form of
%   <localavgrefs> used for supplying literal filenames to
%   dg_makeLocalAvgRefs).
% sessiondir: string
% N: integer
%OUTPUTS
% ...to same <sessiondir> that it reads inputs from, filenames constructed
%   as in 'suffix' option to dg_downsampleAndConvert.  Gracefully handles
%   filenames that already have a suffix indicating they were downsampled.
%OPTIONS
% 'outdir', outdir - writes output files to <outdir> instead of to
%   <sessiondir>.  If <outdir> does not already exist, it gets created.
% 'convert' - converts any .ncs files in <filenames> to .mat format.
% 'ncs' - for each filename in <filenames>, substitute the extension '.ncs'
%   for whatever extension (or lack thereof) each filename has before
%   attempting to read from the file.
% 'noTScheck' - skips the timestamp check; this should only be used in
%   cases where it is already known that there are no TS problems.
% 'overwrite' - raises a warning instead of an error and overwrites the
%   existing <dest> file if there is one.
% 'reconcile', recspec - reconciles timestamps against those of other CSC
%   files according to <recspec>.  If <recspec> is 'all', then every CSC
%   file in <filenames> is used as a source of frame timestamps to
%   reconcile; otherwise, <recspec> should be a cell vector of strings
%   giving the filenames (not including the directory) of files to be
%   reconciled.  Specifying 'reconcile' implicitly also entails
%   'noTScheck', since after reconciliation the timestamps are known to be
%   consistent.
% 'skip' - raises a warning instead of an error and skips the
%   existing <dest> file if there is one.
% 'verbose' - that.

argnum = 1;
convertflag = false;
downsample_opts = {};
ncsflag = false;
outdir = sessiondir;
recspec = '';
TScheckflag = true;
verboseflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'filter'
            downsample_opts{end+1} = 'filter'; %#ok<AGROW>
            argnum = argnum + 1;
            downsample_opts{end+1} = varargin{argnum}; %#ok<AGROW>
        case 'convert'
            convertflag = true;
        case 'ncs'
            ncsflag = true;
        case 'noTScheck'
            TScheckflag = false;
        case 'outdir'
            argnum = argnum + 1;
            outdir = varargin{argnum};
        case 'overwrite'
            downsample_opts{end+1} = varargin{argnum}; %#ok<AGROW>
        case 'reconcile'
            argnum = argnum + 1;
            recspec = varargin{argnum};
            TScheckflag = false;
        case 'skip'
            downsample_opts{end+1} = varargin{argnum}; %#ok<AGROW>
        case 'verbose'
            downsample_opts{end+1} = varargin{argnum}; %#ok<AGROW>
            verboseflag = true;
        otherwise
            error('dg_downsample4LocalAvgRef:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

if iscell(filenames{1})
    % This is a cell vector of cell arrays of strings, and we must flatten
    % it into a cell vector.
    filenamevector = {};
    for refidx = 1:length(filenames)
        for fileidx = 1:length(filenames{refidx})
            filenamevector{end+1} = filenames{refidx}{fileidx}; %#ok<AGROW>
        end
    end
    filenames = unique(filenamevector);
end

if ~isempty(recspec)
    if isequal(recspec, 'all')
        files2reconcile = filenames;
    else
        files2reconcile = recspec;
    end
    timestamps = {};
    % At this point, any files2reconcile that are not .mat files must be
    % Neuralynx files.
    for k = 1:length(files2reconcile)
        [~, ~, ext] = fileparts(files2reconcile{k});
        if verboseflag
            tic;
            fprintf('Reading timestamps to reconcile: %s\n', ...
                files2reconcile{k});
        end
        if strcmpi(ext, '.mat')
            load('-MAT', fullfile(sessiondir, files2reconcile{k}));
            timestamps{k} = dg_Nlx2Mat_Timestamps; %#ok<AGROW>
            clear dg_Nlx2Mat_Timestamps;
            clear dg_Nlx2Mat_Samples;
        else
            try
                timestamps{k} = dg_readCSC( ...
                    fullfile(sessiondir, files2reconcile{k})); %#ok<AGROW>
            catch e
                thefile = dir(fullfile(sessiondir, files2reconcile{k}));
                if thefile.bytes == 2^14
                    % Empty file, skip with warning
                    warning('dg_downsample4LocalAvgRef:empty3', ...
                        'Skipping empty file: %s', ...
                        fullfile(sessiondir, files2reconcile{k}));
                    timestamps{k} = []; %#ok<AGROW>
                    files2skip{end+1} = files2reconcile{k}; %#ok<AGROW>
                    continue
                else
                    errmsg = sprintf( ...
                        'dg_readCSC in dg_downsample4LocalAvgRef:\n%s\n%s', ...
                        dg_thing2str(fullfile(sessiondir, ...
                        files2reconcile{k})) );
                    errmsg = sprintf('%s\n%s\n%s', ...
                        errmsg, e.identifier, e.message);
                    for stackframe = 1:length(e.stack)
                        errmsg = sprintf('%s\n%s\nline %d', ...
                            errmsg, e.stack(stackframe).file, ...
                            e.stack(stackframe).line);
                    end
                    error('dg_downsample4LocalAvgRef:failed_dg_readCSC', ...
                        '%s', errmsg);
                end
            end
        end
        if verboseflag
            elapsed = toc;
            fprintf('Read timestamps in %.1f s\n', elapsed);
        end
        if any(timestamps{k} == 0)
            if skipbadflag
                warning('dg_downsample4LocalAvgRef:timestamps', ...
                    'File %s contains timestamps with value 0.', ...
                    files2reconcile{k});
                timestamps{k} = []; %#ok<AGROW>
                files2skip{end+1} = files2reconcile{k}; %#ok<AGROW>
                continue
            else
                error('dg_downsample4LocalAvgRef:timestamps', ...
                    'File %s contains timestamps with value 0.', ...
                    files2reconcile{k});
            end
        end
        if any(diff(timestamps{k}) <= 0)
            if skipbadflag
                warning('dg_downsample4LocalAvgRef:timestamps2', ...
                    'File %s nonincreasing timestamps.', ...
                    files2reconcile{k});
                timestamps{k} = []; %#ok<AGROW>
                files2skip{end+1} = files2reconcile{k}; %#ok<AGROW>
                continue
            else
                error('dg_downsample4LocalAvgRef:timestamps2', ...
                    'File %s nonincreasing timestamps.', ...
                    files2reconcile{k});
            end
        end
    end
    if verboseflag
            disp('Reconciling timestamps...');
    end
    canonicalstamps = dg_reconcileFrames(timestamps);
    downsample_opts{end+1} = 'reconcile';
    downsample_opts{end+1} = canonicalstamps;
end

[~, name, ext] = fileparts(filenames{1});
if TScheckflag
    if ncsflag || ~isequal(lower(ext), '.mat')
        if ncsflag
            infilename = [name '.ncs'];
        else
            infilename = filenames{1};
        end
        TS = dg_readCSC(fullfile(sessiondir, infilename));
        if any(diff(TS)<0)
            error('dg_downsample4LocalAvgRef:TS1', ...
                '%s contains corrupted timestamps.', infilename);
        end
        for fileidx = 2:length(filenames)
            fprintf('Loading file #%d\n', fileidx);
            if ncsflag
                [~, name] = fileparts(filenames{fileidx});
                infilename = [name '.ncs'];
            else
                infilename = filenames{fileidx};
            end
            TS2 = dg_readCSC(fullfile(sessiondir, infilename));
            if ~isequal(TS, TS2)
                error('dg_downsample4LocalAvgRef:TS2', ...
                    'Timestamps do not match across files.');
            end
        end
    else
        load(fullfile(sessiondir, filenames{1}), '-mat');
        TS = dg_Nlx2Mat_Timestamps;
        if any(diff(TS)<0)
            error('dg_downsample4LocalAvgRef:TS3', ...
                'TS, Elliot...');
        end
        for fileidx = 2:length(filenames)
            fprintf('Loading file #%d\n', fileidx);
            load(fullfile(sessiondir, filenames{fileidx}), '-mat');
            if ~isequal(TS, dg_Nlx2Mat_Timestamps) || any(diff(TS)<0)
                error('dg_downsample4LocalAvgRef:TS4', ...
                    'TS, Elliot...');
            end
        end
    end
end

for fileidx = 1:length(filenames)
    [~,name,ext] = fileparts(filenames{fileidx});
    Nprevtokens = regexp(name, '^(.*)_down(\d+).*$', 'tokens');
    if isempty(Nprevtokens)
        prevname = name;
        Nprev = 1;
    else
        prevname = Nprevtokens{1}{1};
        Nprev = str2double(Nprevtokens{1}{2});
    end
    [~, name, ext] = fileparts(filenames{fileidx});
    if ncsflag
        infilename = [name '.ncs'];
    else
        infilename = filenames{fileidx};
    end
    if convertflag
        ext2 = '.mat';
    else
        ext2 = ext;
    end
    outfilename = lower(sprintf('%s_down%d%s', prevname, N * Nprev, ext2));
    dg_downsample(fullfile(sessiondir, infilename), ...
        fullfile(outdir, outfilename), N, downsample_opts{:});
end


