function dg_downsampleAndConvert(sessiondir, destdir, N, varargin)
% Downsample by factor of N any CSC files in <sessiondir> and save as .mat
% files in <destdir>; convert Nlx events file(s) to .MAT files; ignore
% .nse, .ntt, and .nvt files and any .dat files whose names indicate those
% formats; copy any other files without modification.  Do not process
% subdirectories, collect $200, or do anything else you might want to do.
% Passes optional args on to dg_downsample.  If N = 1, then instead of
% calling dg_downsample, dg_Nlx2Mat is called.  If <destdir> does not
% exist, it is created.
%   CSC files are identified on the basis of having the extension '.ncs' or
% having names that begin with either 'lfp' or 'csc' followed by a string
% of digits.
%OPTIONS
% These options affect dg_downsampleAndConvert only, and are not passed on
%   to dg_downsample:
%   'allfiles' - strictly requires <canonicalstamps> to contain timestamps
%       that are common to literally ALL of the files that are chosen for
%       downsampling, even when large numbers of frames have to be skipped.
%       See <files2downsample> in code.  Implementation is done in
%       dg_reconcileFrames.
%   'allowedmissing', nums - skips any source files whose electrode number
%       is a member <nums>.  It is assumed that the first consecutive
%       string of digits in the file name is the electrode number.
%   'filenames', filenames - <filenames> is a cell array of file names (not
%       including directory) to downsample (and potentially reconcile
%       with each other).  This overrides the default method of identifying
%       CSC files.
%   'nlx' - writes output to a Neuralynx-format CSC file with the same file
%       extension as the input file. Works only on Windows.
%   'reconcile', recspec - reconciles timestamps against those of other CSC
%       files according to <recspec>.  If <recspec> is 'all', then every
%       CSC file found is used as a source of frame timestamps to
%       reconcile; otherwise, <recspec> should be a cell vector of strings
%       giving the filenames (not including the directory) of files to be
%       reconciled.  Any file that does not contain all of the timestamps
%       in the reconciled list is excluded from the downsampling.  Frames
%       with timestamps that are not in the reconciled list are excluded
%       before downsampling.
%   'skipbad' - instead of raising an error when defective files are found,
%       skip them and exclude them from further processing.
%   'skipmisc' - does not copy any files other than event files and CSC
%       files.
%   'suffix' - appends the string "_down<N>" to the output filenames
%       before the extension, where <N> is the downsampling factor.
% All other options are passed on to dg_downsample.  The following have
%   effects in dg_downsampleAndConvert as well as in dg_downsample:
%   'overwrite' - also applies to files that are simply being copied
%   'skip' - also applies to files that are simply being copied
%   'verbose' - no comment
%NOTES
% Unfortunately, there is not much space saving that results from
% downsampling *and* converting, since the .mat files contain double floats
% whereas the Nlx files contain integers.  E.g. downsampling 5X and
% converting to .mat only saves about 25% on space.
%BULK PROCESSING EXAMPLE:
% function bulkdownsample(sessiondir)
% [path, sessionID] = fileparts(sessiondir);
% [path, animalID] = fileparts(path);
% destdir = fullfile('C:\downsampled', animalID, sessionID);
% dg_downsampleAndConvert(sessiondir, destdir, 5, ...
%     'verbose', 'overwrite');


%$Rev: 272 $
%$Date: 2021-02-12 12:46:53 -0500 (Fri, 12 Feb 2021) $
%$Author: dgibson $

if ~ischar(sessiondir)
    error('dg_downsampleAndConvert:sessiondir', ...
        '%s', '<sessiondir> must be a string.');
end
if ~ischar(destdir)
    error('dg_downsampleAndConvert:destdir', ...
        '%s', '<destdir> must be a string.');
end
if ~isnumeric(N) || numel(N) > 1
    error('dg_downsampleAndConvert:N', ...
        '%s', '<N> must be a numeric scalar.');
end

allfilesflag = false;
args2delete = [];
filenames = {};
nlxflag = false;
overwriteflag = false;
recspec = '';
skipflag = false;
skipbadflag = false;
skipmiscflag = false;
suffixflag = false;
trodes2skip = [];
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'allfiles'
            allfilesflag = true;
        case 'allowedmissing'
            args2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            trodes2skip = varargin{argnum};
            args2delete(end+1) = argnum; %#ok<*AGROW>
        case 'filenames'
            args2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            filenames = varargin{argnum};
            args2delete(end+1) = argnum; %#ok<*AGROW>
        case 'nlx'
            nlxflag = true;
            args2delete(end+1) = argnum; %#ok<*AGROW>
        case 'overwrite'
            overwriteflag = true;
        case 'reconcile'
            args2delete(end+1) = argnum;
            argnum = argnum + 1;
            recspec = varargin{argnum};
            args2delete(end+1) = argnum;
        case 'skip'
            skipflag = true;
        case 'skipbad'
            skipbadflag = true;
            args2delete(end+1) = argnum; %#ok<*AGROW>
        case 'skipmisc'
            args2delete(end+1) = argnum;
            skipmiscflag = true;
        case 'suffix'
            args2delete(end+1) = argnum;
            suffixflag = true;
        case 'verbose'
            verboseflag = true;
    end
    argnum = argnum + 1;
end
varargin(args2delete) = [];

filenames2skip = {'vt1.dat' 'vt1.nvt'};
files2skip = {};
regexps2skip = {'^sc\d+\.' '^tt\d+\.'};
files = dir(sessiondir);
files2downsample = cell(0,2);
for k = 1:length(files)
    if ~files(k).isdir && ~ismember(files(k).name, filenames2skip)
        skipthis = false;
        for j = 1:length(regexps2skip)
            if ~isempty(regexpi(files(k).name, regexps2skip{j}))
                skipthis = true;
                break
            end
        end
        if ~isempty(filenames) && ~ismember(files(k).name, filenames)
            skipthis = true;
        end
        if skipthis
            continue
        end
        [p, name, ext] = fileparts(files(k).name); %#ok<*ASGLU>
        if ismember(ext, {'.nse', '.ntt', '.nvt'})
            continue
        end
        switch exist(destdir) %#ok<EXIST>
            case 7  % directory already exists, good.
            case 2  % there is already a file with the same name
                error('dg_downsampleAndConvert:baddestdir', ...
                    'There is an ordinary file with the same name as destdir: %s', ...
                    destdir);
            case 0  % directory does not exist
                mkdir(destdir);
        end
        if isempty(filenames) && strcmpi(ext, '.ncs') || ...
                isempty(filenames) && ...
                ~isempty(regexpi(name, '^(lfp|csc)\d+')) || ...
                ~isempty(filenames) && ismember(files(k).name, filenames)
            % It's a CSC file
            if nlxflag
                destfilename = [name ext];
            else
                destfilename = [name '.mat'];
            end
            files2downsample(end+1, :) = {files(k).name destfilename};
        else
            if ~isempty(regexpi(name, '^events$')) || ...
                    isequal(lower(ext), '.nev')
                % It's an events file
                srcfile = fullfile(sessiondir, files(k).name);
                if nlxflag
                    destfile = fullfile(destdir, [name ext]);
                    if exist(destfile,'file')
                        if overwriteflag
                            dg_copyfile(srcfile, destfile);
                        elseif ~skipflag
                            error('dg_downsampleAndConvert:fileexists1', ...
                                'The destination file %s already exists', ...
                                matfilename );
                        end
                    else
                        dg_copyfile(srcfile, destfile);
                    end
                else
                    matfilename = fullfile(destdir, [name '.mat']);
                    try
                        if exist(matfilename,'file')
                            if overwriteflag
                                dg_Nlx2Mat(srcfile, 'dest', destdir);
                            elseif ~skipflag
                                error('dg_downsampleAndConvert:fileexists2', ...
                                    'The destination file %s already exists', ...
                                    matfilename );
                            end
                        else
                            dg_Nlx2Mat(srcfile, 'dest', destdir);
                        end
                    catch e
                        thefile = dir(srcfile);
                        if thefile.bytes == 2^14
                            % Empty file, skip with warning
                            warning('dg_downsampleAndConvert:empty', ...
                                'Skipping empty file: %s', srcfile);
                            files2skip{end+1} = files(k).name;
                            continue
                        else
                            errmsg = sprintf( ...
                                'dg_Nlx2Mat in dg_downsampleAndConvert:\n%s\n%s', ...
                                dg_thing2str(srcfile), ...
                                dg_thing2str(destdir) );
                            errmsg = sprintf('%s\n%s\n%s', ...
                                errmsg, e.identifier, e.message);
                            for stackframe = 1:length(e.stack)
                                errmsg = sprintf('%s\n%s\nline %d', ...
                                    errmsg, e.stack(stackframe).file, ...
                                    e.stack(stackframe).line);
                            end
                            error('dg_downsampleAndConvert:failedNlx2Mat', ...
                                '%s', errmsg);
                        end
                    end
                end
            else
                % It's a miscellaneous file
                destfilename = fullfile(destdir, files(k).name);
                if exist(destfilename,'file')
                    if overwriteflag
                        dg_copyfile(fullfile(sessiondir, files(k).name), ...
                            destfilename);
                    elseif ~skipflag
                        error('dg_downsampleAndConvert:fileexists3', ...
                            'The destination file %s already exists', ...
                            destfilename );
                    end
                elseif ~skipmiscflag
                    dg_copyfile(fullfile(sessiondir, files(k).name), ...
                        destfilename);
                end
            end
        end
    end
end

if ~isempty(recspec)
    if isequal(recspec, 'all')
        files2reconcile = files2downsample(:,1);
    else
        files2reconcile = recspec;
    end
    % Each element of <timestamps> contains one file's timestamps vector,
    % and is empty if either the file is empty or the file has been flagged
    % as bad pursuant to 'skipbadflag'.
    timestamps = {};
    % At this point, any files2reconcile that are not .mat files must be
    % Neuralynx files.
    for k = 1:length(files2reconcile)
        [p, name, ext] = fileparts(files2reconcile{k});
        if verboseflag
            tic;
            fprintf('Reading timestamps to reconcile: %s\n', ...
                files2reconcile{k});
        end
        if strcmpi(ext, '.mat')
            load('-MAT', fullfile(sessiondir, files2reconcile{k}));
            timestamps{k} = dg_Nlx2Mat_Timestamps;
            clear dg_Nlx2Mat_Timestamps;
            clear dg_Nlx2Mat_Samples;
        else
            try
                timestamps{k} = dg_readCSC( ...
                    fullfile(sessiondir, files2reconcile{k}));
            catch e
                thefile = dir(fullfile(sessiondir, files2reconcile{k}));
                if thefile.bytes == 2^14
                    % Empty file, skip with warning
                    warning('dg_downsampleAndConvert:empty3', ...
                        'Skipping empty file: %s', ...
                        fullfile(sessiondir, files2reconcile{k}));
                    timestamps{k} = [];
                    files2skip{end+1} = files2reconcile{k};
                    continue
                else
                    errmsg = sprintf( ...
                        'dg_readCSC in dg_downsampleAndConvert:\n%s\n%s', ...
                        dg_thing2str(fullfile(sessiondir, ...
                        files2reconcile{k})) );
                    errmsg = sprintf('%s\n%s\n%s', ...
                        errmsg, e.identifier, e.message);
                    for stackframe = 1:length(e.stack)
                        errmsg = sprintf('%s\n%s\nline %d', ...
                            errmsg, e.stack(stackframe).file, ...
                            e.stack(stackframe).line);
                    end
                    error('dg_downsampleAndConvert:failed_dg_readCSC', ...
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
                warning('dg_downsampleAndConvert:timestamps', ...
                    'File %s contains timestamps with value 0.', ...
                    files2reconcile{k});
                timestamps{k} = [];
                files2skip{end+1} = files2reconcile{k};
                continue
            else
                error('dg_downsampleAndConvert:timestamps', ...
                    'File %s contains timestamps with value 0.', ...
                    files2reconcile{k});
            end
        end
        if any(diff(timestamps{k}) <= 0)
            if skipbadflag
                warning('dg_downsampleAndConvert:timestamps2', ...
                    'File %s nonincreasing timestamps.', ...
                    files2reconcile{k});
                timestamps{k} = [];
                files2skip{end+1} = files2reconcile{k};
                continue
            else
                error('dg_downsampleAndConvert:timestamps2', ...
                    'File %s nonincreasing timestamps.', ...
                    files2reconcile{k});
            end
        end
    end
    if verboseflag
            disp('Reconciling timestamps...');
    end
    if allfilesflag
        dg_reconcileFramesOpts = {'allfiles'};
    else
        dg_reconcileFramesOpts = {};
    end
    [canonicalstamps, activefnums, includedfile] = ...
        dg_reconcileFrames(timestamps, dg_reconcileFramesOpts{:});
    if any(~includedfile(activefnums))
        excludedfilestr = '';
        for badfnum = reshape( ...
                activefnums(~includedfile(activefnums)), 1, [] )
            excludedfilestr = sprintf( '%s\n%s', excludedfilestr, ...
                files2reconcile{badfnum} );
        end
        warning('dg_downsampleAndConvert:timestamps3', ...
            'The following files have been excluded due to mismatching timestamps:%s', ...
            excludedfilestr);
    end
    varargin{end+1} = 'reconcile';
    varargin{end+1} = canonicalstamps;
end

% At this point, <activefnums> is an index into <timestamps> that points to
% non-empty cells.  <timestamps{k}> is returned via dg_readCSC from
% <files2reconcile{k}>, so <activefnums> is also an index into
% <files2reconcile>.  <files2reconcile> is identical to <files2downsample>
% unless overridden by <recspec> (i.e. by the literal enumeration form of
% the 'reconcile' option), so we have two cases: (1) 1

if N == 1
    for k = 1:size(files2downsample,1)
        if verboseflag
            fprintf('Converting %s\n', files2downsample{k,1});
        end
        [p,n] = fileparts(files2downsample{k,1});
        if exist(fullfile(destdir, [n '.mat']), 'file')
            if ~skipflag
                error('dg_downsampleAndConvert:fileexists4', ...
                    'The destination file %s already exists', ...
                    fullfile(destdir, files2downsample{k,1}));
            end
            continue
        else
            try
                dg_Nlx2Mat(fullfile(sessiondir, files2downsample{k,1}), ...
                    'dest', destdir);
            catch e
                thefile = dir(fullfile(sessiondir, files2downsample{k,1}));
                if thefile.bytes == 2^14
                    % Empty file, skip with warning
                    warning('dg_downsampleAndConvert:empty2', ...
                        'Skipping empty file: %s', ...
                        fullfile(sessiondir, files2downsample{k,1}));
                    continue
                else
                    errmsg = sprintf( ...
                        'dg_Nlx2Mat arguments in dg_downsampleAndConvert:\n%s\n%s', ...
                        dg_thing2str(fullfile(sessiondir, files2downsample{k,1})), ...
                        dg_thing2str(destdir) );
                    errmsg = sprintf('%s\n%s\n%s', ...
                        errmsg, e.identifier, e.message);
                    for stackframe = 1:length(e.stack)
                        errmsg = sprintf('%s\n%s\nline %d', ...
                            errmsg, e.stack(stackframe).file, ...
                            e.stack(stackframe).line);
                    end
                    error('dg_downsampleAndConvert:failedNlx2Mat', ...
                        '%s', errmsg);
                end
            end
        end
    end
else
    for k = 1:size(files2downsample,1)
        if ismember(files2downsample{k,1}, files2skip)
            continue
        end
        toks = regexpi(files2downsample{k,1}, '^\D*(\d+)', 'tokens');
        trodenum = str2double(toks{1}{1});
        if ismember(trodenum, trodes2skip)
            if verboseflag
                fprintf('Skipping %s\n', files2downsample{k,1});
            end
            continue
        end
        if verboseflag
            fprintf('Downsampling %s\n', files2downsample{k,1});
        end
        outfilename = files2downsample{k,2};
        if suffixflag
            [p,name,ext] = fileparts(outfilename);
            outfilename = sprintf('%s_down%d%s', name, N, ext);
        end
        try
            dg_downsample(fullfile(sessiondir, files2downsample{k,1}), ...
                fullfile(destdir, outfilename), N, varargin{:});
        catch e
            thefile = dir(fullfile(sessiondir, files2downsample{k,1}));
            if thefile.bytes == 2^14
                % Empty file, skip with warning
                warning('dg_downsampleAndConvert:empty4', ...
                    'Skipping empty file: %s', ...
                    fullfile(sessiondir, files2downsample{k,1}));
                continue
            else
                errmsg = sprintf( ...
                    'dg_downsample in dg_downsampleAndConvert:\n%s\n%s', ...
                    dg_thing2str(fullfile(sessiondir, files2downsample{k,1})), ...
                    dg_thing2str(fullfile(destdir, outfilename)) );
                errmsg = sprintf('%s\n%s\n%s', ...
                    errmsg, e.identifier, e.message);
                for stackframe = 1:length(e.stack)
                    errmsg = sprintf('%s\n%s\nline %d', ...
                        errmsg, e.stack(stackframe).file, ...
                        e.stack(stackframe).line);
                end
                error('dg_downsampleAndConvert:failed_dg_downsample', ...
                    '%s', errmsg);
            end
        end
    end
end
end

