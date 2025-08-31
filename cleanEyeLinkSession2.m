function cleanEyeLinkSession2(sessiondir)
% Same as 'cleanEyeLinkSession' but for csc126 only.
% 9-Sep-2024: added error checks in "Make sure all required files are
% present" block.

tempdir = 'tmp';
mkdir('tmp');
outroot = '/annex2/analysis/dgibson/eyelink2';
prezIDs = {'prez' 'Prez'};
debIDs = {'debbie' 'LittleDebbie'};

% Convert animalIDs to standard ones {Prez|Debbie}.
pathparts = strsplit(sessiondir, '/');
if any(ismember(pathparts, prezIDs))
    pathparts{ismember(pathparts, prezIDs)} = 'Prez';
elseif any(ismember(pathparts, debIDs))
    pathparts{ismember(pathparts, debIDs)} = 'Debbie';
end
outdir = fullfile(outroot, pathparts{4:end});
fprintf('outdir = %s\n', outdir);
if ~exist(outdir, 'dir')
    dg_mkdir_p(outdir);
end

dg_Nlx2Mat_Timestamps = [];
dg_Nlx2Mat_Samples = [];
dg_Nlx2Mat_SamplesUnits = [];
% Events file must be first in <filenames> for file version number
% recognition to work (as in <altevtfilenames>).
filenames = { 'events.nev' 'csc126.ncs'  };
altevtfilenames = {'events_0001.nev'
    'events_0002.nev'
    'events_0003.nev'
    'events_0004.nev'};
dirfiles = dir(sessiondir);
dirfilenames = {dirfiles(~[dirfiles.isdir]).name};

% Make sure all required files are present.
suffix = '';
for fidx = 1:length(filenames)
    % <infiles> must be taken verbatim from <dirfilenames>, not from
    % anything that might only match case-insensitively.
    if ismember(filenames{fidx}, lower(dirfilenames))
        infiles{fidx} = dirfilenames{ ...
            ismember(lower(dirfilenames), filenames{fidx}) }; %#ok<*AGROW> 
    elseif fidx == 1 && ...
            any(ismember( altevtfilenames, lower(dirfilenames) ))
        vernum = find(ismember(altevtfilenames, lower(dirfilenames)));
        suffix = sprintf('_%04d', vernum);
        infiles{fidx} = dirfilenames{ismember(lower(dirfilenames), ...
            altevtfilenames)};
    elseif fidx > 1
        targfilename = sprintf('csc126%s.ncs', suffix);
        istarg = ismember(lower(dirfilenames), targfilename);
        if sum(istarg) == 1
            infiles{fidx} = dirfilenames{istarg};
        elseif sum(istarg) > 1
            error('cleanEyeLinkSession2:targfilename', ...
                'Multiple matches for targfile "%s"', targfilename);
        elseif sum(istarg) == 0
            error('cleanEyeLinkSession2:targfilename2', ...
                'There is no "csc126" with suffix "%s" in %s.', ...
                suffix, sessiondir);
        end
    else
        error('cleanEyeLinkSession:missingfile', ...
            'There is no %s in %s', filenames{fidx}, sessiondir);
    end
end

% Convert Nlx files <sessiondir> to .mat files in <tempdir>.
for fidx = 1:length(infiles)
    [~, basename] = fileparts(infiles{fidx});
    outbasename = strrep(basename, suffix, '');
    if fidx == 1
        % Events file; convert straight to output dir:
        fprintf('Reading events file...\n');
        tic;
        dg_Nlx2Mat(fullfile(sessiondir, infiles{fidx}), 'dest', outdir);
        if ~isequal(basename, outbasename)
            movefile( fullfile(outdir, [basename, '.mat']), ...
                fullfile(outdir, [outbasename, '.mat']) );
        end
        toc;
    else
        % Convert to .mat temp file and read into memory:
        tic;
        fprintf('Converting %s\n', infiles{fidx});
        dg_Nlx2Mat(fullfile(sessiondir, infiles{fidx}), 'dest', tempdir);
        toc;
        tempfilename = fullfile(tempdir, [basename '.mat']);
        tic;
        fprintf('Loading %s\n', basename);
        load(tempfilename); %#ok<LOAD>
        toc;
        % Estimate <Fs>, clean, save to output dir:
        framedurs = diff(dg_Nlx2Mat_Timestamps);
        Fs = 1e6 / (median(framedurs)/512);
        tic;
        fprintf('Cleaning %s\n', basename);
        dg_Nlx2Mat_Samples = dg_cleanEyeLink(dg_Nlx2Mat_Samples, Fs);
        toc;
        matfilenames{fidx - 1} = sprintf('%sclean.mat', outbasename);
        tic;
        fprintf('Saving %s\n', outbasename);
        save(fullfile(outdir, matfilenames{fidx - 1}), 'dg_Nlx2Mat_Timestamps', ...
            'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
        toc;
        delete(tempfilename);
    end
end

% Downsample .mat files in <tempdir> 32X into <outdir>.  Note that the
% 'reconcile' option is not invoked here because it CPU-heavy, so
% 'dg_downsample4LocalAvgRef' simply raises an error if there are timestamp
% conflicts or disordered timestamps.
tic;
dg_downsample4LocalAvgRef(matfilenames, outdir, 32, 'outdir', outdir, ...
    'filter', 0.4);
toc;

