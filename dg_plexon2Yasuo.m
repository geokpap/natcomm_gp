function dg_plexon2Yasuo(datadir)
% Copies and renames Multiple Cut Cluster files created by Plexon so that
% they adhere to the Yasuo naming conventions, i.e. they reside in the
% animal directory and have names of the form <sessionID>-T<tnum>.mat (e.g.
% "acq10-T01.mat"), where <tnum> is the tetrode number as saved in the .mat
% file, regardless of the NAME of the .mat file.  Identifies session
% directories on the basis of their containing *.ntt files, and loads every
% .mat file in that directory to find cluster data.    Normally, the animal
% directory is assumed to be the parent directory of <datadir>.  However,
% if the directory name is 'onlytrialcut', then the animal directory is
% assumed to be the grandparent directory, and if the directory name is
% 'ccw', then the animal directory is assumed to be the great-grandparent
% directory. The session ID is the name of the directory immediately after
% the animal directory in <datadir>. If a file already exists with the name
% of an output file, the old file is renamed to <filename>.old<n> where <n>
% is the smallest counting number that does not result in a name collision.
% Recursively runs on every subdirectory of <datadir>.  Skips any empty
% cluster arrays found in .mat files and issues a warning.  The cluster
% data may contain waveforms as well as timestamps and cluster IDs.

%$Rev: 157 $
%$Date: 2012-09-18 21:58:14 -0400 (Tue, 18 Sep 2012) $
%$Author: dgibson $

isSessionDir = false;
matfilenames = {};
files = dir(datadir);
for k = 1:length(files)
    if ~ismember(files(k).name, {'.' '..'})
        if files(k).isdir
            dg_plexon2Yasuo(fullfile(datadir, files(k).name));
        else
            [p,n,x] = fileparts(files(k).name); %#ok<*ASGLU>
            if isequal(lower(x), '.mat')
                matfilenames{end+1} = files(k).name; %#ok<AGROW>
            elseif isequal(lower(x), '.ntt')
                isSessionDir = true;
            end
        end
    end
end

if isSessionDir
    [p,n] = fileparts(datadir);
    switch n
        case 'onlytrialcut'
            [p, sessionID] = fileparts(p);
        case 'ccw'
            p = fileparts(p);
            [p, sessionID] = fileparts(p);
        otherwise
            sessionID = n;
    end
    animaldir = p;
    % Check every .mat file for cluster data, which is identified by a
    % variable whose name is of the for tt<n> where <n> is the tetrode
    % number.
    for matidx = 1:length(matfilenames)
        S = load(fullfile(datadir, matfilenames{matidx}));
        flds = fieldnames(S);
        for fldidx = 1:length(flds)
            toks = regexpi(flds{fldidx}, '^tt(\d+)$', 'tokens');
            if ~isempty(toks) && ~isempty(S.(flds{fldidx}))
                ttnum = str2double(toks{1});
                outfilename = sprintf('%s-T%02d.mat', ...
                    sessionID, ttnum);
                if exist(fullfile(animaldir, outfilename), 'file')
                    % Name collision!
                    n = 1;
                    while exist(fullfile(animaldir, sprintf('%s.old%d', ...
                            outfilename, n )), 'file')
                        n= n + 1;
                    end
                    movefile(fullfile(animaldir, outfilename), ...
                        fullfile(animaldir, sprintf('%s.old%d', ...
                            outfilename, n )));
                end
                save(fullfile(animaldir, outfilename), '-struct', 'S', ...
                    flds{fldidx});
            end
        end
    end
end

