function dg_vertAlignSession(sessiondir, varargin)
%INPUTS
% sessiondir: absolute or relative pathname to a directory containing
%   Neuralynx format *.nse files from which to make vertically aligned
%   versions using 'dg_vertAlignSpike'.
%OUTPUTS
% All output is to new *.nse files in <sessiondir>, named '*_va.nse'
%   where '*' denotes the same prefix that it does in the input *.nse
%   files.
%OPTIONS
% 'outdir', outdir - creates output files in <outdir>.  Default is
%   <sessiondir>.
% 'wrap' - creates a subdirectory in <outdir> for each output file so that
%   each output file is by itself in its own subdirectory.  The name of
%   each subdirectory is the filename prefix, denoted above as '*'.
% 'overwrite' - overwrites any *.nse that have the same names as files that
%   would normally be created by this function.  Default is to leave the
%   existing files intact.
%NOTES
% Uses 'dg_findRealSessionFiles' to skip bogus *.nse files.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

outdir = sessiondir;
overwriteflag = false;
wrapflag = false;
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
        case 'outdir'
            argnum = argnum + 1;
            outdir = varargin{argnum};
        case 'overwrite'
            overwriteflag = true;
        case 'wrap'
            wrapflag = true;
        otherwise
            error('afp_Mat2SE:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

suffix = dg_findRealSessionFiles(sessiondir);
files = dir(fullfile(sessiondir, ['*' suffix '.nse']));
for fidx = 1:length(files)
    filename = files(fidx).name;
    [~, name] = fileparts(filename);
    fullfilename = fullfile(sessiondir, filename);
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    if files(fidx).bytes == 16384
        warning('dg_vertAlignSession:empty', ...
            'File %s is empty.', filename);
        continue
    end
    if wrapflag
        destdir = fullfile(outdir, name);
        if ~exist(destdir, 'dir');
            mkdir(destdir);
        end
    else
        destdir = outdir;
    end
    outfilename = sprintf('%s_va.nse', name);
    if ~overwriteflag
        if exist(fullfile(destdir, outfilename), 'file')
            % Output file already exists; skip.
            continue
        end
    end
    dg_vertAlignSpike(fullfilename, fullfile(destdir, outfilename));
end


