function dg_localAvgRef(sessiondir, CSCfilespec, outfilepath)
% Using a specified set of CSC files in <sessiondir>, constructs a local
% average reference signal from them.  Non-MAT files are assumed to be
% Neuralynx files, and are converted to volts before averaging.  MAT files
% must contain dg_Nlx2Mat_SamplesUnits, and its value must be 'V'.
% Reconciles timestamps across all files before doing the computation.
% If there is a pre-existing file with the same name as the output file, it
% is overwritten.
%INPUTS
% sessiondir: string containing relative or absolute path of directory
%   containing CSC files.
% CSCfilespec: string or cell string array.  If string, it is interpreted
%   as a a system-appropriate expression for matching to the file name.  If
%   cell, it is interpreted as a vector of strings explicitly containing
%   file names including the extension.
% outfilepath: a string containing a relative or absolute file name to
%   which to write the timestamp-reconciled result.  If extension is
%   '.mat', writes in dg_Nlx2Mat format.  Otherwise writes in Neuralynx
%   format.  The output directory must already exist or else writing the
%   output file fails.
%OUTPUT is only to the file at <outfilepath>.
%NOTES
% See also dg_makeLocalAvgRefs, dg_subtractLocalAvgRefs,
% dg_downsample4LocalAvgRef.

%$Rev: 211 $
%$Date: 2015-03-12 21:08:04 -0400 (Thu, 12 Mar 2015) $
%$Author: dgibson $

if ischar(CSCfilespec)
    filelist = dir(fullfile(sessiondir, CSCfilespec));
    CSCfilenames = {filelist.name};
elseif iscell(CSCfilespec)
    CSCfilenames = CSCfilespec;
else
    error('dg_localAvgRef:CSCfilespec', ...
        '<CSCfilespec> must be a cell array or a string.');
end
if isempty(CSCfilenames)
    error('dg_localAvgRef:CSCfilespec2', ...
        '<CSCfilespec> does not match any file.');
end
framesize = 512;
numInSum = 0;

% Since there is no way to avoid reading the samples from a .MAT file, we
% do the timestamps reconciliation and the running sum in a single somewhat
% convoluted loop.  Running sum is computed in units of Volts.
for filenum = 1:length(CSCfilenames)
    fullfilename = fullfile(sessiondir, CSCfilenames{filenum});
    if ~exist(fullfilename, 'file')
        error('dg_localAvgRef:filename', ...
            'The file %s does not exist', fullfilename);
    end
    [p, n, ext] = fileparts(CSCfilenames{filenum}); %#ok<*ASGLU>
    if isequal(upper(ext), '.MAT')
        load(fullfilename, '-mat');
        if ~isequal(dg_Nlx2Mat_SamplesUnits, 'V') %#ok<NODEF>
            error('dg_localAvgRef:units', ...
                '%s has units %s; must be in units of V.', ...
                fullfilename, dg_Nlx2Mat_SamplesUnits);
        end
        if filenum == 1
            reconcTS = dg_Nlx2Mat_Timestamps; %#ok<NODEF>
            runningsum = dg_Nlx2Mat_Samples;
        else
            frames2delfromthis = ...
                ~ismember(dg_Nlx2Mat_Timestamps, reconcTS); %#ok<NODEF>
            frames2delfromsum = ...
                ~ismember(reconcTS, dg_Nlx2Mat_Timestamps);
            reconcTS = intersect(reconcTS, dg_Nlx2Mat_Timestamps);
            dg_Nlx2Mat_Samples = ...
                reshape(dg_Nlx2Mat_Samples, framesize, []);
            dg_Nlx2Mat_Samples(:, frames2delfromthis) = [];
            runningsum(:, frames2delfromsum) = [];
            runningsum = runningsum + dg_Nlx2Mat_Samples;
        end
        clear dg_Nlx2Mat_SamplesUnits dg_Nlx2Mat_Timestamps dg_Nlx2Mat_Samples
    else
        [TS, Samples, Hdr] = dg_readCSC(fullfilename);
        foundADBitVolts = false;
        for k = 1:length(Hdr)
            if regexp(Hdr{k}, '^\s*-ADBitVolts\s+')
                ADBitVoltstr = regexprep(Hdr{k}, ...
                    '^\s*-ADBitVolts\s+', '');
                ADBitVolts = str2double(ADBitVoltstr);
                if isempty(ADBitVolts)
                    error('dg_localAvgRef:badADBitVolts', ...
                        'Could not convert number from:\n%s', ...
                        Hdr{k} );
                else
                    Samples = ADBitVolts * Samples;
                end
                foundADBitVolts = true;
                break
            end
        end
        if ~foundADBitVolts
            error('dg_localAvgRef:noADBitVolts', ...
                'There is no voltage scale in the header of %s', ...
                fullfilename);
        end
        if filenum == 1
            reconcTS = TS;
            runningsum = Samples;
        else
            frames2delfromthis = ...
                ~ismember(TS, reconcTS);
            frames2delfromsum = ...
                ~ismember(reconcTS, TS);
            reconcTS = intersect(reconcTS, TS);
            Samples(:, frames2delfromthis) = [];
            runningsum(:, frames2delfromsum) = [];
            runningsum = runningsum + Samples;
        end
        clear TS Samples
    end
    numInSum = numInSum + 1;
end

[p, n, ext] = fileparts(outfilepath);
dg_Nlx2Mat_Samples = runningsum / numInSum; %#ok<NASGU>
clear runningsum
if isequal(upper(ext), '.MAT')
    dg_Nlx2Mat_SamplesUnits = 'V'; %#ok<NASGU>
    dg_Nlx2Mat_Timestamps = reconcTS; %#ok<NASGU>
    save(outfilepath, 'dg_Nlx2Mat_Timestamps', ...
            'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
else
    % Since Nlx files are strictly integer, we must scale the output and
    % put a scale factor in the header.
    posbitvolts = max(dg_Nlx2Mat_Samples(:)) / 2047;
    negbitvolts = min(dg_Nlx2Mat_Samples(:)) / -2048;
    bitvolts = max(posbitvolts, negbitvolts);
    Hdr = {'######## Neuralynx Data File Header '
        sprintf('## File Name: %s ', outfilepath)
        sprintf('## Time Opened: %s ', datestr(now))
        '## written by dg_localAvgRef '
        ' '
        sprintf('-ADBitVolts %s ', bitvolts);
        ' '
        };
    dg_writeCSC(outfilepath, reconcTS, ...
        dg_Nlx2Mat_Samples / bitvolts, Hdr);
end

