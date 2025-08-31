function dg_subtractLocalAvgRefs(localavgrefs, sessiondir, prefix, varargin)
% Computes a complete set of local average referenced .mat files for
% <sessiondir>, using previously computed local average references. The
% results are saved to files named <infilename>-<refname>.mat, where
% <infilename> is the name of the input CSC file with the extension
% removed, and <refname> is the name of the reference file with the .mat
% extension removed. Output files are written to <sessiondir>.  Files that
% do not appear in <localavgrefs> are ignored.  Cells in <localavgrefs>
% that are empty produce no output files, but still affect the reference
% numbering.  If dg_Nlx2Mat_SamplesUnits in the reference file is 'V' and
% the input CSC file is in Neuralynx .ncs format, an attempt is made to
% convert the input CSC file to volts, and if the attempt fails an error is
% raised.  If the input CSC file is .mat format, and
% dg_Nlx2Mat_SamplesUnits in the reference file is not the same as in the
% input CSC file, an error is raised. Any pre-existing file with the same
% names as output files are overwritten.
%INPUTS
% localavgrefs: a cell vector of strings or cell arrays of strings where
%   each element is a legitimate value for the <CSCfilespec> argument
%   to dg_localAvgRef that specifies the files to average together and
%   save as a "local average reference" channel.  (Note this is identical to the
%   <tackle.localavgrefs> input to lfp_fishStory.)
% sessiondir: absolute or relative path to directory containing CSC files.
% prefix: string to use as local average reference filename prefix, as for
%   files created by dg_makeLocalAvgRefs.  Notes: if the first three
%   letters are csc or CSC, they are not included in the <refname> part of
%   the output filename. The file extension is assumed to be .mat.
%OUTPUTS are only to .mat files.
%OPTIONS
% 'srcdir', srcdir - the raw CSC files are in <srcdir>. Default value is
%   <sessiondir>.
%NOTES
% dg_downsample4LocalAvgRef.m Use this when you need to downsample only the
% files that will be used in your local avg reffing
% 
% dg_localAvgRefIfNeeded.m In bulk processing, you sometimes have a mix of
% sessions that have already been local avg reffed and others that haven't.
% Use this to skip the ones that are already done.
% 
% These are the three core local avg reffing functions:
% 
% dg_localAvgRef.m Calculates one local avg ref channel from the set of CSC
% files you specify.
% 
% dg_makeLocalAvgRefs.m Calls dg_localAvgRef for as many local avg refs you
% wish to define.
% 
% dg_subtractLocalAvgRefs.m This does the actual referencing, i.e. it
% subtracts the previously defined local avg ref signals from each of the
% CSC channels that went the avg, and produces a new local avg reffed file
% for each CSC channel.


%$Rev: 239 $
%$Date: 2016-03-10 19:30:16 -0500 (Thu, 10 Mar 2016) $
%$Author: dgibson $


srcdir = sessiondir;

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'srcdir'
                argnum = argnum + 1;
                srcdir = varargin{argnum};
            otherwise
                error('dg_subtractLocalAvgRefs:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('dg_subtractLocalAvgRefs:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

for refidx = 1:length(localavgrefs)
    CSCfilespec = localavgrefs{refidx};
    if isempty(CSCfilespec)
        continue
    end
    % load local avg ref
    refname = sprintf('%s%d.mat', prefix, refidx);
    reffilepath = fullfile(sessiondir, refname);
    load(reffilepath, '-mat');
    refsamples = dg_Nlx2Mat_Samples;
    refTS = dg_Nlx2Mat_Timestamps;
    refunits = dg_Nlx2Mat_SamplesUnits;
    clear dg_Nlx2Mat_Samples dg_Nlx2Mat_Timestamps;
    % get list of files to process
    if ischar(CSCfilespec)
        filelist = dir(fullfile(srcdir, CSCfilespec));
        CSCfilenames = {filelist.name};
    elseif iscell(CSCfilespec)
        CSCfilenames = CSCfilespec;
    else
        error('dg_subtractLocalAvgRefs:CSCfilespec', ...
            'Each cell in <localavgrefs> must contain a cell array or a string.');
    end
    for fileidx = 1:length(CSCfilenames)
        % subtract ref from file and save result
        [path, infilename, ext] = fileparts(CSCfilenames{fileidx}); %#ok<*ASGLU>
        switch ext
            case {'.ncs' '.NCS'}
                [TS, Samples, Hdr] = dg_readCSC(fullfile( srcdir, ...
                    CSCfilenames{fileidx} ));
                % <refTS> contains the reference list of reconciled time
                % stamps, whereas the raw files still need to be reconciled
                % to <refTS>.  We therefore need to delete any frames in
                % <TS>, <Samples> that do not exist in <refTS>,
                % <refsamples>.
                frames2del = find(~ismember(TS, refTS));
                if ~isempty(frames2del)
                    Samples(:, frames2del) = [];
                    TS(frames2del) = [];
                end
                dg_Nlx2Mat_Samples = reshape(Samples, size(refsamples));
                dg_Nlx2Mat_Timestamps = reshape(TS, size(refTS));
                dg_Nlx2Mat_SamplesUnits = 'AD';
                if isequal(refunits, 'V')
                    foundADBitVolts = false;
                    for k = 1:length(Hdr)
                        if regexp(Hdr{k}, '^\s*-ADBitVolts\s+')
                            ADBitVoltstr = regexprep(Hdr{k}, ...
                                '^\s*-ADBitVolts\s+', '');
                            ADBitVolts = str2double(ADBitVoltstr);
                            if isempty(ADBitVolts)
                                error('dg_subtractLocalAvgRefs:badADBitVolts', ...
                                    'Could not convert number from:\n%s', ...
                                    Hdr{k} );
                            end
                            dg_Nlx2Mat_Samples = ADBitVolts * Samples;
                            dg_Nlx2Mat_SamplesUnits = 'V';
                            foundADBitVolts = true; 
                            break
                        end
                    end
                    if ~foundADBitVolts
                        error('dg_subtractLocalAvgRefs:noADBitVolts', ...
                            'There is no voltage scale in the header of %s', ...
                            CSCfilenames{fileidx});
                    end
                end
            case {'.mat' '.MAT'}
                load(fullfile( srcdir, CSCfilenames{fileidx}), '-mat');
            otherwise
                error('dg_subtractLocalAvgRefs:unkext', ...
                    'Unknown filename extension: %s', ext);
        end
        if isequal(refunits, dg_Nlx2Mat_SamplesUnits)
            dg_Nlx2Mat_Samples = dg_Nlx2Mat_Samples - refsamples;
        else
            error('dg_subtractLocalAvgRefs:badunits', ...
                '%s has different units (%s) from ref (%s)', ...
                CSCfilenames{fileidx}, dg_Nlx2Mat_Samples, refunits);
        end
        if ~isempty(regexpi(refname, '^csc'))
            refname(1:3) = [];
        end
        outfilepath = fullfile(sessiondir, sprintf( ...
            '%s-%s', infilename, refname ));
        save(outfilepath, 'dg_Nlx2Mat_Samples', ...
            'dg_Nlx2Mat_SamplesUnits', 'dg_Nlx2Mat_Timestamps', '-v7.3');
    end
end

