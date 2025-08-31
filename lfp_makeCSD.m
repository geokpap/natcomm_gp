function lfp_makeCSD(datadir, filenames, varargin)
% Creates a series of CSD .mat files from an ordered list of files.
%INPUTS
% datadir: string containing the absolute or relative pathname to the
%   directory that contains the data files.
% filenames: cell array of strings, each of which contains the name of a
%   file in <datadir> to be used as input data.  The strings must be in the
%   same order as the physical contacts that the files represent.  A CSD
%   signal will be calculated for each consecutive triple of <filenames>.
%   The input files must be readable as a group by lfp_read2 (e.g.
%   Neuralynx format .ncs files, or in .mat files as created by
%   lfp_save_noGUI or dg_Nlx2Mat).
%OUTPUTS
% All outputs are stored in new files in <datadir>.  The CSD output
% filenames are of the form "CSD<n>.mat", where <n> starts at 1 and is
% incremented for each additional file.  The output file format is of the
% type created by lfp_save_noGUI or dg_Nlx2Mat.  In addition, a metadata
% file "lfp_makeCSD.mat" is created that contains variables <startTS>,
% <endTS>, <datadir>, and <filenames>.  <startTS>, <endTS> are timestamps
% at start and end of processing respectively.
%OPTIONS
% 'destdir', destdir - save output files to <destdir> instead of <datadir>.
%   If <destdir> is different from <datadir>, then any events files
%   (i.e. named events.* or *.evtsav) are copied to <destdir> from
%   <datadir>.
%NOTES
% The CSD formula used here comes from the discrete space version of the
% one-dimensional Laplacian operator, or second spatial derivative, of the
% voltage.  Given recordings from three uniformly spaced colinear recording
% sites a, b, and c, the first derivatives are (b-a)/d and (c-b)/d where d
% is the distance between the contacts.  The second derivative is thus
%   ((c-b)/d) - (b-a)/d) / d = (c + a - 2*b) / d^2
% There are additional constants of proportionality between the Laplacian
% of voltage and the current source density, including the resistivity of
% the medium, which is typically not measured.  We therefore abandon all
% hope of getting CSD in physically meaningful units like microamps per
% cubic millimeter, and simply compute c + a - 2*b and designate the units
% as 'arbs'.
%   The actual formula in continuous space is CSD = -sigma * Del^2(V) where
% sigma is the electrical conductivity of the tissue, V is the voltage
% field, and Del is the gradient operator.  Consequently, the default
% formula used here is proportional to the negative of the current source
% density.  That means that if you run it on recordings that were recorded
% with inverted polarity (the Neuralynx default), the result is the
% non-inverted CSD.  See Bedard & Destexhe, PHYSICAL REVIEW E 84, 041909
% (2011); and ?Electric Fields Of The Brain?, PL Nunez & R Srinivasan
% (2006) Oxford University Press.

%$Rev: 402 $
%$Date: 2019-08-20 16:50:18 -0400 (Tue, 20 Aug 2019) $
%$Author: dgibson $

global lfp_ActiveFilenums lfp_Samples lfp_SamplesPerFrame

destdir = datadir;
startTS = now(); %#ok<NASGU>

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
        case 'destdir'
            argnum = argnum + 1;
            destdir = varargin{argnum};
        otherwise
            error('lfp_makeCSD:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

if ~isequal(destdir, datadir)
    if ~exist(destdir, 'dir')
        [p, n] = fileparts(destdir);
        if ~exist(p, 'dir')
            error('datadir:mkdir', ...
                'Directory auto-creation is not recursive: %s', destdir);
        end
        mkdir(p, n);
    end
    evtfiles = dir(fullfile(datadir, 'events.*'));
    evtfnames = {evtfiles.name};
    evtfiles = dir(fullfile(datadir, '*.evtsav'));
    evtfnames = unique([evtfnames {evtfiles.name}]);
    for k = 1:length(evtfnames)
        dg_copyfile( fullfile(datadir, evtfnames{k}), ...
            fullfile(destdir, evtfnames{k}) );
    end
end

lfp_read2('free', 'preset', datadir, {'', filenames{1:3}});
a = lfp_Samples{lfp_ActiveFilenums(1)}(:);
b = lfp_Samples{lfp_ActiveFilenums(2)}(:);
c = lfp_Samples{lfp_ActiveFilenums(3)}(:);
for fidx = 3:length(filenames)
    CSD = reshape(c + a - 2*b, lfp_SamplesPerFrame, []);
    CSDname = sprintf('CSD%d', fidx - 2);
    lfp_createWave(@lfp_waverecord, 1, CSD, 'name', CSDname);
    % save CSD
    lfp_save_noGUI(CSDname, 'csc', 'destdir', destdir);
    if fidx < length(filenames)
        % rotate last two files, read new one
        clear CSD;
        lfp_Samples = {};
        a = b;
        b = c;
        lfp_read2('free', 'preset', datadir, {'', filenames{fidx+1}});
        if numel(lfp_Samples{lfp_ActiveFilenums(1)}) ~= numel(b)
            error( 'lfp_makeCSD:numsamp', ...
                '%s has a different number of samples from %s', ...
                filenames(fidx+1), filenames(fidx) );
        end
        c = lfp_Samples{lfp_ActiveFilenums(1)}(:);
    end
end

endTS = now(); %#ok<NASGU>
save( fullfile(destdir, 'lfp_makeCSD.mat'), 'startTS', 'endTS', ...
    'datadir', 'filenames' );

