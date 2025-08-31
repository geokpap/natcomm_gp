function dg_patchDroppedFrames(sessiondir, varargin)
% Constructs a master list of timestamps containing the union of the
% timestamps across all files in <sessiondir>, and creates a new set of
% *.mat files in dg_Nlx2Mat (aka lfp_lib) format that contain frames full
% of a dummy value where the missing frames should have been.
%INPUTS
% sessiondir: a directory containing *.ncs (Neuralynx CSC) files.
%OUTPUTS
% All output goes to files.
%OPTIONS
% 'destdir', destdir - Creates output files in <destdir>.  Default value is
%   <sessiondir>.
% 'dummyval', val - <val> is a floating point value to use to fill missing
%   frames.  Default value is NaN.
% 'verbose': guess.

%$Rev: 243 $
%$Date: 2016-06-24 21:51:26 -0400 (Fri, 24 Jun 2016) $
%$Author: dgibson $

destdir = sessiondir;
framesize = 512;
val = NaN;
verboseflag = false;

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
        case 'dummyval'
            argnum = argnum + 1;
            val = varargin{argnum};
        case 'destdir'
            argnum = argnum + 1;
            destdir = varargin{argnum};
        case 'verbose'
            verboseflag = true;
        otherwise
            error('funcname:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

cscfiles = dir(fullfile(sessiondir, '*.ncs'));
cscfilenames = {cscfiles.name};

% Collect timestamps from all files:
if verboseflag
    fprintf('Constructing master TS set...\n');
end
TS = cell(size(cscfilenames));
for fidx = 1:length(cscfilenames)
    filepath = fullfile(sessiondir, cscfilenames{fidx});
    [TS{fidx}] = dg_readCSC(filepath);
    if verboseflag
        fprintf('Read timestamps from %s\n', cscfilenames{fidx});
    end
    if any(diff(TS{fidx} >=0 ))
        error('dg_patchDroppedFrames:TSoutoforder', ...
            '%s contains disordered timestamps.', filepath);
    end
end
% Remove empty channels:
fnames2del = cellfun(@isempty, TS);
TS(fnames2del) = [];
cscfilenames(fnames2del) = [];
% Under mysterious conditions, the Neuralynx hardware can produce files
% that have timestamps that are off by roughly half a frame.  All of the
% files have timestamp problems of one kind or another under such
% conditions, so it necessary to resample them all with padding inserted as
% needed.
%   First make sure that the timestamps are equivalent to 32.6 kHz
% nominal sample rate for each file:
for fidx = 1:length(TS)
    framedur = median(diff(TS{fidx}));
    if abs(framedur - 512e6 / 32.6e3) > 1e3
        warning('dg_patchDroppedFrames:rate', ...
            '%s has Fs = %.2g', cscfilenames{fidx}, 512e6/framedur);
    end
end
%   Then resample each file into a common time framework consisting of
% equally spaced frames at the true sample rate <Ts> determined from frames
% that have typical durations.  Start by finding the typical durations
% <framedurvals>:
allframedurs = diff(cell2mat(TS));
allframedurs(allframedurs<0) = [];
framedurvals = unique(allframedurs);
framedurvals(framedurvals > 1.1 * framesize * 1e6 / 32.6e3) = [];
% Compute sample period <Ts> weighted by the number of frames of each
% typical duration:
numdurs = NaN(size(framedurvals));
sampdurs = framedurvals/framesize; % all times are in usec
for fdvix = 1:length(framedurvals)
    numdurs(fdvix) = sum(allframedurs == framedurvals(fdvix));
end
Ts = sum(sampdurs .* numdurs) / sum(numdurs);
fd = Ts*framesize; % average frame duration
% and find the starting and ending timestamps:
startTS = min(cellfun(@(a) a(1), TS(~cellfun(@isempty, TS))));
endTS = max(cellfun(@(a) a(end), TS(~cellfun(@isempty, TS))));
% Construct the canonical analytic timestamp vector.
numframes = ceil((endTS - startTS) / fd);
allTS = round(startTS + (0 : numframes) * fd);

% Create patched, resampled output files:
for fidx = 1:length(cscfilenames)
    if verboseflag
        fprintf('Reading samples from %s\n', cscfilenames{fidx});
    end
    [~, Samples, header] = dg_readCSC(fullfile( sessiondir, ...
        cscfilenames{fidx} ));
    if isempty(TS{fidx})
        continue
    end
    % Get scale factor:
    for k = 1:length(header)
        if regexp(header{k}, '^\s*-ADBitVolts\s+')
            ADBitVoltstr = regexprep(header{k}, ...
                '^\s*-ADBitVolts\s+', '');
            ADBitVolts = str2double(ADBitVoltstr);
        end
    end
    % Construct output file in dg_Nlx2Mat_Samples:
    dg_Nlx2Mat_Samples = NaN(framesize, length(allTS));
    % In an uncorrupted file, TS{fidx} is monotonic increasing, as is <TS>.
    % We therefore use a pointer to the first element of allTS that is >=
    % TS{fidx}(frmidx) to bracket the search region: 
    ptr = find(allTS >= TS{fidx}(1), 1); % 
    % <frmidx> points at the frame in TS we are seeking.
    for frmidx = 1:length(TS{fidx})
        if verboseflag && mod(frmidx, 10000) == 0
            fprintf('Frame %d/%d\n', frmidx, length(TS{fidx}));
        end
        if ptr == 1 && TS{fidx}(frmidx) == allTS(ptr)
            allTSidx = ptr;
        else
            while ptr <= length(allTS) && TS{fidx}(frmidx) > allTS(ptr)
                ptr = ptr + 1;
            end
            if allTS(ptr) - TS{fidx}(frmidx) < ...
                    TS{fidx}(frmidx) - allTS(ptr - 1)
                allTSidx = ptr;
            else
                allTSidx = ptr - 1;
            end
        end
        TSerr = allTS(allTSidx) - TS{fidx}(frmidx);
        if abs(TSerr) < Ts/2
            % copy frame verbatim
            dg_Nlx2Mat_Samples(:, allTSidx) = Samples(:, frmidx);
        else
            % shift samples as required; positive <TSerr> means the matched
            % element in allTS is greater than the timestamp in the file.
            errpts = round(TSerr / Ts);
            allTSidx2 = framesize * (allTSidx - 1);
            dg_Nlx2Mat_Samples( (allTSidx - 1) * framesize ...
                - errpts + (1:512) ) = Samples(:, frmidx);
        end
    end
    % Substitute different default value if necessary
    if ~isnan(val)
        dg_Nlx2Mat_Samples(isnan(dg_Nlx2Mat_Samples)) = val;
    end
    % Create output file:
    clear('Samples');
    if isempty(ADBitVolts)
        warning('dg_Nlx2Mat:badADBitVolts', ...
            'Could not convert number from:\n%s', ...
            header{k} );
        dg_Nlx2Mat_SamplesUnits = 'AD'; %#ok<NASGU>
    else
        dg_Nlx2Mat_Samples = ADBitVolts ...
            * dg_Nlx2Mat_Samples; %#ok<NASGU>
        dg_Nlx2Mat_SamplesUnits = 'V'; %#ok<NASGU>
    end
    dg_Nlx2Mat_Timestamps = allTS; %#ok<NASGU>
    [~, name] = fileparts(cscfilenames{fidx});
    save( fullfile(destdir, [name '.mat']), ...
        'dg_Nlx2Mat_Timestamps', ...
        'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3' );
    clear('dg_Nlx2Mat_Timestamps', ...
        'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits');
end

