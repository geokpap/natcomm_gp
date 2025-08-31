function dg_downsample36(src, dest, N, varargin)
% Read CSC file from src, write to dest after filtering and downsampling by
% a factor of N.  The first sample remains the first sample after
% downsampling, but there is no attempt to reduce turn-on and turn-off
% transients in the filtering.  Frame size is kept the same, and the end is
% padded with zeros if needed to fill the last frame.
%   Filter is butter(4, fc, 'low'), where fc = 0.9/N, run forward and
% backward to yield 8th-order filtering with zero phase shift.
% The <src> and <dest> file formats are assumed to be Neuralynx .ncs unless
% the file name extension is .mat, in which case they are assumed to be the
% format written by dg_Nlx2Mat.
%OPTIONS
% 'int16' - writes output .mat file with sample data in int16
%   representation.  Note that this will result in a file of all zeros if
%   typical physiological data are read from a Nlx file without specifying
%   'nounits'.
% 'nounits' - leaves sample data in unscaled (i.e. raw) A/D units when
%   reading Nlx format.
% 'overwrite' - raises a warning instead of an error and overwrites the
%   existing <dest> file if there is one.
% 'skip' - raises a warning instead of an error and skips the
%   existing <dest> file if there is one.
% 'verbose' - this is an option which optionally can be specified as an
%   option in the argument list when calling the function from the command
%   line or from a file or in any other manner.
%NOTES
%   This is a copy of dg_downsample as it was in Revision 36 of dg_lib,
% before I added the 'reconcile' option, which brought with it a whole
% cascade of subtle bugs in the batch-processing mode which is required for
% handling extremely large files.  I have backed it out of the version
% control system here for emergency use in case any *more* such bugs
% surface and somebody just needs to downsample one file.
%   The rest of these notes are from the original Rev 36 version.
%
%   No attempt is made to set ADBitVolts (or any other Nlx header
% parameter) correctly when producing a Neuralynx format file from a .mat
% file, but that's only due to laziness on the part of the programmer. :-/
%   A limit related to the size of the largest contiguous block of
% available memory is imposed on the amount of data slurped in at one time
% when reading from Neuralynx files.  If it is necessary to process the
% file in multiple slurps, then to minimize the startup transients on the
% slurp boundaries, successive slurps are made to overlap by two frames,
% and the beginning/end frames at the boundary are thrown away after
% filtering.

%$Rev: 116 $
%$Date: 2011-05-26 19:09:07 -0400 (Thu, 26 May 2011) $
%$Author: dgibson $

int16flag = false;
nounitsflag = false;
overwriteflag = false;
skipflag = false;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'int16'
            int16flag = true;
        case 'nounits'
            nounitsflag = true;
        case 'overwrite'
            overwriteflag = true;
        case 'skip'
            skipflag = true;
        case 'verbose'
            verboseflag = true;
        otherwise
            error('dg_downsample:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~exist(src, 'file')
    error('dg_downsample:nofile', ...
        'The input file %s does not exist', src);
end
if exist(dest, 'file')
    if overwriteflag
        warning('dg_downsample:fileexists', ...
            'The output file %s already exists; overwriting', src);
    elseif skipflag
        warning('dg_downsample:fileexists2', ...
            'The output file %s already exists; skipping', src);
        return
    else
        error('dg_downsample:fileexists', ...
            'The output file %s already exists', src);
    end
end
if verboseflag
    str = (sprintf('dg_downsample(%s, %s, %d', src, dest, N));
    for k = 1:length(varargin)
        str = [str ', ' dg_thing2str(varargin{k})];
    end
    str = [str ')'];
    disp(str);
end

% Set up constant things:
[srcdir, srcname, srcext] = fileparts(src);
[destdir, destname, destext] = fileparts(dest);
fc = 0.9/N;
[z, p, k] = butter(4, fc, 'low');
[sos,g]=zp2sos(z,p,k);
h2=dfilt.df2sos(sos,g);
tempfilenum = 1;
tempbase = sprintf('dg_downsample_temp%d_', tempfilenum);
while ~isempty(dir([tempbase '*']))
    tempfilenum = tempfilenum + 1;
    tempbase = sprintf('dg_downsample_temp%d_', tempfilenum);
end

% Read and process:
if verboseflag
    disp(sprintf('Reading file %s', src));
end
if isequal(lower(srcext), '.mat')
    load('-MAT', src);
    ts = dg_Nlx2Mat_Timestamps;
    clear dg_Nlx2Mat_Timestamps;
    samples = dg_Nlx2Mat_Samples;
    clear dg_Nlx2Mat_Samples;
    if exist('dg_Nlx2Mat_SamplesUnits', 'var')
        units = dg_Nlx2Mat_SamplesUnits;
    else
        units = 'arbs';
    end
    dg_downsample_filtfilt;
else
    % For 512-sample frames, it's 4128 bytes per record:
    recsize = 4128;
    userview = memory;
    bytesinblock = userview.MaxPossibleArrayBytes;
    [ch, header] = Nlx2MatCSC(src, [0, 1, 0, 0, 0], 1, 1);
    numrecs = length(ch);
    clear ch;
    startrec = 1;
    if bytesinblock < 24
        error('dg_downsample:outofmem', ...
            'There are only %d bytes available; pack memory', bytesinblock);
    end
    maxbytes2use = bytesinblock/10;
    % calculate endrec and numbatches:
    if maxbytes2use < numrecs * recsize
        batchsize = fix(maxbytes2use/recsize);
        endrec = batchsize;
        % There will actually be slightly more than ceil(numrecs/batchsize)
        % batches because two frames (records) are thrown away at each
        % batch boundary.
        % Analysis of extreme cases: 
        %   batchsize = 2 makes the overlap impossible, so therefore we
        %       must always have batchsize > 2.
        %   maxbytes2use = numrecs * recsize; this triggers the single
        %   batch method
        %   maxbytes2use = numrecs * recsize - 1
        %       batchsize = 3 => numbatches = numrecs - 2
        %       recsize = 1 =>  batchsize = maxbytes2use =>
        %           numbatches = 2 + ceil(-(numrecs-4)/(numrecs-3)), or
        %           numbatches = 2 - floor((numrecs-4)/(numrecs-3))
        %           numbatches = 2 - 0 for numrecs > 4
        %           numbatches = 2 - 0 for numrecs = 4
        %           impossible for numrecs = 3, i.e. for maxbytes2use = 2
        % That last case determines the minimum bytesinblock required.
        numbatches = 2 + ceil((numrecs - 2*(batchsize-1))/(batchsize-2));
    else
        endrec = numrecs;
        numbatches = 1;
    end
    if verboseflag
        disp(sprintf('Reading in %d batches', numbatches));
    end
    done = false;
    numsubsamples = 0;
    while ~done
        if verboseflag
            disp(sprintf( ...
                'Reading frames %.0f-%.0f', ...
                startrec, endrec ));
        end
        % To ensure maximum confusion, Nlx2MatCSC uses 0-based array
        % indexing, not 1-based indexing:
        [ts, ch, fs, nv, samples] ...
            = Nlx2MatCSC( src, [1, 1, 1, 1, 1], 0, 2, ...
            [startrec-1 endrec-1] );
        if isequal(lower(destext), '.mat') 
            if nounitsflag
                units = 'AD';
            else
                % Convert from AD units to volts if possible
                if exist('ADBitVolts', 'var')
                    if ~isempty(ADBitVolts)
                        samples = ADBitVolts ...
                            * samples;
                        units = 'V';
                    end
                else
                    for k = 1:length(header)
                        units = 'AD';
                        if regexp(header{k}, '^\s*-ADBitVolts\s+')
                            ADBitVoltstr = regexprep(header{k}, ...
                                '^\s*-ADBitVolts\s+', '');
                            ADBitVolts = str2double(ADBitVoltstr);
                            if isempty(ADBitVolts)
                                warning('dg_downsample:badADBitVolts', ...
                                    'Could not convert number from:\n%s', ...
                                    header{k} );
                            else
                                samples = ADBitVolts ...
                                    * samples;
                                units = 'V';
                            end
                        end
                    end
                end
            end
        end
        dg_downsample_filtfilt;
        if startrec == 1
            % compute sample period <sampleT>, initial sample indices:
            framesize = size(samples,1);
            framedurs = ts(2:end)-ts(1:end-1);
            minfd = min(framedurs);
            framedurs = framedurs(framedurs < 2*minfd);
            medianfd = median(framedurs);
            sampleT = medianfd/framesize;
            startsubsample = 1;     % absolute sample index in src file
            firstsampinbatch = 1;   % absolute sample index in src file
        end
        % compute <endsubsample>, absolute index in src file of the
        % range of samples to subsample:
        if endrec == numrecs
            % last batch:
            endsubsample = numrecs * framesize;
        else
            % omit last frame in inner batches:
            endsubsample = (size(samples,2) - 1) * framesize ...
                + firstsampinbatch - 1;
        end
        if numbatches == 1
            done = true;
        else
            % downsample the current batch, picking up from where previous
            % batch left off, and save in file in cwd:
            allTS = repmat(ts, framesize, 1) + ...
                repmat((0 : framesize - 1)' * sampleT, 1, length(ts));
            if verboseflag
                disp(sprintf( ...
                    'Downsampling %dX, TS %.0f-%.0f', ...
                     N, startrec, endrec));
            end
            subsampleidx = (startsubsample : N : endsubsample);
            tempsamples = samples(subsampleidx - firstsampinbatch + 1);
            tempallTS = allTS(subsampleidx - firstsampinbatch + 1);
            tempfilename = sprintf('%s%07d.mat', tempbase, startrec);
            if verboseflag
                disp(sprintf( ...
                    'Saving file %s, TS %.0f-%.0f', ...
                     tempfilename, tempallTS(1), tempallTS(end)));
            end
            if int16flag
                tempsamples = int16(tempsamples);
            end
            save(tempfilename, 'tempsamples', 'tempallTS');
            numsubsamples = numsubsamples + numel(tempsamples);
            if endrec == numrecs
                done = true;
            else
                % update indices for next batch:
                startsubsample = subsampleidx(end) + N;
                startrec = endrec - 1;  % repeat last two frames
                firstsampinbatch = framesize * (startrec - 1) + 1;
                endrec = min(startrec + batchsize - 1, numrecs);
            end
        end
    end
    if numbatches > 1
        % Reconstruct entire signal with timestamps into row vectors from
        % batches saved in tempfiles, deleting each file after reading:
        tempfiles = dir([tempbase '*.mat']);
        samples = NaN(1, numsubsamples);
        allTS = NaN(1, numsubsamples);
        lastsampidx = 0;
        for k = 1:length(tempfiles)
            load(tempfiles(k).name);
            samplerange = lastsampidx + (1:length(tempsamples));
            samples(samplerange) = tempsamples;
            allTS(samplerange) = tempallTS;
            lastsampidx = lastsampidx + length(tempsamples);
            delete(tempfiles(k).name);
        end
    end
end

if isequal(lower(srcext), '.mat') || numbatches == 1
    % Entire file is already in memory; process it.
    % compute sample period <sampleT>:
    framesize = size(samples,1);
    framedurs = ts(2:end)-ts(1:end-1);
    minfd = min(framedurs);
    framedurs = framedurs(framedurs < 2*minfd);
    medianfd = median(framedurs);
    sampleT = medianfd/framesize;
    % downsample into row vector:
    allTS = repmat(ts, framesize, 1) + ...
        repmat((0 : framesize - 1)' * sampleT, 1, length(ts));
    if verboseflag
        disp(sprintf('Downsampling %dX', N));
    end
    samples = reshape(samples(1:N:end), 1, []);
    allTS = reshape(allTS(1:N:end), 1, []);
end

% pad if needed and reshape into frames:
extrasamples = rem(numel(samples), framesize);
padding = [];
if extrasamples ~= 0
    warning( 'dg_downsample:padding', ...
        'padded "%s" with zeros to fill last frame', destname );
    padding = zeros(1, framesize-extrasamples);
end
samples = reshape([samples padding], framesize, []);
allTS = reshape([allTS padding], framesize, []);
ts = allTS(1,:);

% Create output file:
if isequal(lower(destext), '.mat')
    dg_Nlx2Mat_Timestamps = ts;
    dg_Nlx2Mat_Samples = samples;
    dg_Nlx2Mat_SamplesUnits = units;
    if int16flag
        dg_Nlx2Mat_Samples = int16(dg_Nlx2Mat_Samples);
    end
    save(dest, ...
        'dg_Nlx2Mat_Timestamps', 'dg_Nlx2Mat_Samples', ...
        'dg_Nlx2Mat_SamplesUnits' );
else
    if isequal(lower(srcext), '.mat')
        maxsamp = max(samples(:));
        minsamp = min(samples(:));
        if maxsamp >= -minsamp * 2047/2048
            samples = round(2047 * samples / maxsamp);
        else
            samples = round(2048 * samples / -minsamp);
        end
        ch = zeros(size(ts));
        samplerate = framesize / medianfd;
        fs = repmat(samplerate, size(ts));
        nv = repmat(framesize, size(ts));
        header = {'######## Neuralynx Data File Header'
            sprintf('## Converted by dg_downsample from: %s ', src)
            sprintf( '## Time Opened (dd-mmm-yyyy HH:MM:SS): %s ', ...
            datestr(now, 0) )};
    else
        ch = repmat(ch(1), size(ts));
        fs = repmat(fs(1), size(ts));
        nv = repmat(framesize, size(ts));
    end
    if verboseflag
        disp(sprintf('Writing file %s', dest));
    end
    Mat2NlxCSC_411(dest, 0, 1, 1, length(ts), ones(1,6), ...
        ts, ch, fs, nv, samples, header);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dg_downsample_filtfilt
        samplesize = size(samples);
        if verboseflag
            disp(sprintf('Filtering first time'));
        end
        samples = filter(h2, samples(end:-1:1));
        if verboseflag
            disp(sprintf('Filtering second time'));
        end
        samples = reshape(filter(h2, samples(end:-1:1)), samplesize);
    end

end
