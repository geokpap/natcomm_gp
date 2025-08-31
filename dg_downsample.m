function dg_downsample(src, dest, N, varargin)
% Read CSC file from src, write to dest after filtering and downsampling by
% a factor of N.  The first sample remains the first sample after
% downsampling, but there is no attempt to reduce turn-on and turn-off
% transients in the filtering.  Frame size is kept the same, and the end is
% padded with zeros if needed to fill the last frame.
%   Filter is butter(4, fc, 'low'), where fc = 0.9/N, run forward and
% backward to yield 8th-order filtering with zero phase shift.
% The <src> and <dest> file formats are assumed to be Neuralynx .ncs unless
% the file name extension is .mat, in which case they are assumed to be the
% format written by dg_Nlx2Mat.  Note that there are no unix/Mac Nlx file
% creator functions, so on those platforms the output must ALWAYS be .mat.
%   Special handling for empty <src> files is to do nothing, not even
% create an empty outputfile, but issue a warning.
%OPTIONS
% 'filter', factor - set filter fc = <factor>/N, where <factor> must be
% less than or equal to 0.9.
% 'int16' - writes output .mat file with sample data in int16
%   representation.  Note that this will result in a file of all zeros if
%   typical physiological data are read from a Nlx file without specifying
%   'nounits'.
% 'no_subTScheck' - skips the check for error 'dg_downsample:subTScheck'.
% 'nounits' - leaves sample data in unscaled (i.e. raw) A/D units when
%   reading Nlx format.  Has no effect on .mat format.
% 'overwrite' - raises a warning instead of an error and overwrites the
%   existing <dest> file if there is one.
% 'reconcile', ts - removes frames from raw file as needed to match
%   the timestamps in <ts>.  If <src> is a .mat file containing samples in
%   vector form, then the frame size is assumed to be 512, and it is an
%   error if the number of samples is not divisible by 512.  Timestamps are
%   assumed to be integer-valued and so must match exactly.  <ts> would
%   typically be the <canonicalstamps> value returned by
%   dg_reconcileFrames.
% 'skip' - raises a warning instead of an error and skips the
%   existing <dest> file if there is one.
% 'verbose' - this is an option which optionally can be specified as an
%   option in the argument list when calling the function from the command
%   line or from a file or in any other manner.
%NOTES
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
%   It's a (slightly) comforting thought that since N is an integer, and
% the final output file has the same frame size as the input file, all of
% the timestamps in the output file will be members of the timestamps of
% the input file.  In the simple case, they are just every Nth timestamp
% starting with the first (so, 1 + N*(0:numsubframes-1)).  In the case of
% the 'reconcile' option, there will be hiccups, so the offsets will
% change, but the general pattern of "usually every Nth timestamp" remains,
% and it is strictly true that they are RECONCILED timestamps numbers ( 1 +
% N*(0:numsubframes-1)).
%   I have belatedly (19-Mar-2015) realized that it would have avoided a
% lot of headaches if I had constrained the non-overlappd part of the batch
% size to be an integral multiple of <N>, so that the first subsample in
% the non-overlapping batch would always be the first sample in the frame,
% and therefore the first subsample in every batch except for the first
% would always be the first sample of the second frame read in the batch.

%$Rev: 305 $
%$Date: 2023-09-08 14:19:01 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

fac = 0.9;
framesize = NaN;
int16flag = false;
% a variable into which to load .mat files, to keep nested functions happy:
m = [];
no_subTScheckflag = false;
nounitsflag = false;
overwriteflag = false;
rawsamples = [];
rawts = [];
recTS = []; % reconciled timestamps
skipflag = false;
verboseflag = false;

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'filter'
            argnum = argnum + 1;
            fac = varargin{argnum};
        case 'int16'
            int16flag = true;
        case 'no_subTScheck'
            no_subTScheckflag = true;
        case 'nounits'
            nounitsflag = true;
        case 'overwrite'
            overwriteflag = true;
        case 'reconcile'
            argnum = argnum + 1;
            recTS = varargin{argnum};
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

if fac > 0.9
    error('dg_downsample:factor', ...
        'The ''filter'' factor cannot be > 0.9');
end
if ~exist(src, 'file')
    error('dg_downsample:nofile', ...
        'The input file %s does not exist', src);
end
if exist(dest, 'file')
    if overwriteflag
        warning('dg_downsample:fileexists', ...
            'The output file %s already exists; overwriting', dest);
    elseif skipflag
        warning('dg_downsample:fileexists2', ...
            'The output file %s already exists; skipping', dest);
        return
    else
        error('dg_downsample:fileexists', ...
            'The output file %s already exists', dest);
    end
end
if fix(N) ~= N
    error('dg_downsample:N', ...
        'The downsampling factor must be an integer.');
end
if ~isempty(recTS) && no_subTScheckflag
    warning('dg_downsample:no_subTS', ...
        ['You have specified timestamp reconciliation, but have also ' ...
        'specified no_subTScheck, which disables the verification of ' ...
        'the final timestamps.']);
end

if verboseflag
    str = (sprintf('dg_downsample(%s, %s, %d', src, dest, N));
    for k1 = 1:length(varargin)
        str = [str ', ' dg_thing2str(varargin{k1})]; %#ok<*AGROW>
    end
    str = [str ')'];
    disp(str);
end

% Set up constant things:
[srcdir, srcname, srcext] = fileparts(src); %#ok<*ASGLU>
[destdir, destname, destext] = fileparts(dest);
if ~ispc && ~isequal(lower(destext), '.mat')
    error('dg_downsample:badoutput', ...
        'Neuralynx format files cannot be created on non-Windows platforms.');
end
fc = fac/N;
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
    disp(sprintf('Reading file %s', src)); %#ok<*DSPS>
end
% Get lowpass filtered data ready for downsampling in <rawsamples>:
if isequal(lower(srcext), '.mat')
    readmatfile;
else
    readnlxfile;
end
% 'readnlxfile' in batch mode does not leave any <rawsamples> behind, but
% there should always be <subsamples>:
if exist('rawsamples', 'var') && isempty(rawsamples) || ...
        ~exist('subsamples', 'var') || isempty(subsamples)
    warning('dg_downsample:empty', ...
        'The input file is empty; no output file was produced.');
    return
end

% pad if needed and reshape into frames:
extrasamples = rem(numel(subsamples), framesize);
padding = [];
if extrasamples ~= 0
    warning( 'dg_downsample:padding', ...
        'padded "%s" with zeros to fill last frame', destname );
    padding = zeros(1, framesize-extrasamples);
end
subsamples = reshape([subsamples padding], framesize, []);
allsubTS = reshape([allsubTS padding], framesize, []);
subTS = allsubTS(1,:);
if ~isequal(subTS, targetTS)
    error('dg_downsample:TS', ...
        'Timestamp discrepancy indicating program error.');
end

% Create output file:
if isequal(lower(destext), '.mat')
    dg_Nlx2Mat_Timestamps = subTS; %#ok<NASGU>
    dg_Nlx2Mat_Samples = subsamples;
    dg_Nlx2Mat_SamplesUnits = units; %#ok<NASGU>
    if int16flag
        dg_Nlx2Mat_Samples = int16(dg_Nlx2Mat_Samples); %#ok<NASGU>
    end
    save(dest, ...
        'dg_Nlx2Mat_Timestamps', 'dg_Nlx2Mat_Samples', ...
        'dg_Nlx2Mat_SamplesUnits', '-v7.3' );
else
    if isequal(lower(srcext), '.mat')
        maxsamp = max(subsamples(:));
        minsamp = min(subsamples(:));
        if maxsamp >= -minsamp * 2047/2048
            subsamples = round(2047 * subsamples / maxsamp);
        else
            subsamples = round(2048 * subsamples / -minsamp);
        end
        header = {'######## Neuralynx Data File Header'
            sprintf('## Converted by dg_downsample from: %s ', src)
            sprintf( '## Time Opened (dd-mmm-yyyy HH:MM:SS): %s ', ...
            datestr(now, 0) )};
    end
    if verboseflag
        disp(sprintf('Writing file %s', dest));
    end
    dg_writeCSC(dest, subTS, subsamples, header);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function local_filtfilt
        % Lowpass filter <rawsamples> in place.
        % Keeps shape of <rawsamples>.
        samplesize = size(rawsamples);
        if verboseflag
            disp(sprintf('Filtering first time'));
        end
        rawsamples = filter(h2, rawsamples(end:-1:1));
        if verboseflag
            disp(sprintf('Filtering second time'));
        end
        rawsamples = reshape(filter(h2, rawsamples(end:-1:1)), samplesize);
    end

    function readmatfile
        % Life is simple because a .mat file had to fit in memory in one
        % batch. This also lowpass filters it, leaving the result in
        % <rawsamples>.  Times here are in seconds.
        m = load('-MAT', src);
        rawts = m.dg_Nlx2Mat_Timestamps;
        rawsamples = m.dg_Nlx2Mat_Samples;
        % Guarantee that rawsamples is in samples x frames format:
        if size(rawsamples,2) ~= length(rawts)
            rawsamples = reshape(rawsamples, [], length(rawts));
        end
        framesize = size(rawsamples, 1);
        if isfield(m, 'dg_Nlx2Mat_SamplesUnits')
            units = m.dg_Nlx2Mat_SamplesUnits;
        else
            units = 'arbs';
        end
        clear m;
        % Compute <targetTS>, the intended frame timestamps in the final
        % downsampled output file:
        if isempty(recTS)
            recTS = rawts;
        end
        numsubframes = ceil(length(recTS)/N);
        targetTS = recTS(1 + N*(0:numsubframes-1));
        if ~isempty(recTS) && ~isequal(rawts, recTS)
            if length(recTS) > length(rawts)
                error('dg_downsample:badrecTS', ...
                    'The timestamp list specified with ''reconcile'' is longer than that from the file %s', ...
                    src );
            end
            if size(rawsamples,2) ~= 512
                rawsamples = reshape(rawsamples, 512, []);
            end
            [c, frames2del] = setdiff(rawts, recTS);
            rawsamples(:, frames2del) = [];
        end
        local_filtfilt;
        downsampleSingleBatch;
    end

    function readnlxfile
        % Life is miserable, because raw Nlx files can be arbitrarily
        % large, and so the processing may need to be batched.  That
        % implies that if recTS is empty, we need to read the timestamps
        % from the input file in order to to calculate <targetTS>
        % (reconciled timestamps intended for final output file).  Then we
        % need to create the temporary files containing the batched
        % downsampled data, allowing enough overlap so that the startup
        % transients can be deleted on both ends, and finally stitch them
        % back together and leave them in <rawsamples>.
        % Times here are in microseconds.
        ts2 = dg_readCSC(src);
        if isempty(ts2)
            return
        end
        if isempty(recTS)
            recs2read = 1:numel(ts2);
            recTS = ts2;
        else
            if any(~ismember(recTS, ts2))
                error('dg_downsample:recTS', ...
                    '%s is missing %d records specified by ''reconcile'' option', ...
                    src, sum(~ismember(recTS, ts2)) );
            end
            recs2read = find(ismember(ts2, recTS));
        end
        numrecs2read = length(recs2read);
        clear ts2;
        numsubframes = ceil(length(recTS)/N);
        targetTS = recTS(1 + N*(0:numsubframes-1));
        if ispc
            userview = memory;
            bytesinblock = userview.MaxPossibleArrayBytes;
        else
            % There is no 'memory' function in unix/Mac, so we simply
            % assume that there is 1 GB available:
            bytesinblock = 2^30;
        end
        [TS, Samples, header] = dg_readCSC(src, 'header');
        % We assume that all Neuralynx files are framesize = 512:
        framesize = 512;
        recsize = 8 * framesize + 32;
        % Numbers in recs2read denote record numbers in the source file,
        % not frame numbers in the reconciled data. recs2read(startrecidx)
        % is the first record in the current batch, and always refers to
        % the record in the src file (but using the Matlab indexing
        % convention of "first index is 1").  Likewise for
        % recs2read(endrecidx), which is the last in the current batch.
        startrecidx = 1;
        if bytesinblock < 24
            error('dg_downsample:outofmem', ...
                'There are only %d bytes available; pack memory', ...
                bytesinblock);
        end
        maxbytes2use = bytesinblock/10;
        % calculate endrecidx and s:
        if maxbytes2use < numrecs2read * recsize
            batchsize = fix(maxbytes2use/recsize);
            endrecidx = batchsize; % record in the src file
            % There will actually be slightly more than
            % ceil(numrecs2read/batchsize) batches because two frames (records)
            % are thrown away at each batch boundary. Analysis of extreme
            % cases:
            %   batchsize = 2 makes the overlap impossible, so therefore we
            %       must always have batchsize > 2.
            %   maxbytes2use = numrecs2read * recsize; this triggers the single
            %   batch method maxbytes2use = numrecs2read * recsize - 1
            %       batchsize = 3 => numbatches = numrecs2read - 2 recsize = 1
            %       =>  batchsize = maxbytes2use =>
            %           numbatches = 2 + ceil(-(numrecs2read-4)/(numrecs2read-3)), or
            %           numbatches = 2 - floor((numrecs2read-4)/(numrecs2read-3))
            %           numbatches = 2 - 0 for numrecs2read > 4 numbatches = 2 -
            %           0 for numrecs2read = 4 impossible for numrecs2read = 3, i.e.
            %           for maxbytes2use = 2
            % That last case determines the minimum bytesinblock required.
            numbatches = 2 + ceil((numrecs2read - 2*(batchsize-1))/(batchsize-2));
        else
            endrecidx = numrecs2read;
            numbatches = 1;
        end
        if verboseflag
            disp(sprintf('Reading in ~%d batches', numbatches));
        end
        done = false;
        batchnum = 1;
        numsubsamples = 0;
        startidx = 1;
        while ~done
            if verboseflag
                disp(sprintf( ...
                    '%d/~%d: Reading records %.0f-%.0f', ...
                    batchnum, numbatches, ...
                    recs2read(startrecidx), recs2read(endrecidx) ));
            end
            % Using extraction mode 3 would work against the purpose of
            % batching, which is to limit the amount of memory being used.
            % Therefore we must endure the headache of searching each batch
            % of recs2read for any jumps in the record sequence so that we
            % can read each sequence of consecutive records using one call
            % to dg_readCSC in extraction mode 2.
            numrecsinbatch = endrecidx - startrecidx + 1;
            rawts = zeros(1, numrecsinbatch);
            rawsamples = zeros(framesize, numrecsinbatch);
            % subbatches works exactly like startrecidx, i.e. it contains
            % indices into <recs2read>.
            numsubbatches = sum( ...
                diff(recs2read(startrecidx:endrecidx)) ~= 1 ) + 1;
            if numsubbatches == 1
                subbatches = [startrecidx endrecidx];
            else
                subbatches = zeros(numsubbatches, 2);
                % <brkidx> points to the last element of recs2read in each
                % subbatch. 
                brkidx = find( diff(recs2read(startrecidx:endrecidx)) ...
                    ~= 1 ) + startrecidx - 1;
                subbatches(1,:) = [startrecidx brkidx(1)];
                for subbatchnum = 2 : (numsubbatches - 1)
                    subbatches(subbatchnum,:) = ...
                        [brkidx(subbatchnum-1)+1, brkidx(subbatchnum)];
                end
                subbatches(end,:) = [brkidx(end)+1, endrecidx];
            end
            for subbatchnum = 1 : size(subbatches,1)
                relrecidx = subbatches(subbatchnum,:) - startrecidx + 1;
                % Note that <rawts> is in integral microseconds.
                [rawts(relrecidx(1):relrecidx(2)), ...
                    rawsamples(:, relrecidx(1):relrecidx(2))] = ...
                    dg_readCSC(src, 'mode', 2, ...
                    recs2read(subbatches(subbatchnum,:)));
            end
            framedurs = rawts(2:end)-rawts(1:end-1);
            if any(framedurs<1)
                warning('dg_downsample:framedurs', ...
                    'There is/are disordered timestamps at record(s) %s in %s', ...
                    dg_thing2str(find(framedurs<1)+1), src);
            end
            % Set units:
            if isequal(lower(destext), '.mat')
                if nounitsflag
                    units = 'AD';
                else
                    % Convert from AD units to volts if possible
                    if exist('ADBitVolts', 'var')
                        if ~isempty(ADBitVolts)
                            rawsamples = ADBitVolts ...
                                * rawsamples;
                            units = 'V';
                        end
                    else
                        units = 'AD';
                        for k3 = 1:length(header)
                            if regexp(header{k3}, '^\s*-ADBitVolts\s+')
                                ADBitVoltstr = regexprep(header{k3}, ...
                                    '^\s*-ADBitVolts\s+', '');
                                ADBitVolts = str2double(ADBitVoltstr);
                                if isempty(ADBitVolts)
                                    warning('dg_downsample:badADBitVolts', ...
                                        'Could not convert number from:\n%s', ...
                                        header{k3} );
                                else
                                    rawsamples = ADBitVolts ...
                                        * rawsamples;
                                    units = 'V';
                                end
                                break
                            end
                        end
                    end
                end
            end
            local_filtfilt;
            % We now have lowpass filtered data from startrecidx through
            % endrecidx in rawsamples, with bad frames deleted.
            if batchnum == 1
                % compute sample period <sampleT>:
                if size(rawsamples,1) ~= framesize
                    error('dg_downsample:framesize', ...
                        '%s has wrong frame size (%d)', ...
                        src, size(rawsamples,1));
                end
                minfd = min(framedurs(framedurs>0));
                framedurs = framedurs(framedurs < 2*minfd);
                medianfd = median(framedurs);
                sampleT = medianfd/framesize;
            else
                % delete first frame in non-first batches:
                rawsamples(:,1) = [];
                rawts(1) = [];
            end
            if endrecidx ~= numrecs2read
                % delete last frame in non-last batches:
                rawsamples(:,end) = [];
                rawts(end) = [];
            end
            if numbatches == 1
                downsampleSingleBatch;
                done = true;
            else
                % <allrawTS> contains timestamps of all samples in current
                % batch, in the same shape array as <rawsamples>.
                allrawTS = repmat(rawts, framesize, 1) + ...
                    repmat(round((0 : framesize - 1)' * sampleT), ...
                    1, length(rawts));
                if verboseflag
                    disp(sprintf( ...
                        'Downsampling %dX, records %.0f-%.0f', ...
                        N, startrecidx, endrecidx));
                end
                % do the downsampling. tempsamples is the final result of
                % the current batch, to save to tempfile with matching
                % timestamps tempallTS.
                tempsamples = reshape(rawsamples(startidx:N:end), 1, []);
                tempallTS = reshape(allrawTS(startidx:N:end), 1, []);
                if int16flag
                    tempsamples = int16(tempsamples);
                end
                tempfilename = sprintf('%s%07d.mat', tempbase, startrecidx);
                if verboseflag
                    disp(sprintf('Saving file %s', tempfilename));
                end
                % save temp file
                save(tempfilename, 'tempsamples', 'tempallTS');
                numsubsamples = numsubsamples + numel(tempsamples);
                % Updating startidx is something very annoying: it is the
                % index of the first sample in the reconciled filtered
                % batch data that is to be used as a subsample.  This is
                % non-trivial because it can only be calculated based on
                % the current startidx and the actual number of samples in
                % the current batch.  If there are S subsamples in the
                % batch, then if there were no overlap between batches, the
                % next value of startidx would be given by
                %   startidx = mod(startidx + S*N -1, framesize) + 1
                % (The +/- 1 is because Matlab foolishly considers the
                % first element's index to be 1 and the last's to be
                % <framesize>, whereas 'mod' returns values from 0 to
                % <framesize>-1.) Mercifully, by the time we get here,
                % there actually is no overlap between the batches because
                % the overlap frames have already been deleted.
                S = length(tempsamples);
                startidx = mod(startidx + S*N -1, framesize) + 1;
                if endrecidx == numrecs2read
                    done = true;
                else
                    % update indices for next batch:
                    startrecidx = endrecidx - 1;  % repeat last two frames
                    endrecidx = min(startrecidx + batchsize - 1, numrecs2read);
                    batchnum = batchnum + 1;
                end
            end % if numbatches == 1; else...
        end
        if numbatches > 1
            % Reconstruct entire signal with timestamps into row vectors
            % <subsamples> and <allsubTS> from batches saved in tempfiles,
            % deleting each file after reading:
            tempfiles = dir([tempbase '*.mat']);
            clear rawsamples;
            subsamples = NaN(1, numsubsamples);
            allsubTS = NaN(1, numsubsamples);
            lastsampidx = 0;
            for k4 = 1:length(tempfiles)
                % This reloads previously incarnated variables, so it does
                % not have to be assigned to a variable:
                load(tempfiles(k4).name);
                % samplerange is an index into subsmaples and allsubTS.
                samplerange = lastsampidx + (1:length(tempsamples));
                % tempsamples, tempallTS are loaded from tempfiles(k4).name
                % and are already trimmed properly to go into subsamples,
                % allsubTS:
                subsamples(samplerange) = tempsamples;
                allsubTS(samplerange) = tempallTS;
                if ~no_subTScheckflag
                    subTScheck = allsubTS(1:framesize:samplerange(end));
                    if ~all(subTScheck == targetTS(1:length(subTScheck)))
                        errmsg = sprintf( ...
                            'Timestamp discrepancy in temp file %s\n', ...
                            tempfiles(k4).name);
                        erridx = find(subTScheck ~= ...
                            targetTS(1:length(subTScheck)), 1);
                        errmsg = [errmsg sprintf( ...
                            'First mismatch: subTScheck=%g, targetTS=%g\n', ...
                            subTScheck(erridx), targetTS(erridx) )];
                        errmsg = [errmsg sprintf( ...
                            'Number of mismatches: %d\n', sum( subTScheck ...
                            ~= targetTS(1:length(subTScheck)) )) ];
                        erridx = find(subTScheck ~= ...
                            targetTS(1:length(subTScheck)), 1, 'last');
                        errmsg = [errmsg sprintf( ...
                            'Last mismatch: subTScheck=%g, targetTS=%g', ...
                            subTScheck(erridx), targetTS(erridx) )];
                        error('dg_downsample:subTScheck', '%s', errmsg);
                    end
                end
                lastsampidx = lastsampidx + length(tempsamples);
            end
            for k4 = 1:length(tempfiles)
                delete(tempfiles(k4).name);
            end
        end
    end

    function downsampleSingleBatch
        % Entire file is already in memory; downsample it.
        % compute sample period <sampleT>:
        framedurs = rawts(2:end)-rawts(1:end-1);
        if any(framedurs<0)
            badframes = find(framedurs<0);
            badTSstr = '';
            for k2 = 1:length(badframes)
                badTSstr = sprintf( '%s%.0f: %.0f, %.0f: %.0f\n', badTSstr, ...
                    badframes(k2), rawts(badframes(k2)), ...
                    badframes(k2)+1, rawts(badframes(k2)+1) );
            end
            error('dg_downsample:badTS', ...
                'Time goes backwards between (framenumber: timestamp)\n%s', ...
                badTSstr);
        end
        minfd = min(framedurs);
        framedurs = framedurs(framedurs < 2*minfd);
        medianfd = median(framedurs);
        sampleT = medianfd/framesize;
        % downsample into row vector:
        allrawTS = repmat(rawts, framesize, 1) + ...
            repmat((0 : framesize - 1)' * sampleT, 1, length(rawts));
        if verboseflag
            disp(sprintf('Downsampling %dX', N));
        end
        subsamples = reshape(rawsamples(1:N:end), 1, []);
        clear rawsamples;
        allsubTS = reshape(allrawTS(1:N:end), 1, []);
    end

end
