function spikewaves = dg_lookupSpikeWaveforms(varargin)
%spikewaves = dg_lookupSpikeWaveforms(fn1)
%spikewaves = dg_lookupSpikeWaveforms(fn1, fn2)
% Reads timestamps and cluster IDs from the 'Multiple Cut Cluster (*.MAT)'
% file <fn1>, finds the same spikes in the raw Neuralynx tetrode file
% <fn2>, and appends the spike waveforms to the data from <fn1>.  <fn1> is
% assumed to match the regexp '-T\d\d+' (e.g. "acq10-T01.MAT" or
% "bertha-T567.MAT") and to be located in the animal directory in the Yasuo
% standard data directory tree format.  The part of the filename after the
% "-T" is the tetrode number and the part before is the session ID.
%   Timestamps and cluster IDs can also be submitted directly as a numeric
% array value for <fn1>, in which case no cluster file is read.  <fn2> can
% be omitted or empty, in which case a "best effort" is made to find the
% the Neuralynx file that matches <fn1>; in this case <fn1> must be a
% file pathname, not a numeric array.
%   If more than one matching file is found, the last one found is used.
% If the Neuralynx file contains an "-ADBitVolts" line with 4 elements in
% the header, then the waveform data are converted from A/D units into
% volts.  If one wire consistently contains peak amplitudes that are less
% than one-tenth the peak amplitudes on the good wires, then the warning
% 'dg_lookupSpikeWaveforms:badwire' is raised, in which case you should
% consider using the 'badwire' option.
%OPTIONS
% 'badwire', badwirelist - ignores wires (numbered 1 through 4) that appear
%   in the array <badwirelist>.
% 'ktfilenames' - <fn1> matches the regexp '_\d\d.*' (e.g. "thing_01.MAT"
%   or "acq2_12.MAT"), and resides in the 'onlytrialcut\ccw' or
%   'onlytrialcut' subfolder of the session folder.  The tetrode number is
%   exactly the first two characters after the underscore ('_'), and the
%   rest of the filename is ignored.
% 'plexonnames' - <fn1> matches the regexp '_\d\d.*' as in 'ktfilenames',
%   but the numerical part is ONE LESS than the tetrode number, and the
%   <fn1> is assumed to be located in the session directory, not the animal
%   directory.  E.g., all of the following pathnames will be interpreted as
%   referring to tetrode #3 in session 'acq27' of animal 'sumrat':
%       C:\sumrat\acq27\acq27_002.MAT
%       C:\sumrat\acq27\acq32_002.MAT
%       C:\RData\sumrat\acq27\john-doe_002.MAT
%       my_cooked_data_dir\sumrat\acq27\_002.MAT
% 'verbose'

%$Rev: 160 $
%$Date: 2012-10-31 14:43:11 -0400 (Wed, 31 Oct 2012) $
%$Author: dgibson $

argnum = 1;
args2delete = [];
badwirelist = [];
ktfflag = false;
plexonflag = false;
verboseflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'badwire'
            args2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            badwirelist = varargin{argnum};
            args2delete(end+1) = argnum;
        case 'ktfilenames'
            ktfflag = true;
            args2delete(end+1) = argnum;
        case 'plexonnames'
            plexonflag = true;
            args2delete(end+1) = argnum;
        case 'verbose'
            verboseflag = true;
            args2delete(end+1) = argnum;
    end
    argnum = argnum + 1;
end
varargin(args2delete) = [];

fn1 = varargin{1};
if length(varargin) >=2
    fn2 = varargin{2};
end

if length(varargin) < 2 || isempty(fn2)
    if isa(fn1, 'numeric')
        error('dg_lookupSpikeWaveforms:badfn1', ...
            '<fn2> must be specified when <fn1> is numeric data.' );
    end
    fn2 = '';
end

% <clustfileregexp> matches the part after the sessionID, not including ext
if ktfflag || plexonflag
    clustfileregexp = '_\d\d.*';
else
    clustfileregexp = '-T\d\d+';
end

if isa(fn1, 'numeric')
    spikewaves = fn1;
else
    if plexonflag
        [p, name, ext] = fileparts(fn1);
        [animaldir, sessionID] = fileparts(p);
    else
        [animaldir, name, ext] = fileparts(fn1);
        if isempty(animaldir)
            animaldir = pwd;
        end
    end
    if ~isequal(upper(ext), '.MAT')
        error('dg_lookupSpikeWaveforms:badfn1ext', ...
            '<fn1> must be a ''Multiple Cut Cluster (*.MAT)'' file.' );
    end
    
    if isempty(fn2)
        delim1 = regexpi(name, clustfileregexp);
        if isempty(delim1)
            error('dg_lookupSpikeWaveforms:badfn1name', ...
                '<fn1> does not match pattern for cluster file.' );
        end
        if ktfflag
            clustfiledir = fileparts(fn1);
            [p, dirname] = fileparts(clustfiledir);
            % clustfiledir is either 'onlytrialcut\ccw' or 'onlytrialcut'
            % subdir of sessiondir:
            if isequal(dirname, 'ccw')
                p = fileparts(p);
            end
            sessiondir = p;
            trodestr = name(delim1+1:delim1+2);
            trodenum = str2double(trodestr);
        elseif plexonflag
            sessiondir = fullfile(animaldir, sessionID);
            trodestr = name(delim1+1:end);
            trodenum = str2double(trodestr) + 1;
        else
            % Attempt to parse name portion of fn1.  Examples:
            %   ACQ06-T01.MAT (same interpretation as ACQ06.T01 T-file)
            sessionID = name(1 : delim1 - 1);
            sessiondir = fullfile(animaldir, sessionID);
            trodestr = name(delim1 + 2 : end);
            % chop any non-digits from end:
            delim2 = regexp(trodestr, '\D');
            if ~isempty(delim2)
                trodestr = trodestr(1:delim2(1));
            end
            trodenum = str2double(trodestr);
        end
        
        % Attempt to find matching Neuralynx file.
        nlxname{1} = sprintf('TT%d', trodenum);
        nlxname{2} = sprintf('T%02d', trodenum);
        nlxext{1} = '.ntt';
        nlxext{2} = '.dat';
        for n = 1:2
            for x = 1:2
                myfn = fullfile(sessiondir, [nlxname{n} nlxext{x}]);
                if exist(myfn, 'file')
                    fn2 = myfn;
                    if verboseflag
                        disp(sprintf('Found %s', fn2)); %#ok<*DSPS>
                    end
                else
                    if verboseflag
                        disp(sprintf('%s does not exist.', myfn));
                    end
                end
            end
        end
        if isempty(fn2)
            error('dg_lookupSpikeWaveforms:nofn2', ...
                'Could not find Neuralynx file.' );
        end
    end
    
    % Load the cluster data from the .MAT file:
    S = load(fn1);
    fields = fieldnames(S);
    spikewaves = S.(fields{1});
    clear S;
    if size(spikewaves,2) > 2
        warning('dg_lookupSpikeWaveforms:gotwaves', ...
            'The file %s already contains waveform data; ignoring.', ...
            fn1);
        spikewaves(:,3:end) = [];
    end
end

% Read raw timestamps and waveforms; convert timestamps to seconds.
if exist(fn2, 'file')
    [TS, points, header] = dg_readSpike(fn2);
else
    error('dg_lookupSpikeWaveforms:nofile', ...
        'No such file: %s', fn2);
end
TS = TS * 1e-6;

% Blown wire test
goodwires = setdiff(1:4, badwirelist);
themaxes = max(abs(points), [], 1);
themedians = squeeze(median(themaxes, 3));
if any(themedians(goodwires) < median(themedians(goodwires))/10)
    badwire = goodwires( ...
        themedians(goodwires) < median(themedians(goodwires))/10 );
    warning('dg_lookupSpikeWaveforms:badwire', ...
        'Probable bad wire #%s in %s', mat2str(badwire), fn2);
end

for k = 1:length(header)
    if regexp(header{k}, '^\s*-ADBitVolts\s+')
        ADBitVolts = sscanf(regexprep( header{k}, ...
            '^\s*-ADBitVolts\s+', '' ), '%g');
        if isempty(ADBitVolts)
            warning('dg_lookupSpikeWaveforms:badADBitVolts', ...
                'Could not convert number from:\n%s', ...
                header{k} );
        elseif length(ADBitVolts) ~= 4
            warning('dg_lookupSpikeWaveforms:badADBitVolts2', ...
                'There are not 4 gains in:\n%s', ...
                header{k} );
        else
            if any(points(:) == -2048 | points(:) == 2047)
                warning('dg_lookupSpikeWaveforms:clipping', ...
                    'There are clipped samples.' );
            end
            for w = 1:4
                points(:,w,:) = points(:,w,:)*ADBitVolts(w);
            end
        end
        break
    end
end

% Find each spike in spikewaves, and add its waveforms to spikewaves.
spikewaves = [ spikewaves zeros( size(spikewaves,1), ...
    size(points,1) * size(points,2) ) ];
for spidx = 1:size(spikewaves,1)
    ptsidx = find(abs(TS - spikewaves(spidx,1)) < 1.01e-4);
    if isempty(ptsidx)
        warning('dg_lookupSpikeWaveforms:spikemissing', ...
            'Could not find spike at %.6f s', spikewaves(spidx,1));
        continue
    end
    if length(ptsidx) > 1
        warning('dg_lookupSpikeWaveforms:spikedoubled', ...
            'More than one spike at %.6f s; using first', ...
            spikewaves(spidx,1));
        ptsidx = ptsidx(1);
    end
    spikewaves(spidx, 3:end) = reshape(points(:, :, ptsidx), 1, []);
end

