function lfp_fishAnalysis(tackle, sessiondir, varargin)
%lfp_fishAnalysis(tackle, sessiondir)
%LFP ANALYSIS
% The LFP signal is cleaned and smoothed as follows.  If there is a file
% named CSC_LOCALAVGREF1.MAT in <sessiondir>, then it is assumed to also
% exist in the fragment subdirectories, and it is subtracted from each
% other CSC channel before analyzing. The harmonics of 60 Hz are removed
% using lfp_rmlnspec, the baseline DC level is removed, and the resulting
% trace is smoothed using a Hanning window of width <smoothtime> (5 ms
% default). Extrema in the smoothed wave are then found using
% dg_findFlattops, and are marked if they are significantly different from
% the preceding OR following extremum, AND are significantly different from
% baseline (plotted as gray horizontal line). It is assumed that the
% variance of the LFP trace is much greater than the variance in the
% estimated baseline value, and so significance is determined entirely by
% the confidence limits on the LFP trace.
%MULTIUNIT ACTIVITY ANALYSIS
% MUA is extracted from each file whose name matches 'tt*.ntt' by simply
%   using all the records in the file as a multiunit spike train.  The
%   spike train is then converted to a spike density function using the
%   'multispike' option to lfp_spike2wave and smoothing width specified by
%   <smoothtime>.  The rest of the analysis is the same as for the cleaned
%   and smoothed LFP.
%SINGLE UNIT ACTIVITY ANALYSIS
% Same as MUA analysis, except that all '*.mat' files are assumed to be
%   Multiple Cut Cluster files, and each cluster is analyzed separately.
%INPUTS
%   tackle:  as for lfp_fishStory, except some fields may not be used.
%   sessiondir: must contain 'baseline' and 'block1', 'block2', etc.
%       fragment subdirectories.  Each fragment must contain the same set
%       of CSC files as the root session directory, and CSC files must have
%       names that begin with 'csc_'.
%OUTPUT
% No return values.  Optionally produces plots and/or files of various
% types (see OPTIONS).
%OPTIONS
%   'baseline', blspec - override the baseline block result with a baseline
%       computed according to <blspec>, which may be:
%           'block': (default) compute from baseline block, with SD assumed
%           to be negligibly small.
%           'none': no baseline subtraction of any kind.
%           {evtIDs win}: <blspec> is a two-element cell array, where the
%               first element, <evtIDs>, is a list of numerical event IDs
%               that is used like lfp_Alignmentref to align the baseline
%               window, and the second element, <win>, is a two-element
%               time window specified in seconds relative to the first
%               instance of <evtIDs> in the trial.  In this case, SD is
%               computed over all samples in <win>, and significance is
%               judged based on sqrt(baseSD^2 + trialSD^2).
%       Applies both to LFP and spike analyses.  WARNING: there is no
%       protection against poorly chosen baseline windows that contain
%       stimulation artifact (which can appear before the stimulation
%       marker due to backwards filtering) or EPs.
%   'copycol', [from to] - superimposes a copy of the mean trace from
%       column number <from> on the plot in column number <to>.
%   'endTime', endTime - of analysis period relative to stim.  Default is
%       1.995 seconds.  Note that this can cause trouble if the analysis
%       time window is long enough to run into other constraints (e.g.
%       recorded segment boundaries), and will probably cause a crash in
%       lfp_selectByDuration.
%   'filenamestr', filenamestr - <filenamestr> is an operating system
%       dependent string that may contain wildcards (e.g. '*') to use when
%       searching for files to read.  Default value is 'tt*.ntt' for MUA,
%       and '*.mat' for SUA; default behavior for LFPs is to use the entire
%       set of files determined by tackle.localavgrefs as described in
%       lfp_fishStory and dg_localAvgRef, and is completely overridden by
%       specifying 'filenamestr' except that the 'csc_localavgref' files
%       are still excluded.
%   'LFPonly' - skips all unit activity analysis.
%   'MUAonly' - skips LFP analysis and single unit activity analysis.
%   'noplot' - suppresses creation of Matlab Figure windows.  Also see
%       'saveresult' option.
%   'numSEMs', n - sets confidence limits for plots and significance
%       judgments to be <n> times the Standard Error of the Mean.  Default
%       value is n=2.  <n> can have a fractional value.
%   'read' - reads the previously saved non-smoothed waveforms and applies
%       the specified smoothing.
%   'readsmooth' - reads the previously saved smoothed waveforms.
%   'save' - saves the calculated waveforms immediately before doing the
%       smoothing, with '_clean' or '_MUA' appended to filename.  Does not
%       save unsmoothed single unit wave data. (The '_MUA' file is actually
%       saved for every MUA analysis regardless of whether 'save' is
%       specified or not, and it is of type 'Single Cut Cluster (*.MAT)'
%       rather than CSC.)
%   'saveresult', resultdir - saves the final results, i.e. everything that
%       would normally be plotted in the Figure windows, as a pair of .mat
%       files in <resultdir>, where each .mat file's name is of the form
%       <LFP|MUA><SessionID>.mat, where <SessionID> is the name of the
%       parent directory of <sessiondir>.  These output files can be
%       plotted using lfp_plotFishAnalysisFile.
%   'savesmooth' - saves the final smoothed waveforms, with '_clean<N>'
%       or '_MUA<N>' appended to filename, where <N> is the number of
%       milliseconds of smoothing.
%   'smoothtime', smoothtime - width of smoothing window in milliseconds;
%       default = 5.
%   'spikedir', spikedir - looks for raw tetrode spike files in <spikedir>.
%   'startTime', startTime - of analysis period relative to stim, in
%       seconds.  Default is tackle.artifactTimeout.  Note that this can
%       cause trouble if the analysis time window is long enough to run
%       into other constraints (e.g. recorded segment boundaries), and will
%       probably cause a crash in lfp_selectByDuration.
%   'SUAregexp', SUAregexp - uses the regular expression <SUAregexp> to
%       find Multiple Cut Cluster files in conjunction with <filenamestr>
%       (see 'filenamestr' option); that is, a filename must match both.
%   'SUAonly' - skips LFP analysis and multiunit activity analysis.
%   'sum' - on the plot at the right hand end of each row, superimposes a
%       plain blue trace plotting the sum of the mean signals from the
%       other columns.
%   'sumcol', sumcol - same as 'sum' except that it is column number
%       <sumcol> that is left out of the sum and on which the summed
%       waveform is superimposed.
%   'verbose' - modifies output to command window in such a manner as to
%       produce an effect commensurate with the name of the option.

%$Rev: 289 $
%$Date: 2012-12-10 17:42:50 -0500 (Mon, 10 Dec 2012) $
%$Author: dgibson $

global lfp_SpikeNames lfp_FileNames lfp_ActiveFilenums

argnum = 1;
filenamestr = '';
opts.blspec = 'block';
opts.copycol = [0 0];
opts.endTime = 1.995; % of analysis period re: stim to exclude artifact
opts.LFPonlyflag = false;
opts.MUAonlyflag = false;
opts.numSEMs = 2;
opts.readflag = false;
opts.readsmoothflag = false;
opts.saveflag = false;
opts.savesmoothflag = false;
opts.smoothtime = 5; % half-width of Hanning window in milliseconds
opts.spikedir = sessiondir;
% of analysis period re: stim to exclude artifact:
opts.startTime = tackle.artifactTimeout * 1e-6;
opts.SUAonlyflag = false;
opts.sumcol = 0;
opts.sumflag = false;
opts.verboseflag = false;
plotflag = true;
resultdir = '';
opts.SUAregexp = '';
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'baseline'
            argnum = argnum + 1;
            opts.blspec = varargin{argnum};
        case 'copycol'
            argnum = argnum + 1;
            opts.copycol = varargin{argnum};
            if ~isnumeric(opts.copycol) || ...
                    ~isequal(size(opts.copycol), [1 2])
                error('lfp_fishAnalysis:copycol', ...
                    '''copycol'' [from to] must be a two-element numeric row vector.');
            end
        case 'endTime'
            argnum = argnum + 1;
            opts.endTime = varargin{argnum};
            if ~isnumeric(opts.endTime)
                error('lfp_fishAnalysis:endTime', ...
                    '''endTime'' must be numeric.');
            end
        case 'filenamestr'
            argnum = argnum + 1;
            filenamestr = varargin{argnum};
            if ~ischar(filenamestr)
                error('lfp_fishAnalysis:filenamestr', ...
                    '''filenamestr'' must be a string.');
            end
        case 'noplot'
            plotflag = false;
        case 'LFPonly'
            opts.LFPonlyflag = true;
        case 'MUAonly'
            opts.MUAonlyflag = true;
        case 'numSEMs'
            argnum = argnum + 1;
            opts.numSEMs = varargin{argnum};
            if ~isnumeric(opts.numSEMs)
                error('lfp_fishAnalysis:numSEMs', ...
                    '''numSEMs'' must be numeric.');
            end
        case 'read'
            opts.readflag = true;
        case 'readsmooth'
            opts.readsmoothflag = true;
        case 'save'
            opts.saveflag = true;
        case 'saveresult'
            argnum = argnum + 1;
            resultdir = varargin{argnum};
            fprintf('Saving result files to %s\n', resultdir);
        case 'savesmooth'
            opts.savesmoothflag = true;
        case 'smoothtime'
            argnum = argnum + 1;
            opts.smoothtime = varargin{argnum};
            if ~isnumeric(opts.smoothtime)
                error('lfp_fishAnalysis:smoothtime', ...
                    '''smoothtime'' must be numeric.');
            end
        case 'spikedir'
            argnum = argnum + 1;
            opts.spikedir = varargin{argnum};
            if ~ischar(opts.spikedir)
                error('lfp_fishAnalysis:spikedir', ...
                    '''spikedir'' must be a string.');
            end
        case 'SUAregexp'
            argnum = argnum + 1;
            opts.SUAregexp = varargin{argnum};
            if ~ischar(opts.SUAregexp)
                error('lfp_fishAnalysis:SUAregexp', ...
                    '''SUAregexp'' requires a string value.');
            end
        case 'SUAonly'
            opts.SUAonlyflag = true;
        case 'startTime'
            argnum = argnum + 1;
            opts.startTime = varargin{argnum};
            if ~isnumeric(opts.startTime)
                error('lfp_fishAnalysis:startTime', ...
                    '''startTime'' must be numeric.');
            end
        case 'sum'
            opts.sumflag = true;
        case 'sumcol'
            argnum = argnum + 1;
            opts.sumcol = varargin{argnum};
            if ~isnumeric(opts.sumcol)
                error('lfp_fishAnalysis:sumcol', ...
                    '''sumcol'' must be numeric.');
            end
        case 'verbose'
            opts.verboseflag = true;
        otherwise
            error('lfp_fishAnalysis:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if opts.readflag && opts.readsmoothflag
    warning('lfp_fishAnalysis:read1', ...
        'Reading smoothed waveforms, ignoring ''read'' option.');
    opts.readflag = false;
end
if opts.readflag && opts.saveflag
    warning('lfp_fishAnalysis:read2', ...
        'Reading unsmoothed waveforms, ignoring ''save'' option.');
    opts.saveflag = false;
end
if opts.readsmoothflag && (opts.saveflag || opts.savesmoothflag)
    warning('lfp_fishAnalysis:read3', ...
        'Reading smoothed waveforms, ignoring ''save'' option(s).');
    opts.saveflag = false;
    opts.savesmoothflag = false;
end

blockdirs = dir(fullfile(sessiondir, 'block*'));
if isempty(blockdirs)
    error('lfp_fishAnalysis:csc', ...
        'There are no "block" directories in %s', sessiondir);
end

prevblockdirs = [];

% LFP Analysis
if ~opts.MUAonlyflag && ~opts.SUAonlyflag
    for refidx = 1:length(tackle.localavgrefs)
        reffilename = sprintf('csc_localavgref%d.mat', refidx );
        if isempty(filenamestr)
            CSCfilespec = tackle.localavgrefs{refidx};
        else
            CSCfilespec = filenamestr;
        end
        CSCfiles{refidx} = getCSCFileList(sessiondir, ...
            CSCfilespec);
        if exist(fullfile(sessiondir, reffilename), 'file')
            % Remove the reference channel from the list to analyze:
            for k = 1:length(CSCfiles{refidx})
                if strcmpi(CSCfiles{refidx}(k).name, reffilename)
                    CSCfiles{refidx}(k) = [];
                    break
                end
            end
        else
            % There is no such reference channel file
            reffilename = '';
        end
        if isempty(CSCfiles{refidx})
            warning('lfp_fishAnalysis:csc', ...
                'There are no CSC files for localavgref #%d in %s', ...
                refidx, sessiondir);
            continue
        end
        LFPisgoodfile = true(size(CSCfiles{refidx}));
        rowidx = 0;
        for fileidx = 1:length(CSCfiles{refidx})
            if CSCfiles{refidx}(fileidx).bytes == 16384
                % skip empty files
                LFPisgoodfile(fileidx) = false;
                continue
            end
            rowidx = rowidx + 1;
            [LFPEPdata{refidx}(:,:,rowidx,:), ...
                LFPextrema{refidx}(rowidx, :), blockdirs] ...
                = analyzeOneChannel('LFP', tackle, sessiondir, ...
                CSCfiles{refidx}(fileidx).name, opts, reffilename);
            if ~isempty(prevblockdirs) && ...
                    ~isequal(blockdirs, prevblockdirs)
                warning('lfp_fishAnalysis:blockdirs1', ...
                    '<blockdirs> has changed in LFP at (%d, %d).', ...
                    fileidx, refidx);
            end
        end
        prevblockdirs = blockdirs;
        LFPfilenames{refidx} = {CSCfiles{refidx}(LFPisgoodfile).name};
    end
end

% MUA Analysis
if ~opts.LFPonlyflag &&  ~opts.SUAonlyflag
    % The apparent stupidity of the variable <spikefiles{1}> is actually
    % necessitated by the fact that LFPs come in local-reference groups,
    % and consequently they need to have file lists where each top-level
    % cel represents one group (see lfp_plotFishAnalysis).
    if isempty(filenamestr)
        namestr = 'tt*.ntt';
    else
        namestr = filenamestr;
    end
    spikefiles{1} = dir(fullfile(opts.spikedir, namestr));
    if isempty(spikefiles{1})
        warning('lfp_fishAnalysis:csc', ...
            'There are no spike files in %s', opts.spikedir);
        return
    end
    MUAisgoodfile = true(size(spikefiles{1}));
    rowidx = 0;
    for fileidx = 1:length(spikefiles{1})
        if spikefiles{1}(fileidx).bytes == 16384
            % skip empty files
            MUAisgoodfile(fileidx) = false;
            continue
        end
        rowidx = rowidx + 1;
        [MUAEPdata{1}(:,:,rowidx,:) MUAextrema{1}(rowidx,:), blockdirs] = ...
            analyzeOneChannel('MUA', tackle, sessiondir, ...
            spikefiles{1}(fileidx).name, opts);
        if ~isempty(prevblockdirs) && ...
                ~isequal(blockdirs, prevblockdirs)
            warning('lfp_fishAnalysis:blockdirs2', ...
                '<blockdirs> has changed in MUA at (%d, %d).', ...
                fileidx, refidx);
        end
    end
    MUAfilenames{1} = {spikefiles{1}(MUAisgoodfile).name};
end

% Single Unit Analysis
if ~opts.LFPonlyflag && ~opts.MUAonlyflag
    if isempty(filenamestr)
        namestr = '*.mat';
    else
        namestr = filenamestr;
    end
    spikedirlist = dir(fullfile(opts.spikedir, namestr));
    spikedirlist(cell2mat({spikedirlist.isdir})) = [];
    spikefnames = {spikedirlist.name};
    if isempty(opts.SUAregexp)
        SUAfiles = spikedirlist;
    else
        ismatch = regexp(spikefnames, opts.SUAregexp);
        SUAfiles = spikedirlist(~cellfun(@isempty, ismatch));
    end
    if isempty(SUAfiles)
        warning('lfp_fishAnalysis:csc', ...
            'There are no single unit files in %s', opts.spikedir);
        return
    end
    SUAisgoodfile = true(size(SUAfiles));
    channelNames = cell(0,1);
    rowidx = 0;
    for fileidx = 1:length(SUAfiles)
        if SUAfiles(fileidx).bytes == 16384
            % skip empty files
            SUAisgoodfile(fileidx) = false;
            continue
        end
        if opts.readsmoothflag
            % We assume there is always a baseline directory, and so it can
            % be used to find out how many clusters this tetrode has
            unitfiles = dir(fullfile(sessiondir, 'baseline', 'tt*.mat'));
            [p,basename] = fileparts(SUAfiles(fileidx).name); %#ok<ASGLU>
            unitfiles(cellfun( @isempty, ...
                regexp({unitfiles.name}, [basename 'c\d+s\d+.mat']) )) ...
                = [];
            numclusters = length(unitfiles);
        else
            numclusters = length(lfp_SpikeNames);
        end
        clustnum = 0; % first cluster is number 0
        while clustnum == 0 || clustnum < numclusters
            rowidx = rowidx + 1;
            [SUAEPdata{1}(:,:,rowidx,:) SUAextrema{1}(rowidx,:), blockdirs] = ...
                analyzeOneChannel('SUA', tackle, sessiondir, ...
                SUAfiles(fileidx).name, opts, '', clustnum+1);
            % on clustnum 0, lfp_SpikeNames was left over from the previous
            % file (or empty), so numclusters was wrong and needs updating.
            if clustnum == 0 && ~opts.readsmoothflag
                numclusters = length(lfp_SpikeNames);
            end
            if ~isempty(prevblockdirs) && ...
                    ~isequal(blockdirs, prevblockdirs)
                warning('lfp_fishAnalysis:blockdirs3', ...
                    '<blockdirs> has changed in SUA at (%d, %d).', ...
                    fileidx, refidx);
            end
            channelNames{end+1} = lfp_FileNames{lfp_ActiveFilenums(1)};
            clustnum = clustnum + 1;
        end
    end
    SUAfilenames{1} = channelNames;
end

% Plot Figures
if plotflag
    if ~opts.MUAonlyflag && ~opts.SUAonlyflag
        lfp_plotFishAnalysis(LFPEPdata, LFPextrema, LFPfilenames, opts, ...
            blockdirs);
    end
    if ~opts.LFPonlyflag && ~opts.SUAonlyflag
        lfp_plotFishAnalysis(MUAEPdata, MUAextrema, MUAfilenames, opts, ...
            blockdirs);
    end
    if ~opts.LFPonlyflag && ~opts.MUAonlyflag
        lfp_plotFishAnalysis(SUAEPdata, SUAextrema, SUAfilenames, opts, ...
            blockdirs);
    end
end

% Save results
if ~isempty(resultdir)
    p = fileparts(sessiondir);
    [p, sessionID] = fileparts(p); %#ok<ASGLU>
    if ~opts.MUAonlyflag && ~opts.SUAonlyflag
        EPdata = LFPEPdata; %#ok<NASGU>
        extrema = LFPextrema; %#ok<NASGU>
        files = CSCfiles; %#ok<NASGU>
        isgoodfile = LFPisgoodfile; %#ok<NASGU>
        save(fullfile(resultdir, sprintf('LFP%s.mat', sessionID)), ...
            'EPdata', 'extrema', 'opts', ...
            'files', 'isgoodfile', 'blockdirs');
    end
    if ~opts.LFPonlyflag && ~opts.SUAonlyflag
        EPdata = MUAEPdata; %#ok<NASGU>
        extrema = MUAextrema; %#ok<NASGU>
        files = spikefiles; %#ok<NASGU>
        isgoodfile = MUAisgoodfile; %#ok<NASGU>
        save(fullfile(resultdir, sprintf('MUA%s.mat', sessionID)), ...
            'EPdata', 'extrema', 'opts', ...
            'files', 'isgoodfile', 'blockdirs');
    end
    if ~opts.LFPonlyflag && ~opts.MUAonlyflag
        EPdata = SUAEPdata; %#ok<NASGU>
        extrema = SUAextrema; %#ok<NASGU>
        spikefiles{1} = SUAfiles;
        files = spikefiles; %#ok<NASGU>
        isgoodfile = SUAisgoodfile; %#ok<NASGU>
        save(fullfile(resultdir, sprintf('SUA%s.mat', sessionID)), ...
            'SUAfilenames', 'EPdata', 'extrema', 'opts', ...
            'files', 'isgoodfile', 'blockdirs');
    end
end
end


function [EPdata, extrema, blockdirs] = analyzeOneChannel(datatype, ...
    tackle, sessiondir, filename, opts, reffilename, clustidx)
% Note that this can still be fooled by sufficiently strong 60 Hz, i.e. if
% even after smoothing etc. the 60 Hz is large enough to produce
% significantly different artifactual local extrema on rapidly rising or
% falling phases of the physiological extrema.
%INPUTS
% filename: assumed to be a .mat file with '.mat' extension.
% reffilename: when datatype is 'LFP', this is the reference channel to
%   subtract from <filename>; if empty, no subtraction is done.  Not used
%   when datatype is 'MUA' or 'SUA'.
% clustidx: when datatype is 'SUA', contains the cluster number (i.e. index
%   into lfp_Spikes) to be analyzed.  Not used otherwise.
%OUTPUTS
% EPdata: the average and error waveform data as returned by lfp_disp, i.e.
%   in samples X {timestamp|average|error} format, with a third dimension
%   representing different blocks, s.t. EPdata(:,:,colidx) contains the
%   data returned by lfp_disp for blockdirs(colidx).
% extrema: row vector containing one element for each block, and each
%   element has a 'sigmaxTS' field and a 'sigminTS' field, each of which is
%   a (possibly empty) vector of latencies from stim to extrema that are
%   significantly different from baseline AND significantly different from
%   at least one neighboring extremum.
%NOTES
% <clustidx> is initialized to 1 when <datatype> is 'LFP' or 'MUA', and
% therefore does not need to be supplied when analyzing those data types.
global lfp_TrialIndex lfp_Samples lfp_AlignmentRef ...
    lfp_ActiveFilenums lfp_XLimAll
[p, namestem] = fileparts(filename); %#ok<ASGLU>
extrema = [];
meanbaseline = []; %#ok<NASGU>
blevtIDs = [];
blwin = [];
if isequal(opts.blspec, 'block')
    switch datatype
        case 'LFP'
            readTheCSC(fullfile(sessiondir, 'baseline'), ...
                filename, reffilename, opts, 0 );
            clustidx = 1;
        case 'MUA'
            readThe_UA('multi', fullfile(sessiondir, 'baseline'), ...
                filename, tackle, opts, 0 );
            clustidx = 1;
        case 'SUA'
            readThe_UA(clustidx, fullfile(sessiondir, 'baseline'), ...
                filename, tackle, opts, 0 );
        otherwise
            warning('lfp_fishAnalysis:datatype', ...
                'Unrecognized datatype, skipping');
            return
    end
    
    % Here we take advantage of the fact that there is absolutely "nothing
    % happening" during the baseline period to compute a single mean of ALL
    % the samples rather than a mean waveform aligned on meaningless
    % non-events:
    meanbaseline = nanmean(lfp_Samples{lfp_ActiveFilenums(1)}( ...
        lfp_TrialIndex(1,3):lfp_TrialIndex(end,4) ));
elseif isequal(opts.blspec, 'none')
    meanbaseline = 0;
elseif iscell(opts.blspec) && length(opts.blspec) == 2 ...
        && length(opts.blspec{2}) == 2
    blevtIDs = opts.blspec{1};
    blwin = opts.blspec{2};
    meanbaseline = 0;
else
    error('lfp_fishAnalysis:blspec', ...
        'Unrecognized format for ''baseline'', <blspec>');
end
if opts.verboseflag
    fprintf('Found baseline for %s\n', namestem);
end

blockdirs = dir(fullfile(sessiondir, 'block*'));
ncols = length(blockdirs);
for colnum = 1:ncols
    % <colnum> refers to the array of result plots, and is distinguished
    % from <blocknum>, which is extracted from the actual name of the block
    % directory and blocknum might be different from colnum (e.g. if there
    % is a block missing):
    blocknum =  str2double(regexprep(blockdirs(colnum).name, ...
        '^block(\d+)$', '$1'));
    if blockdirs(colnum).isdir && blocknum > 0
        switch datatype
            % The readThe* functions subtract <meanbaseline>:
            case 'LFP'
                readTheCSC(fullfile(sessiondir, blockdirs(colnum).name), ...
                    filename, reffilename, opts, meanbaseline );
            case 'MUA'
                readThe_UA('multi', fullfile(sessiondir, blockdirs(colnum).name), ...
                    filename, tackle, opts, meanbaseline );
            case 'SUA'
                readThe_UA(clustidx, fullfile(sessiondir, blockdirs(colnum).name), ...
                    filename, tackle, opts, meanbaseline );
            otherwise
                warning('lfp_fishAnalysis:datatype2', ...
                    'Unrecognized datatype, skipping');
                return
        end
        if opts.verboseflag
            fprintf('Read data for %s block %d\n', namestem, blocknum);
        end
        lfp_AlignmentRef = tackle.blkIDs(blocknum);
        lfp_XLimAll = [opts.startTime opts.endTime];
        lfp_selectByDuration;
        [junk, data] = lfp_disp([], lfp_ActiveFilenums(1), ...
            [], 'avg', 'err2', 'noplot'); %#ok<ASGLU>
        if colnum == 1
            EPdata = NaN(size(data{1},1), size(data{1},2), ncols);
            EPdata(:,:,1) = data{1};
        else
            if size(data{1},1) == size(EPdata,1) && ...
                    size(data{1},2) == size(EPdata,2)
                EPdata(:,:,colnum) = data{1};
            else
                if size(data{1},2) ~= size(EPdata,2)
                    error('lfp_fishAnalysis:oops', ...
                        'This is logically impossible.');
                end
                testidx = round(size(EPdata,1)/2);
                [v, matchidx] = min( ...
                    abs(data{1}(:,1) - EPdata(testidx,1)) ); %#ok<ASGLU>
                if size(data{1},1) > size(EPdata,1)
                    warning('lfp_fishAnalysis:trim1', ...
                        'Trimming %s colnum %d to match previous lengths', ...
                        filename, colnum);
                    data{1} = data{1}( ...
                        (1:size(EPdata,1)) + matchidx - testidx, : );
                    EPdata(:,:,colnum) = data{1};
                else
                    warning('lfp_fishAnalysis:trim2', ...
                        'Trimming previous columns'' data to match %s colnum %d', ...
                        filename, colnum);
                    EPdata = EPdata( ...
                        (1:size(data{1},1)) + testidx - matchidx, :, : );
                    EPdata(:,:,colnum) = data{1};
                end
            end
        end
        vmaxidx = dg_findFlattops(EPdata(:,2,colnum));
        vminidx = dg_findFlattops(-EPdata(:,2,colnum));
        
        % There necessarily must be a series of alternating mins and maxes.
        % Find the ones that are statistically significantly different with
        % respect to at least one neighbor. The first and last of the
        % series can only be significant with respect to one neighbor, but
        % they can also be significant with respect to the start or end
        % samples respectively.  The internal elements can be significant
        % with respect to neighbors on either side.
        if isempty(vmaxidx) && isempty(vminidx)
            warning('fishAnal:noextrema', ...
                '%s has no extrema in block %d', namestem, blocknum);
            extrema(1,blocknum).sigmaxTS = []; %#ok<*AGROW>
            extrema(1,blocknum).sigminTS = [];
            continue
        else
            % Compute significance between each pair of extrema and
            % start or end sample on outside elements.  For this
            % purpose, we define significance as "the mean at the
            % second point is outside the CLs at the first point".
            maxIsFirst = false;
            if isempty(vminidx) || ~isempty(vmaxidx) && vmaxidx(1) < vminidx(1)
                maxIsFirst = true;
            end
            extrmidx = sort([1 reshape(vmaxidx,1,[]) ...
                reshape(vminidx,1,[]) size(EPdata,1)]);
            % Since we are including the start and end samples in the
            % list of points <extrmidx>, there is one more significance
            % evaluation than there are extrema.  To prevent headaches, we
            % make the <issig> array the same size as <extrmidx> and define
            % the first element to be false.
            issig = false(size(extrmidx));
            issig(2:end) = abs( EPdata(extrmidx(2:end),2,colnum) ...
                - EPdata(extrmidx(1:end-1),2,colnum) ) ...
                > opts.numSEMs*EPdata(extrmidx(1:end-1),3,colnum)/2;
            % The first extremum in each significant comparison has the
            % same claim to fame as the second, so we truncate, shift and
            % OR. This result in a logical index the same length as
            % <extrmidx>, which gives us <sigidx2> as an index into an
            % index into <data> (hence the "2" in its name).
            sigidx2 = find(issig(1:end-1) | issig(2:end));
            % But please do not report the first sample as a significant
            % extremum even if it is significantly different from the first
            % actual extremum:
            if ~isempty(sigidx2) && sigidx2(1) == 1
                sigidx2(1) = [];
            end
            % Throw away any that are not significantly different from
            % baseline, which we here define as the ones whose CLs include
            % the baseline value.  If blspec is not 'block', this is where
            % the infinite horror of noisy baselines rears its ugly head
            % and requires us to calculate new CLs and to subtract the new
            % baseline from the EP.  In the case where blspec is 'block',
            % then <meanbaseline> has already been subtracted from the
            % channels in lfp_Samples, and we are looking for the ones
            % where the CLs include zero.
            SEM = EPdata(extrmidx(sigidx2),3,colnum)/2;
            if isempty(blevtIDs)
                refvalue = 0;
            else
                lfp_AlignmentRef = blevtIDs;
                lfp_XLimAll = blwin;
                lfp_selectByDuration;
                [junk, data] = lfp_disp([], lfp_ActiveFilenums(1), ...
                    [], 'ovr', 'noplot'); %#ok<ASGLU>
                blvals = data(:);
                refvalue = mean(blvals);
                blSEMsqrd = var(blvals) / numel(blvals);
                SEM = sqrt(blSEMsqrd * ones(size(SEM)) + SEM .* SEM);
                EPdata(:,2,colnum) = EPdata(:,2,colnum) - refvalue;
            end
            CLwidth = opts.numSEMs * SEM;
            badones = abs(EPdata(extrmidx(sigidx2),2,colnum)) ...
                <= CLwidth;
            sigidx2(badones) = [];
            % Now convert back to separate lists of min and max, recalling
            % that the first element in <extrmidx> is actually the first
            % sample and not an extremum, so <maxIsFirst> implies that
            % EVEN-numbered entries in <extrmidx> are maxes.
            if maxIsFirst
                sigmaxidx2 = sigidx2(mod(sigidx2, 2) == 0);
                sigminidx2 = sigidx2(mod(sigidx2, 2) > 0);
            else
                sigmaxidx2 = sigidx2(mod(sigidx2, 2) > 0);
                sigminidx2 = sigidx2(mod(sigidx2, 2) == 0);
            end
            % Meta-post-finally, we retrieve the relative times of the
            % significant extrema:
            extrema(1,blocknum).sigmaxTS = ...
                EPdata(extrmidx(sigmaxidx2),1,colnum); %#ok<*AGROW>
            extrema(1,blocknum).sigminTS = ...
                EPdata(extrmidx(sigminidx2),1,colnum);
        end
    end
end
if opts.verboseflag
    fprintf('Finished Analysis for %s\n', namestem);
end
end


function readTheCSC(presetdir, filename, reffilename, opts, baseline)
global lfp_ActiveFilenums lfp_FileNames lfp_SamplePeriod lfp_Samples
if opts.readflag
    [p, name, ext] = fileparts(filename); %#ok<ASGLU>
    filename = [name '_clean' ext];
elseif opts.readsmoothflag
    [p, name, ext] = fileparts(filename); %#ok<ASGLU>
    filename = sprintf('%s_clean%d%s', name, opts.smoothtime, ext);
end
files2read = {'lfp_fishstory.mat'};
files2read{end+1} = filename;
if  ~isempty(reffilename)
    files2read{end+1} = reffilename;
end
lfp_read2('preset', presetdir, files2read);
if ~(opts.readflag || opts.readsmoothflag)
    if ~isempty(reffilename)
        lfp_createWave(@lfp_wavediff, lfp_ActiveFilenums(1:2), ...
            'replace', lfp_ActiveFilenums(1), 'name', sprintf( ...
            '%s_clean', lfp_FileNames{lfp_ActiveFilenums(1)} ));
    end
    lfp_createWave(@lfp_rmlnspec, lfp_ActiveFilenums(1), 60, ...
        'replace', lfp_ActiveFilenums(1), ...
        'name', lfp_FileNames{lfp_ActiveFilenums(1)} );
    lfp_Samples{lfp_ActiveFilenums(1)} = ...
        lfp_Samples{lfp_ActiveFilenums(1)} - baseline;
    if opts.saveflag
        lfp_save_noGUI(lfp_FileNames{lfp_ActiveFilenums(1)}, 'csc');
    end
end
if ~opts.readsmoothflag
    lfp_createWave(@lfp_wavesmooth, lfp_ActiveFilenums(1), ...
        round(opts.smoothtime*1e-3/lfp_SamplePeriod), ...
        'replace', lfp_ActiveFilenums(1), 'name', sprintf('%s%d', ...
        lfp_FileNames{lfp_ActiveFilenums(1)}, opts.smoothtime ));
    if opts.savesmoothflag
        lfp_save_noGUI(lfp_FileNames{lfp_ActiveFilenums(1)}, 'csc');
    end
end
end


function readThe_UA(mode, cscdir, filename, tackle, opts, baseline)
% mode:  if 'multi', then treat <filename> as containing multiunit
%   activity, i.e. not clustered (or ignore the cluster IDs), so this is
%   just one channel of data.  If numeric, it designates the <clustidx>
%   from <filename> to convert to spike density function.
global lfp_ActiveFilenums lfp_FileNames lfp_SamplePeriod lfp_Samples ...
    lfp_SpikeNames lfp_TimeStamps lfp_SamplesPerFrame lfp_Events
if isequal(mode, 'multi')
    clustidx = 1;
else
    clustidx = mode;
end
files2read = {'lfp_fishstory.mat'};
[p, name] = fileparts(filename); %#ok<ASGLU>
if opts.readsmoothflag
    [p, name] = fileparts(filename); %#ok<ASGLU>
    if isequal(mode, 'multi')
        files2read{end+1} = sprintf('%s_MUA%d.mat', name, opts.smoothtime);
    else
        % tt17C6s5
        files2read{end+1} = sprintf('%sc%ds%d.mat', name, clustidx-1, ...
            opts.smoothtime);
    end
else
    % In order to create a spike density function CSC channel, we must find
    % and read a token CSC file to give us the CSC timestamp structure.
    % Reconciling timestamps is not an issue here because reconciliation
    % will be done by lfp_read2 when the newly created CSC files are
    % actually used.
    cscfilename = '';
    if ~isempty(tackle.localavgrefs)
        for refidx = 1:length(tackle.localavgrefs)
            reffilename = sprintf('csc_localavgref%d.mat', refidx );
            if exist(fullfile(cscdir, reffilename), 'file')
                cscfilename = reffilename;
                break
            else
                CSCfiles = getCSCFileList(cscdir, ...
                    tackle.localavgrefs{refidx});
                if ~isempty(CSCfiles)
                    cscfilename = CSCfiles(1).name;
                end
            end
        end
    end
    if isempty(cscfilename)
        CSCfiles = dir(fullfile(cscdir, 'csc_*'));
        if isempty(CSCfiles)
            CSCfiles = dir(fullfile(cscdir, 'lfp_*'));
        end
        if isempty(CSCfiles)
            CSCfiles = dir(fullfile(cscdir, '*.ncs'));
        end
        if isempty(CSCfiles)
            warning('lfp_fishAnalysis:noCSC', ...
                'Could not find CSC file to set timestamps in %s', ...
                cscdir);
            return
        else
            cscfilename = CSCfiles(1).name;
        end
    end
    files2read{end+1} = cscfilename;
end
lfp_read2('preset', cscdir, files2read);
if ~opts.readflag && ~opts.readsmoothflag && isequal(mode, 'multi')
    % Create the 'Single Cut Cluster (*.MAT)' MUA file.
    TS = 1e-6 * dg_readSpike(fullfile(opts.spikedir, filename));
    MUAfile = [name '_MUA.mat'];
    % delete spikes that are outside of this directory's time range
    TSmin = min(lfp_Events(1,1), lfp_TimeStamps(1));
    TSmax = max(lfp_Events(end,1), ...
        lfp_TimeStamps(end) + lfp_SamplePeriod * lfp_SamplesPerFrame);
    TS(TS < TSmin | TS > TSmax) = []; %#ok<NASGU>
    save(fullfile(opts.spikedir, MUAfile), 'TS');
    clear('TS');
end
if ~opts.readsmoothflag
    % Annoyingly, we need to lfp_add the spike file as the opposite type of
    % file to what we are analyzing, i.e. to do multiunit analysis, we load
    % the data AS IF it were a single cluster, whereas to do SINGLE unit
    % analysis, we must analyze each cluster separately, and thus load the
    % file as what it truly is, namely a Multiple Cut Cluster file.
    if isequal(mode, 'multi')
        lfp_add('preset', opts.spikedir, {MUAfile}, ...
            'Single Cut Cluster (*.MAT)', false);
        namefmt = '%s%d';
    else
        lfp_add('preset', opts.spikedir, {filename}, ...
            'Multiple Cut Cluster (*.MAT)', false);
        namefmt = '%ss%d';
    end
    % convert the spike train(s) to CSC(s), smooth, and
    % convert to Hz.
    lfp_createWave( @lfp_spike2wave, clustidx, 'multispike', ...
        'replace', lfp_ActiveFilenums(1) );
    lfp_createWave(@lfp_wavesmooth, lfp_ActiveFilenums(1), ...
        round(opts.smoothtime*1e-3/lfp_SamplePeriod), ...
        'replace', lfp_ActiveFilenums(1), 'name', sprintf(namefmt, ...
        lfp_SpikeNames{clustidx}, opts.smoothtime ));
    lfp_Samples{lfp_ActiveFilenums(1)} = ...
        lfp_Samples{lfp_ActiveFilenums(1)} / lfp_SamplePeriod;
    lfp_SamplesUnits{lfp_ActiveFilenums(1)} = 'Hz'; %#ok<NASGU>
    lfp_Samples{lfp_ActiveFilenums(1)} = ...
        lfp_Samples{lfp_ActiveFilenums(1)} - baseline;
    if opts.savesmoothflag
        lfp_save_noGUI(lfp_FileNames{lfp_ActiveFilenums(1)}, 'csc', ...
            'destdir', cscdir);
    end
end
end


function CSCfiles = getCSCFileList(sessiondir, CSCfilespec)
%INPUTS
% sessiondir: directory to search
% CSCfilespec: OS-compatible filename matching string
%OUTPUT
% CSCfiles: as returned by calling matlab 'dir' function
if iscell(CSCfilespec)
    CSCfiles = [];
    for specidx = 1:length(CSCfilespec)
        CSCfiles = [ CSCfiles
            dir(fullfile(sessiondir, CSCfilespec{specidx}))]; %#ok<AGROW>
    end
else
    CSCfiles = dir(fullfile(sessiondir, CSCfilespec));
end
end

