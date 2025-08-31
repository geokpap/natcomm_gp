function [histcnts, histbinctrs, nsuprathr] = ...
    dg_Simon2Mat(infilename, cscfilename, ...
    evtfilename, varargin)
% Version 2: ignores all the timestamp estimates except for purposes of
%   calculating <medFs>.  Constructs synthetic timestamps on the assumption
%   that there are absolutely no samples missing and that Fs is absolutely
%   fixed at medFs.
%INPUTS
% infilename: absolute or relative path to Simon-format input file.
% cscfilename, evtfilename: absolute or relative path to dg_Nlx2Mat-format
%   output files.
%OUTPUTS
% Output is to the two output files, plus:
% histcnts: N x 2 cell array containing the histogram returned from
%   dg_S2M_findvalley.  Column 1: half-wave rectified detrended sample
%   values; column 2: differenced detrended sample values.  N is the number
%   of trials in the file.  Empty if dg_S2M_findvalley was not run.
% histbinctrs: N x 2 cell array containing the center values of the bins
%   returned from dg_S2M_findvalley; matches histcnts cell-for-cell.
% nsuprathr: N x 3 numeric array containing the values of numsampgt,
%   numdiffgt, numbothgt; the value is NaN for counts that were not
%   computed.
%OPTIONS
% 'calcstimTS' - invokes code to calculate the actual stimTS based on the
%   stimulation artifact waveform.
% 'estFs' - invokes code to estimate the actual sampling frequency based on
%   the numbers of samples recorded and the estimated stimTS values in the
%   file.
% 'usebothTS' - same as 'calcstimTS', except that the estimated stimTS
%   values in the file are used to determine a 20 ms long interval to use
%   for the entire stimulation artifact waveform analysis ('calcstimTS'
%   uses the entire trial).
% 'verbose' - yup, eh?
%NOTES
% The Simon format contains LFP data for one channel, broken up into
% records of roughly 1 s long each.  Each record comprises:
%   (int) LFP ID
%   (int) Sampling Rate from DAQ card's spec sheet
%   (long long) estimated CPU clock time of Stimulation in microseconds
%   (long long) estimated CPU clock time of LFP signal start in microseconds
%   (long long) estimated time of LFP signal end
%   (int) N = total number of sample points
%   (int) Number of samples before Stimulation
%   (int) Future use
%   (int)*N The N sample values, one (int) per sample
%  Each record provides a range of samples that are contiguous with its
% predecessor, i.e. the continuous recording can be reconstructed by
% concatening the sample values from each record.  The actual sampling rate
% can only be estimated by computing the actual elapsed time between the
% stimTS timestamps of at least two successive records that have stimTS
% separated by more than 0.2 s (to eliminate cases where the time
% difference is negative and cases where the time difference is small
% enough that the irreducible Windows jitter introduces an error greater
% than 1%).  The "Number of samples before Stimulation" field is also
% subject to jitter.
%   The trial start event marker is TTL 1, the Stim marker is TTL 2, the
% trial end event is TTL 3.

%$Rev: 258 $
%$Date: 2017-06-21 17:36:42 -0400 (Wed, 21 Jun 2017) $
%$Author: dgibson $

calcstimTSflag = false;
estFsflag = false;
usebothTSflag = false;
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
        case 'calcstimTS'
            calcstimTSflag = true;
        case 'estFs'
            estFsflag = true;
        case 'usebothTS'
            usebothTSflag = true;
        case 'verbose'
            verboseflag = true;
        otherwise
            error('dg_Simon2Mat:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end


% Read the file

fid = fopen(infilename);
if isequal(fid, -1)
    error('dg_Simon2Mat:open', ...
        'Could not open file %s', infilename);
end
histcnts = {};
histbinctrs = {};
nsuprathr = zeros(0, 3);
try
    % Unpack binary data into a struct vector:
    recs = struct('LFP_ID', {}, 'Fs', {}, 'stimTS', {}, 'startTS', {}, ...
        'end_est', {}, 'N', {}, 'preN', {}, 'mt', {}, 'samples', {});
    idx = 1;
    LFP_ID = fread(fid, 1, 'int');
    while ~isempty(LFP_ID)
        recs(idx).LFP_ID = LFP_ID;
        recs(idx).Fs = fread(fid, 1, 'int');
        recs(idx).stimTS = fread(fid, 1, 'bit64');
        recs(idx).startTS = fread(fid, 1, 'bit64');
        recs(idx).end_est = fread(fid, 1, 'bit64');
        recs(idx).N = fread(fid, 1, 'int');
        recs(idx).preN = fread(fid, 1, 'int');
        recs(idx).mt = fread(fid, 1, 'int');
        recs(idx).samples = reshape(fread(fid, recs(idx).N, 'int'), 1, []);
        LFP_ID = fread(fid, 1, 'int');
        
        if ~(calcstimTSflag || usebothTSflag)
            recs(idx).actual_preN = recs(idx).preN;
            nsuprathr(end+1, :) = NaN(1, 3);
        else
            % Find the first stimulation artifact in the trial.  The strategy
            % is to compute a histogram of sample values and one of
            % inter-sample slopes, find a suitable valley to use as a
            % threshold, and ...
            % For each of <numsampgt> and <numdiffgt>, there are three
            % possibilities:  too small (cannot use), just right (can use
            % solo), too large (try combining both).
            if usebothTSflag
                % confine analysis time window to +/- 10 ms around the
                % nominal stim time specified by preN:
                numpts = round(0.01 * recs(idx).Fs);
                startidx = max(1, recs(idx).preN + 1 - numpts);
                endidx = recs(idx).preN + 1 + numpts;
            else
                % use the whole trial:
                startidx = 1;
                endidx = numel(recs(idx).samples);
            end
            detrended = detrend(recs(idx).samples(startidx:endidx));
            numsampgt = NaN;
            numdiffgt = NaN;
            numbothgt = NaN;
            % An older formula for <sampthresh> was:
            %    sampvalthresh = prctile(detrended, 99.95);
            detrendedgt0 = detrended;
            detrendedgt0(detrended <= 0) = 0;
            % analyze the actual sample values, setting <sampstatus>
            % according to where we end up in the code
            [sampthresh, histcnts{end+1, 1}, histbinctrs{end+1, 1}] = ...
                dg_S2M_findvalley(detrendedgt0); %#ok<AGROW>
            if isempty(sampthresh)
                sampstatus = 'bad';
            else
                numsampgt = sum(detrended > sampthresh);
                if numsampgt >= 12
                    if numsampgt < 32
                        % Just right!'
                        sampstatus = 'good';
                        recs(idx).actual_preN = ...
                            find(detrended > sampthresh, 1 ) ...
                            + startidx - 1;
                    else
                        % Too many
                        sampstatus = 'large';
                    end
                else
                    % Too few
                    sampstatus = 'small';
                end
            end
            % analyze the first difference of <detrended>, setting
            % <diffstatus> according to where we end up in the code
            diffstatus = ''; % initialize
            if ~isequal(sampstatus, 'good')
                diffsamples = diff(detrended);
                diffsamples(diffsamples < 0) = 0;
                [diffthresh, histcnts{end+1, 2}, histbinctrs{end+1, 2}] = ...
                    dg_S2M_findvalley(diffsamples); %#ok<AGROW>
                if isempty(diffthresh)
                    diffstatus = 'bad';
                else
                    numdiffgt = sum(diffsamples > diffthresh);
                    if numdiffgt >= 3
                        if numdiffgt < 20
                            % Just right!'
                            diffstatus = 'good';
                            recs(idx).actual_preN = ...
                                find(diffsamples > diffthresh, 1 ) ...
                                + startidx - 1;
                        else
                            % Too many
                            diffstatus = 'large';
                        end
                    else
                        % Too few
                        diffstatus = 'small';
                    end
                end
            end
            % If actual_preN has not yet been successfully calculated as
            % reported by <sampstatus, diffstatus>, try ANDing the two
            % crtieria, and if that also fails then mark the trial "bad"
            if ~any(ismember({sampstatus, diffstatus}, 'good'))
                if all(ismember({sampstatus, diffstatus}, 'bad')) ...
                        || any(ismember({sampstatus, diffstatus}, 'small'))
                    % ANDing them isn't possible; this is a "bad" trial and
                    % must be marked as such.
                    recs(idx).actual_preN = 0;
                else
                    numbothgt = sum( detrended(2:end) > sampthresh ...
                        & diffsamples > diffthresh );
                    if numbothgt < 3 || numbothgt > 20
                        % mark trial as "bad"
                        recs(idx).actual_preN = 0;
                    else
                        recs(idx).actual_preN = find( detrended > sampthresh ...
                            & diffsamples > diffthresh ) + startidx - 1;
                    end
                end
            end
            nsuprathr(end+1, :) = [numsampgt, numdiffgt, numbothgt]; %#ok<AGROW>
        end
        if verboseflag
            fprintf('idx: %d; {sampstatus diffstatus}: {%s %s}\n', ...
                idx, sampstatus, diffstatus);
            fprintf( '{numsampgt, numdiffgt, numbothgt}: %s\n', ...
                dg_thing2str({numsampgt, numdiffgt, numbothgt}) );
        end
        idx = idx + 1;
    end
    fclose(fid);
catch e
    fclose(fid);
    rethrow(e);
end
fprintf( 'Read %s; median sample per record = %d\n', infilename, ...
    median([recs.N]) );

if estFsflag
    % Estimate Fs for each pair of records:
    Fs = NaN(length(recs) - 1, 1); % sample freq in MHz between idx and idx+1
    deltaT = NaN(length(recs) - 1, 1); % sample freq in MHz
    numsamp = NaN(length(recs) - 1, 1); % number of samples between timestamps
    for idx = 1:length(Fs)
        deltaT(idx) = recs(idx+1).stimTS - recs(idx).stimTS;
        numsamp(idx) = recs(idx).N - recs(idx).preN + recs(idx+1).preN;
        % Compute Fs only from pairs of trials that meet minimum duration;
        % leave the rest as NaN:
        if deltaT(idx) > 0.2
            Fs(idx) = numsamp(idx) / deltaT(idx);
        end
    end
    
    % Calculate the median Fs and Tsamp:
    nanFs = isnan(Fs);
    if all(nanFs)
        error('dg_Simon2Mat:noFs', ...
            'There are no pairs of trials that can be used to estimate Fs.');
    end
    if sum(nanFs) > length(Fs) / 2
        warning('dg_Simon2Mat:funkyFs', ...
            'Fs could not be estimated from most pairs of trials.');
    end
    coarsefac = 1.5;
    coarseFs = median(Fs(~nanFs));
    medFs = median(Fs( ...
        Fs > coarseFs/coarsefac & Fs < coarseFs*coarsefac & ~nanFs ));
else
    medFs = recs(1).Fs * 1e-6; % given in Hz, but has to be in MHz
end

Tsamp = 1/medFs;

% Pack the samples into the frame-formatted arrays dg_Nlx2Mat_Samples and
% dg_Nlx2Mat_Timestamps.
totnumsamp = sum([recs.N]);
framesize = 512;
dg_Nlx2Mat_SamplesUnits = 'AD'; %#ok<NASGU>
dg_Nlx2Mat_Samples = NaN(framesize, ceil(totnumsamp/framesize));
dg_Nlx2Mat_Samples(1:totnumsamp) = [recs.samples];
sampTS = NaN(size(dg_Nlx2Mat_Samples));
% <sampTS> must not start at 0 because lfp_getTrialID assumes it doesn't.
sampTS(1:totnumsamp) = recs(1).startTS + (0 : totnumsamp - 1) * Tsamp;
dg_Nlx2Mat_Timestamps = round(sampTS(1, :));

% Create evtTS, dg_Nlx2Mat_TTL, dg_Nlx2Mat_EventStrings, inserting trial
% start and trial end events in the process.  If recs(idx).actual_preN is
% 0, that means the stim artifacts were too hard to find, so no stim event
% is recorded.
evtTS = [];
dg_Nlx2Mat_TTL = [];
dg_Nlx2Mat_EventStrings = {};
trialstartidx = 1;
for idx = 1:length(recs)
    % append this trial's start, stim, and end times:
    if recs(idx).actual_preN ~= 0
        evtTS = [ evtTS ...
            sampTS(trialstartidx) ...
            sampTS(trialstartidx + recs(idx).actual_preN) ...
            sampTS(trialstartidx + recs(idx).N - 1) ]; %#ok<AGROW>
        dg_Nlx2Mat_TTL = [ dg_Nlx2Mat_TTL ...
            1 ...
            2 ...
            3 ];  %#ok<AGROW>
        dg_Nlx2Mat_EventStrings = [dg_Nlx2Mat_EventStrings
            { 'RecID: 4098 Port: 0 TTL Value: 0x8001'
            'RecID: 4098 Port: 0 TTL Value: 0x8002'
            'RecID: 4098 Port: 0 TTL Value: 0x8003' } ]; %#ok<AGROW>
    else
        evtTS = [ evtTS ...
            sampTS(trialstartidx) ...
            sampTS(trialstartidx + recs(idx).N - 1) ]; %#ok<AGROW>
        dg_Nlx2Mat_TTL = [ dg_Nlx2Mat_TTL ...
            1 ...
            3 ];  %#ok<AGROW>
        dg_Nlx2Mat_EventStrings = [dg_Nlx2Mat_EventStrings
            { 'RecID: 4098 Port: 0 TTL Value: 0x8001'
            'RecID: 4098 Port: 0 TTL Value: 0x8003' } ]; %#ok<AGROW>
    end
    trialstartidx = trialstartidx + recs(idx).N;
end

% The last frame probably contains NaNs.  If so, truncate:
if any(isnan(dg_Nlx2Mat_Samples(:, end)))
    dg_Nlx2Mat_Samples(:, end) = [];  %#ok<NASGU>
    dg_Nlx2Mat_Timestamps(end) = []; %#ok<NASGU>
end

% Save files.
save(cscfilename, 'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_Timestamps', ...
    'dg_Nlx2Mat_SamplesUnits');
dg_Nlx2Mat_Timestamps = evtTS;  %#ok<NASGU>
save(evtfilename, 'dg_Nlx2Mat_Timestamps', 'dg_Nlx2Mat_TTL', ...
    'dg_Nlx2Mat_EventStrings');



