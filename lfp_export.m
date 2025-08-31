function lfp_export(filename, filenums, varargin)
%lfp_export(filename, filenums, varargin)
% Saves trial-formatted wave data and events as a MAT-file, containing
% variables that are cell arrays with one cell per trial:
%   <wavedata> - each cell contains an array of timepoints x channels, one
%       channel per saved filenum.  Use <filenames> to map channels to file
%       names (which will often be in lexical rather than numerical order).
%   <events> - each cell contains a two-column array with timestamps
%       in seconds relative to the first sample time in <wavedata>.
% and also one scalar variable:
%   <sampleperiod> in seconds
% and one cell string array:
%   <filenames> - cell <k> contains the filename corresponding to column
%       <k> of <wavedata>.
% If <filenums> is empty, all filenums are saved.  <varargin> may contain
% lfp_export options and/or options for the "save" function; any options
% recognized as lfp_export options are removed from <varargin>, and the
% rest are passed on to the final "save" call.
%OPTIONS
% 'selectedTrials' - saves only trials that have lfp_enabledTrials(trial) =
%   true.
% 'useXLim' - calls lfp_findCommonTime to truncate wave data to a single
%   common time interval that exists for all trials, and truncates further
%   if necessary to confine the common interval to that specified by
%   lfp_XLimAll.  This has the result that the length of the wavedata is
%   the same for every trial.
% 'leader', leadtime - instead of beginning at the start of the trial, each
%   trial's data begin at a point <leadtime> seconds before the start of
%   trial.  If a given trial does not have data recorded that far back, a
%   warning is issued and the trial is skipped.
% 'numeric' - sorts <filenums> in numeric order before saving.
%
% Any options other than the above are passed through to the Matlab 'save'
% function.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 2 || isempty(filenums)
    filenums = lfp_ActiveFilenums;
end
if (length(filenums) ~= numel(filenums))
    error('<filenums> must be a vector');
end
argnum = 1;
args2delete = [];
trials = 1:size(lfp_TrialIndex, 1);
leadtime = 0;
numericflag = false;
xlimflag = false;
while argnum <= length(varargin)
    if isequal(class(varargin{argnum}), 'char')
        switch varargin{argnum}
            case 'leader'
                args2delete(end+1) = argnum;
                argnum = argnum + 1;
                args2delete(end+1) = argnum;
                leadtime = varargin{argnum};
            case 'numeric'
                numericflag = true;
                args2delete(end+1) = argnum;
            case 'selectedTrials'
                trials = lfp_enabledTrials(find(lfp_SelectedTrials));
                args2delete(end+1) = argnum;
            case 'useXLim'
                if isempty(lfp_AlignmentRef)
                    error('lfp_AlignmentRef must not be empty with ''useXLim''');
                end
                xlimflag = true;
                args2delete(end+1) = argnum;
        end
    end
    argnum = argnum + 1;
end
varargin(args2delete) = [];

if numericflag
    filenums = sort(filenums);
end

for fidx = 1:length(filenums)
    filenames{fidx} = lfp_FileNames{filenums(fidx)};
end

if xlimflag
    [interval, rawtrialinfo] = lfp_findCommonTime(trials);
    xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
    if ~isempty(lfp_XLimAll)
        interval(1) = max(interval(1), xlimpoints(1));
        interval(2) = min(interval(2), xlimpoints(2));
    end
    badtrialidx = find(rawtrialinfo(:,3) == 0);
    if ~isempty(badtrialidx)
        warning('lfp_export:badtrials', ...
            'Skipping trials with no ref event:\n%s', ...
            dg_canonicalSeries(trials(badtrialidx)) );
        trials(badtrialidx) = [];
        rawtrialinfo(badtrialidx,:) = [];
    end
end

for trialidx = 1:length(trials)
    trial = trials(trialidx);
    if leadtime == 0
        start4export = lfp_TrialIndex(trial,1);
        eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
    else
        start4export = lfp_Events(lfp_TrialIndex(trial,1),1) - leadtime;
        evtidx = find(lfp_Events(1:lfp_TrialIndex(trial,1),1) ...
            > start4export );
        if length(evtidx) == 0
            warning('lfp_export:noevents', ...
                'Trial %d contained no events, skipping.', trial );
            continue
        end
        eventrange = evtidx(1) : lfp_TrialIndex(trial,2);
    end
    events{trialidx} = lfp_Events(eventrange,:);
    if xlimflag
        [samplerange, refpoint] = lfp_getSampleRange(trial);
    else
        if leadtime == 0
            samplerange = lfp_TrialIndex(trial,3):lfp_TrialIndex(trial,4);
        else
            startsample = lfp_time2index(start4export);
            if abs(lfp_index2time(startsample) - start4export) ...
                    / lfp_SamplePeriod > .51
                warning('lfp_export:insuffleader', ...
                    'Trial %d did not contain enough leader data, skipping.', ...
                    trial );
                continue
            end
            samplerange = startsample:lfp_TrialIndex(trial,4);
        end
    end
    if isempty(samplerange)
        warning('lfp_export:nodata', ...
            'No wave data were selected for trial %d', trial );
        wavedata{trialidx} = [];
    else
        reftime = lfp_index2time(samplerange(1));
        events{trialidx}(:,1) = events{trialidx}(:,1) - reftime;
        wavedata{trialidx} = zeros(length(samplerange), length(filenums));
        for fidx = 1:length(filenums)
            wavedata{trialidx}(:,fidx) = lfp_Samples{filenums(fidx)}(samplerange);
        end
    end
end

sampleperiod = lfp_SamplePeriod;
save(filename, 'wavedata', 'events', 'sampleperiod', 'filenames', ...
    varargin{:});
