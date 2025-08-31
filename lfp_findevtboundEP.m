function ts = lfp_findevtboundEP(filenum, evtbounds, offsets, EP, varargin)
%lfp_findevtboundEP(filenum, evtbounds, offsets, EP); for lfp_createEvents
% If more than one lag results in the same maximum xcov, then the first one
% is used (as determined by Matlab 'max' function).
%OPTIONS
% 'fishbad' - adds any trials that provoke warning
%   'lfp_findevtboundEP:fish' to lfp_BadTrials.
% 'logxcovs' - log 'em
% 'minxcov', minxcov - only marks the trial if the value of the normalized
%   xcov at the best lag is greater than or equal to <minxcov>, the
%   normalization factor being sqrt(EPvar*trialvar) * length(EP), where
%   EPvar and trialvar are the variances of the EP and the selected trial
%   data respectively. Note that this normalization can result in xcovs
%   that are greater than 1 if the trial data have a greater amplitude than
%   the EP in the region of best match, or if the trial data have a lower
%   variance than the EP.  Normalized xcovs below 0.5 tend to suck.
% 'nofish' - do not mark trials that provoke warning 
%   'lfp_findevtboundEP:fish'
% 'plothist' - histogram 'em

%$Rev: 126 $
%$Date: 2010-06-01 19:14:25 -0400 (Tue, 01 Jun 2010) $
%$Author: dgibson $

lfp_declareGlobals;

fishbadflag = false;
fishflag = true;
logxcovsflag = false;
minxcov = [];
plothistflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'fishbad'
            fishbadflag = true;
        case 'minxcov'
            argnum = argnum + 1;
            minxcov = varargin{argnum};
        case 'nofish'
            fishflag = false;
        case 'logxcovs',
            logxcovsflag = true;
        case 'plothist'
            plothistflag = true;
        otherwise
            error('lfp_findevtboundmax:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

offsetsamples = round(offsets/lfp_SamplePeriod);
markedsamples = cell(size(lfp_TrialIndex,1), 1);
bestxcovs = NaN(size(lfp_TrialIndex,1), 1); % normalized
EPvar = var(EP);
for trial = 1:size(lfp_TrialIndex,1)
    startevtidx = find(...
        ismember(...
        lfp_Events(lfp_TrialIndex(trial,1):lfp_TrialIndex(trial,2), 2), ...
        evtbounds{1} )) ...
        + lfp_TrialIndex(trial,1) - 1;
    if isempty(startevtidx)
        warning('lfp_findevtboundmax:evtbounds1', ...
            'Trial %d has no ''evtbounds'' start event', ...
            trial );
        continue
    else
        startevtidx = startevtidx(1);
    end
    endevtidx = find(...
        ismember(...
        lfp_Events(startevtidx:lfp_TrialIndex(trial,2), 2), evtbounds{2} )) ...
        + startevtidx - 1;
    if isempty(endevtidx)
        warning('lfp_findevtboundmax:evtbounds2', ...
            'Trial %d has no ''evtbounds'' end event', ...
            trial );
        continue
    else
        endevtidx = endevtidx(1);
    end
    startsample = lfp_time2index(lfp_Events(startevtidx,1)) ...
        + offsetsamples(1);
    endsample = lfp_time2index(lfp_Events(endevtidx,1)) ...
        + offsetsamples(2);
    [c,lags] = xcov(EP, ...
        lfp_Samples{filenum}(startsample:endsample), 'none');
    [v,ix] = max(c);
    bestxcovs(trial) = v / ( (length(EP) - 1) * sqrt( EPvar * ...
            var(lfp_Samples{filenum}(startsample:endsample)) ));
    if lags(ix) > 0
        warning('lfp_findevtboundEP:fish', ...
            'Fishy positive value for best lag on trial %d: %d samples', ...
            trial, lags(ix));
        if fishbadflag
            lfp_BadTrials = union(lfp_BadTrials, trial);
        end
    end
    if fishflag || lags(ix) <= 0
        if isempty(minxcov) || bestxcovs(trial) >= minxcov
            markedsamples{trial} = startsample - lags(ix);
        else
            warning('lfp_findevtboundEP:fowl', ...
                'Best xcov for trial %d was %.3g, skipping', ...
                trial, bestxcovs(trial));
        end
        if logxcovsflag
            lfp_log(sprintf('Trial %d: lag=%d, norm xcov=%.3g', ...
                trial, lags(ix), bestxcovs(trial) ));
        end
    end
end
ts = lfp_index2time(cell2mat(markedsamples));
if plothistflag
    figure; 
    hist(bestxcovs,100);
    xlabel('Normalized xcov');
end


