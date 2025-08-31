function result = dg_ComputePerTrial(measure, segspecs, trialdata, ...
    clusternum, condition)

%dg_ComputePerTrial Computes various measures per trial over specified time
%segments.
% result = dg_ComputePerTrial(measure, segspecs, trialdata, ...
%    clusternum, condition)
% Returns a matrix containing an element for each time segment specification 
% in segspecs computed for each trial in trialdata.  
%
% measure: is a string that may have one of these values, used to control
%   what measure will be computed:
%       'duration'
%           Seconds.  Simply (stop-start).  clusternum is ignored.
%       'rate'
%           Spikes per second.  
%       'spikes'
%           Raw spike count.  The criterion for including a spike is 
%           (start < timestamp <= stop).
% segspecs: should be an
%   array of structures of the type returned by dg_ParseSegmentSpec.  The
%   offsets are assumed to be in mS.
% trialdata:  should be an array of trialdata of the type returned by
%   dg_ReadMouseFormat.  result(m,n) will contain the result for
%   trialdata(m) computed over segspecs(n), or NaN if an error was
%   encountered during  computation.  
% clusternum: is the ID number of the cluster in the
%   trialdata for which result will be computed.  
% condition: is a string which,
%   when evaluated by eval in the context of the innermost loop, must yield
%   a true result, or else the result will be reported as NaN and a warning
%   raised.  (See dg_ParseSelectionSpec, which returns such a string.  Note
%   that this whole mechanism is kludgy and should be replaced by a
%   recursive data structure and matching interpreter.)
% Note that event TTL ID 1 is treated as special, because it is the "Disc
% On" ("Start Recording") event, and therefore always is implicitly present
% and has timestamp 0, whereas for other events a 0 timestamp means that
% the event did not occur.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

result = NaN + ones(length(trialdata), length(segspecs));
for segspecnum = (1:length(segspecs))
    for trialnum = (1:length(trialdata))
        try
            trialselected = eval(condition);
        catch
            error('dg_ComputePerTrial:BadSelection', ...
                'Could not evaluate the selection condition on trial %d', ...
                trialnum);
        end
        if trialselected
            try
                result(trialnum, segspecnum) = compute(measure, ...
                    segspecs(segspecnum), trialdata(trialnum), clusternum);
            catch
                err = lasterror;
                if strncmp(err.identifier, 'dg_ComputePerTrial:', 15)
                    result(trialnum, segspecnum) = NaN;
                    warning(err.identifier, '%s on seg spec %d trial %d', ...
                        err.message, segspecnum, trialnum);
                else
                    rethrow(err);
                end
            end
        else
            result(trialnum, segspecnum) = NaN;
            warning('Trial %d rejected by selection condition', trialnum);
        end
    end
end


function result = compute(measure, segspec, trial, clusternum)
% Returns one measure computed for unit clusternum during trial over the 
% time period specified by segspec.  
% measure:  same as for dg_ComputePerTrial.  
% segspec:  same as for dg_ComputePerTrial.  
% trial:  one element of trialdata.  
% clusternum:  same as for dg_ComputePerTrial.  
% Returns empty array ([]) if error.  Start or stop times that work out to
% less than zero are replaced by zero.
result = [];
for id = segspec.start.idlist'
    if trial.events(id) || id == 1
        break;
    end
end
if trial.events(id) == 0 && id ~= 1
    error('dg_ComputePerTrial:noStart', 'no start event found');
end
startstamp = double(trial.events(id)) + 10 * segspec.start.offset;
startstamp = max([startstamp 0]);
for id = segspec.stop.idlist'
    if trial.events(id) || id == 1
        break;
    end
end
if trial.events(id) == 0 && id ~= 1
    error('dg_ComputePerTrial:noStop', 'no stop event found');
end
stopstamp = double(trial.events(id)) + 10 * segspec.stop.offset;
stopstamp = max([stopstamp 0]);
if (stopstamp <= startstamp)
    error('dg_ComputePerTrial:nonPosDur', ...
        'nonpositive segment duration: %d to %d', startstamp, stopstamp);
else
    if strcmp(measure, 'rate') || strcmp(measure, 'spikes')
        spikes = trial.spikes(1:trial.header.SSize(clusternum), clusternum);
        spikes_included = spikes > startstamp & spikes <= stopstamp;
    end
    switch measure
        case 'duration'
            result = (stopstamp-startstamp)/10000;
        case 'rate'
            result = 10000 * sum(spikes_included)/(stopstamp-startstamp);
        case 'spikes'
            result = sum(spikes_included);
        otherwise
            error('dg_ComputePerTrial:BadMeasure', ...
                'Bad value for measure: %s.', measure);
    end
end