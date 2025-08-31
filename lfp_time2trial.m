function [trialnum, trialID, percent] = lfp_time2trial(ts, varargin)
%[trialnum, trialID, percent] = lfp_time2trial(ts)
%   Finds the trialnum, the Unique Trial ID, and the percent of trial
%   duration at which a given timestamp <ts> occurs.  The trial that <ts>
%   belongs to is defined as the trial having the greatest start time that
%   is less than or equal to <ts>.  <percent> represents the percentage of
%   the trial that has elapsed at <ts>.  If <percent> is greater than 100,
%   then <ts> occurs during the intertrial interval following the trial. If
%   it is much greater than 100, then there is most likely a gap in the
%   recording between the end of the trial and <ts>.  <ts> can be an array,
%   in which case trialnum, trialID, and percent are column vectors.
%   Returns all zeros if <ts> does not belong to any trial.
%OPTIONS
% 'evtbounds', {startIDs stopIDs} - works as for lfp_disp to return the
%   trialnum for each element of <ts> that falls within the event bounds
%   and returns negative trialnum for each element that falls outside the
%   event bounds. A warning is raised for each trial that contains an
%   element of <ts> but does not contain both bounding events, and all <ts>
%   belonging to those trials are returned negative.

%$Rev: 178 $
%$Date: 2010-12-08 21:13:54 -0500 (Wed, 08 Dec 2010) $
%$Author: dgibson $

global lfp_TrialIndex lfp_Events;

evtbounds = {};
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'evtbounds'
            argnum = argnum + 1;
            evtbounds = varargin{argnum};
        otherwise
            error('lfp_time2trial:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
trialnum = NaN(length(ts), 1);
if nargout > 1
    trialID = cell(length(ts), 1);
    if nargout > 2
        percent = NaN(length(ts), 1);
    end
end
for tsidx = 1:numel(ts)
    prevtrials = find(lfp_Events(lfp_TrialIndex(:,1),1) <= ts(tsidx));
    if isempty(prevtrials)
        trialnum(tsidx) = 0;
        trialID{tsidx} = 0;
        percent(tsidx)= 0;
    else
        trialnum(tsidx) = prevtrials(end);
        if nargout > 1
            trialID{tsidx} = lfp_getTrialID(trialnum(tsidx));
            if nargout > 2
                percent(tsidx) = round(100 * (ts(tsidx) - ...
                    lfp_Events(lfp_TrialIndex(trialnum(tsidx),1), 1) ) ...
                    / (lfp_Events(lfp_TrialIndex(trialnum(tsidx),2), 1) ...
                    - lfp_Events(lfp_TrialIndex(trialnum(tsidx),1), 1) ));
            end
        end
        if ~isempty(evtbounds)
            if tsidx == 1 || trialnum(tsidx-1) ~= trialnum(tsidx)
                % find event bounds
                trial = trialnum(tsidx);
                startevtidx = lfp_TrialIndex(trial,1);
                endevtidx = lfp_TrialIndex(trial,2);
                if ~isempty(evtbounds)
                    startevtidx = find(...
                        ismember(...
                        lfp_Events(startevtidx:endevtidx, 2), ...
                        evtbounds{1} )) + startevtidx - 1;
                    if isempty(startevtidx)
                        warning('lfp_time2trial:evtbounds1', ...
                            'Trial %d has no ''evtbounds'' start event', ...
                            trial );
                        startevtidx = NaN;
                    else
                        startevtidx = startevtidx(1);
                    end
                    if ~isnan(startevtidx)
                        endevtidx = find(...
                            ismember(...
                            lfp_Events(startevtidx:endevtidx, 2), ...
                            evtbounds{2} )) + startevtidx - 1;
                        if isempty(endevtidx)
                            warning('lfp_time2trial:evtbounds2', ...
                                'Trial %d has no ''evtbounds'' end event', ...
                                trial );
                            endevtidx = NaN;
                        else
                            endevtidx = endevtidx(1);
                        end
                    end
                end
            end
            if isnan(startevtidx) || isnan(endevtidx) || ...
                    ts(tsidx) < lfp_Events(startevtidx, 1) || ...
                    ts(tsidx) > lfp_Events(endevtidx, 1)
                trialnum(tsidx) = -trialnum(tsidx);
            end
        end
    end
end
