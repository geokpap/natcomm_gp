function IEI = lfp_IEI(evtid1, evtid2, varargin)
%IEI = lfp_IEIhist(evtid1, evtid2)
% For the first event with <evtid1> in a currently enabled trial,
% compute the time interval to the immediately following <evtid2>,
% regardless of whether that <evtid2> belongs to a trial or not.  Only the
% first instance of <evtid1> within a single trial is used, and a warning
% is issued if there are multiple instances.  An error is raised if any
% selected trial does not contain an <evtid1>.  <IEI> is the vector of
% inter-event intervals.  If no <evtid2> is found following a particular
% <evtid1>, a warning is issued and the interval is considered to be zero
% (note that all zero values should therefore be at the end of <IEI>).
%INPUTS
% evtid1: scalar or vector of event IDs.  Vectors are handled as for
%   <lfp_AlignmentRef>, i.e. using 'ismember'.
% evtid2: scalar event ID.
%OPTIONS
% 'last' - when multiple instances of <evtid1> are found within a single
%   trial, use the last one and do not raise a warning.

%$Rev: 414 $
%$Date: 2022-01-14 21:08:58 -0500 (Fri, 14 Jan 2022) $
%$Author: dgibson $

lfp_declareGlobals;

lastflag = false;

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
        case 'last'
            lastflag = true;
        otherwise
            error('lfp_IEI:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
IEI = zeros(length(trials), 1);
for trialidx = 1:length(trials)
    trial = trials(trialidx);
    idx1 = find(ismember( ...
        lfp_Events(lfp_TrialIndex(trial,1):lfp_TrialIndex(trial,2), 2), ...
        evtid1 )) + lfp_TrialIndex(trial,1) - 1;
    if isempty(idx1)
        error('lfp_IEI:noref', ...
            'There is no event %d in trial %d)', evtid1, trial );
    else
        if length(idx1) > 1
            if lastflag
                idx1 = idx1(end);
            else
                warning('lfp_IEI:multiref', ...
                    'There is more than one event %d in trial %d)', ...
                    evtid1, trial );
                idx1 = idx1(1);
            end
        end
    end
    idx2 = idx1 + 1;
    while idx2 <= size(lfp_Events,1)
        if lfp_Events(idx2, 2) == evtid2
            IEI(trialidx) = lfp_Events(idx2, 1) - lfp_Events(idx1, 1);
            break
        else
            idx2 = idx2 + 1;
        end
    end
    if idx2 > size(lfp_Events,1)
        warning('lfp_IEI:noevt2', ...
            'Failed to find evtid2 following evtid1 at timestamp %.6f', ...
            lfp_Events(idx1, 1) );
    end
end

