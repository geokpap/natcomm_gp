function ts = lfp_findevtboundmax(filenum, evtbounds, offsets, varargin)
%ts = lfp_findevtboundmax(filenum); for lfp_createEvents
% Finds the absolute maximum-valued sample in lfp_Samples{filenum} in the
% interval defined by <evtbounds> for each trial.  Duplicate values are
% handled as for the Matlab function 'max'. As for the lfp_disp option
% 'evtbounds', <evtbounds> is of the form {startIDs endIDs}, where
% <startIDs> is a list of alternative event IDs to use as the start of the
% time range, and <endIDs> is a list of alternative event IDs to use as the
% end of the time range.  The chronologically first event that is a member
% of the list is used, subject to the constraint that the end event must be
% after the start event.  Trials that are lacking a start or end event are
% ignored, but raise a warning.  Trials in which the maximum value is at
% the start or end event raise a warning and are re-processed after
% replacing the offending endpoint with the closest minimum or pair of
% equal values that is within the original <evtbounds> interval.  If there
% is no such replacement endpoint, then the trial is ignored.
%OPTIONS
% 'invert' - finds the minimum instead of the maximum.

%$Rev: 126 $
%$Date: 2010-06-01 19:14:25 -0400 (Tue, 01 Jun 2010) $
%$Author: dgibson $

lfp_declareGlobals;

invertflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'invert'
            invertflag = true;
        otherwise
            error('lfp_findevtboundmax:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

offsetsamples = round(offsets/lfp_SamplePeriod);
markedsamples = cell(size(lfp_TrialIndex,1), 1);
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
    if invertflag
        mywave = -lfp_Samples{filenum}(startsample:endsample);
    else
        mywave = lfp_Samples{filenum}(startsample:endsample);
    end
    [v, markedsamples{trial,1}] = max(mywave);
    markedsamples{trial,1} = markedsamples{trial,1} + startsample - 1;
    redo = false;
    if any(ismember(markedsamples{trial,1}, startsample))
        warning('lfp_findevtboundmax:maxstart', ...
                'Trial %d has max at start', trial);
        num2trim = 1;
        while num2trim + 1 < length(mywave) ...
                && mywave(num2trim) > mywave(num2trim + 1)
            num2trim = num2trim + 1;
        end
        if num2trim == length(mywave)
            warning('lfp_findevtboundmax:mono1', ...
                'Trial %d is monotonic decreasing', trial);
            markedsamples{trial,1} = [];
        else
            startsample = startsample + num2trim;
            redo = true;
        end
    end
    if any(ismember(markedsamples{trial,1}, endsample))
        warning('lfp_findevtboundmax:maxend', ...
            'Trial %d has max at end', trial);
        num2trim = 1;
        while num2trim + 1 < length(mywave) ...
                && mywave(end - num2trim) > mywave(end - num2trim - 1)
            num2trim = num2trim + 1;
        end
        if num2trim == length(mywave)
            warning('lfp_findevtboundmax:mono2', ...
                'Trial %d is monotonic increasing', trial);
            markedsamples{trial,1} = [];
        else
            endsample = endsample - num2trim;
            redo = true;
        end
    end
    if redo
        if invertflag
            mywave = -lfp_Samples{filenum}(startsample:endsample);
        else
            mywave = lfp_Samples{filenum}(startsample:endsample);
        end
        [v, markedsamples{trial,1}] = max(mywave);
        markedsamples{trial,1} = markedsamples{trial,1} + startsample - 1;
    end
end
ts = lfp_index2time(cell2mat(markedsamples));

