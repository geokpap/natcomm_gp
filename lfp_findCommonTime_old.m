function [interval, rawtrialinfo] = lfp_findCommonTime(trials, varargin)
%[interval, rawtrialinfo] = lfp_findCommonTime(trials)
%   Returns the time interval shared by all trialnums in <trials> expressed
%   in samples relative to lfp_AlignmentRef, plus the <rawtrialinfo> table on
%   which <interval> is based, containing [startsample endsample refsample]
%   for each trialnum in <trials>.  If there is no lfp_AlignmentRef, then
%   refpoint for that trial is reported as 0, and the trial is excluded
%   from the calculation of <interval>.  <interval> is expressed in units
%   of sample counts, with negative numbers indicating samples before the
%   reference event.
%lfp_findCommonTime(..., 'recseg')
%   The bounds are the enclosing recorded segments rather than start and
%   end of trial.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
argnum = 1;
recsegflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'recseg'
            recsegflag = true;
        otherwise
            error('lfp_findCommonTime:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

rawtrialinfo = [];

for trial = trials
    trialevents = lfp_Events( ...
        lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), ...
        : );
    reftime = trialevents( ...
        find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
        1 );
    if length(reftime) == 0
        refpoint = 0;
    else
        refpoint = lfp_time2index(reftime(1));
    end
    if recsegflag
        rawtrialinfo(end+1,:) = ...
            [lfp_TrialRec(trial,1) lfp_TrialRec(trial,2) refpoint];
    else
        rawtrialinfo(end+1,:) = ...
            [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4) refpoint];
    end
end

trialinfo = rawtrialinfo(find(rawtrialinfo(:,3) ~= 0), :);
pointsbefore = min(trialinfo(:,3) - trialinfo(:,1));
pointsafter = min(trialinfo(:,2) - trialinfo(:,3));
interval = [ -pointsbefore pointsafter];

