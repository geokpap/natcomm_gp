function lfp_checkRecSegments
% Logs and warns about any gap in recording that occurs during a trial,
% which by definition means any gap which has either endpoint
% during a trial.
% Only applies to enabled trials, so putting a trial on
% lfp_BadTrials or de-selecting it will stop the reporting for that trial.
% Three different warnings may be issued, which denote:
%   lfp_checkRecSegments:badrec1 - a gap starts more than one frame before
%   the end of the trial
%   lfp_checkRecSegments:badrec2 - a gap starts within one frame before
%   the end of the trial
%   a gap starts more than one frame before
%   lfp_checkRecSegments:badrec3 - a gap ends after the start of the trial
% In the case of a gap that is contained entirely within the trial, only
% badrec1 or badrec2 will be issued.

%$Rev: 61 $
%$Date: 2009-04-28 17:40:05 -0400 (Tue, 28 Apr 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if isempty(lfp_ActiveFilenums)
    error('lfp_checkRecSegments:noCSC', ...
        'There must be at least one CSC file loaded.' );
end

difs = diff(lfp_TimeStamps);
framedur = lfp_SamplesPerFrame * lfp_SamplePeriod;
gaps = find((difs < framedur/2) | (difs > framedur*1.5));

for trial = lfp_enabledTrials(1:size(lfp_TrialIndex,1))
    % endsampTS is timestamp of last sample recorded before gap.
    % firstsampTS is timestamp of first sample recorded after gap.
    endsampTS = ...
        lfp_TimeStamps(gaps) + (lfp_SamplesPerFrame - 1) * lfp_SamplePeriod;
    firstsampTS = lfp_TimeStamps(gaps + 1);
    % It's OK if the end-of-trial sample is endsampTS, since Rec Off is
    % used as lfp_NominalTrialEnd, but it's bad if the start-of-trial
    % sample is firstsampTS:
    isbadgapstart = endsampTS >= lfp_Events(lfp_TrialIndex(trial,1)) ...
        & endsampTS < lfp_Events(lfp_TrialIndex(trial,2));
    isbadgapend = firstsampTS >= lfp_Events(lfp_TrialIndex(trial,1)) ...
        & firstsampTS <= lfp_Events(lfp_TrialIndex(trial,2));
    isbadgap =  isbadgapstart | isbadgapend;
    if any(isbadgap)
        msg = sprintf('Trial %d contains recording gap(s):', trial);
        for gapidx = reshape(find(isbadgap), 1, [])
            msg = sprintf('%s\n%.6f - %.6f', msg, ...
                endsampTS(gapidx), firstsampTS(gapidx));
            if isbadgapstart(gapidx)
                delay = endsampTS(gapidx) - ...
                    lfp_Events(lfp_TrialIndex(trial,2));
                msg = sprintf('%s (%.6f re trial end)', msg, delay);
                if delay >= framedur
                    warning('lfp_checkRecSegments:badrec1', msg);
                else
                    warning('lfp_checkRecSegments:badrec2', msg);
                end
            else
                msg = sprintf('%s (%.6f re trial start)', msg, ...
                    firstsampTS(gapidx) - ...
                    lfp_Events(lfp_TrialIndex(trial,1)));
                warning('lfp_checkRecSegments:badrec3', msg);
            end
        end
        lfp_log(msg);
    end
end
