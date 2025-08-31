function [durations, sortidx] = lfp_selectByDuration(varargin)
%[durations, sortidx] = lfp_selectByDuration
%   Selects trials based on their duration, so that the total time that
%   would go into an average (which is the product of the number of
%   selected trials with the longest common recorded time interval within
%   lfp_XLimAll) is maximized.  It is assumed that lfp_XLimAll contains only
%   trial ends or trial beginnings, not some of each (an error is raised if
%   this assumption is false).
% OUTPUTS:
%   <durations> is the sorted list of trial durations; <sortidx> is the
%   index into the original list of trials before sorting (which is 
%   the same as the trialnums unless the 'and' option is chosen).
% OPTIONS:
%   'and' - restricts selection to those trials that are already selected.
%   'noexpand' - does not allow lfp_XLimAll to expand trials beyond their
%   bounds, only to limit them to shorter intervals, so we get the longest
%   common within-trial interval within lfp_XLimAll.
%   'fullwindow' - select only trials whose duration (as returned by
%   lfp_getTrialTimes) is at least the specified lfp_XLimAll duration minus
%   two sample periods (to allow for rounding error at both ends). Use of
%   this option bypasses the normal total-time maximization calculation.
%
% SIDE EFFECTS:
%   Sets lfp_SelectionRule
%   Logs new value of lfp_SelectionRule

%$Rev: 406 $
%$Date: 2019-11-22 16:10:40 -0500 (Fri, 22 Nov 2019) $
%$Author: dgibson $

global lfp_XLimAll lfp_SelectedTrials lfp_SelectionRule lfp_SamplePeriod ...
    lfp_BadTrials

if isempty(lfp_XLimAll)
    error('lfp_XLimAll is empty.');
end

argnum = 1;
andflag = false;
fullwindowflag = false;
noexpandflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'and'
            andflag = true;
        case 'fullwindow'
            fullwindowflag = true;
        case 'noexpand'
            noexpandflag = true;
        otherwise
            error('lfp_selectByDuration:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if andflag
    trials = lfp_enabledTrials(find(lfp_SelectedTrials));
else
    trials = 1:length(lfp_SelectedTrials);
    trials = setdiff(trials, lfp_BadTrials);
end

% Find endpoints of each trial subject to lfp_XLimAll and lfp_TrialRec, and
% corresponding durations:
trialinfo = zeros(length(trials), 2);
if noexpandflag
    window = [];
else
    window = lfp_XLimAll;
end
for trialidx = 1:length(trials)
    trial = trials(trialidx);
    [reltrialstart, reltrialend] = lfp_getTrialTimes(trial, window);
    trialinfo(trialidx, :) = [reltrialstart, reltrialend];
end
startsAreShort = any(trialinfo(:,1) > lfp_XLimAll(1));
endsAreShort = any(trialinfo(:,2) < lfp_XLimAll(2));
if (startsAreShort && endsAreShort)
    error('lfp_XLimAll contains trial beginnings and trial endings.');
end
durations = trialinfo(:,2) - trialinfo(:,1);

if fullwindowflag
    istooshort = durations < diff(lfp_XLimAll) - ...
        2 * lfp_SamplePeriod;
    if andflag
        lfp_SelectedTrials(trials(istooshort)) = false;
    else
        lfp_SelectedTrials(:) = false;
        lfp_SelectedTrials(trials) = ~istooshort;
    end
else
    % Compute total duration common to all trials, all but the shortest
    % trial, all but the shortest two trials, etc.  The total common
    % duration is equal to the number of trials times the duration of the
    % shortest trial that is included.
    [durations, sortidx] = sort(durations);
    totaldur = (length(durations) : -1 : 1)' .* durations;
    % Find the maximum total duration and select the trials that were
    % included in it:
    [maxdur, maxduridx] = max(totaldur); %#ok<ASGLU>
    lfp_SelectedTrials(:) = false;
    lfp_SelectedTrials(trials(sortidx(maxduridx:end))) = true;
end

% Set and log lfp_SelectionRule:
argstr = '';
if ~isempty(varargin)
    argstr = dg_thing2str(varargin{1});
end
for argnum = 2:length(varargin)
    argstr = sprintf('%s, %s', argstr, dg_thing2str(varargin{argnum}));
end
if andflag
    lfp_SelectionRule = sprintf('lfp_selectByDuration(%s) from %s', ...
        argstr, lfp_SelectionRule);
else
    lfp_SelectionRule = sprintf('lfp_selectByDuration(%s)', argstr);
end
lfp_log(sprintf('Applied selection: %s', lfp_SelectionRule));
