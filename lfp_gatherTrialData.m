function [data, interval, win] = ...
    lfp_gatherTrialData(trials, filenums, win, varargin)
%[data, interval, win] = ...
%    lfp_gatherTrialData(filenums, trials, win)
% Returns a rectangular array of data, along with metadata that were used
% in constructing it.  See also lfp_getSamples.
%INPUTS
% trials:  verbatim list of trials to gather from; i.e., this value should
%   already be edited pursuant to lfp_enabledTrials as (if) desired.  It is
%   an error for it to be empty.
% filenums:  verbatim list of trials to gather from; i.e., this value
%   should already be edited pursuant to lfp_SelectedFiles as (if) desired.
%   May include NaN values, in which case the corresponding rows of <data>
%   are all NaN.
% win:  determines <interval> together with lfp_XLimAll.
%OUTPUTS
% data:  If <filenums> is a scalar, then <data> is in trials X timepoints
%   format; if <trials> is a scalar, then <data> is in filenums X
%   timepoints format; otherwise, it is in filenums X timepoints X trials
%   format.
% interval:  actual start and end points of <data> relative to the
%   reference sample
% win:  the value actually used, which is <lfp_XLimAll> if <win> is
%   initially empty.
%OPTIONS
% 'logical' - returns data as 'logical' data type.
%NOTES
%   This code was originally derived from lfp_CSCraster_gatherdata, which
% lives in lfp_CSCraster (note substitution of "window" with "win", which
% was because 'window' is actually a Matlab function.  It was originally
% meant to replace lfp_CSCraster_gatherdata, but lfp_CSCraster has several
% idiosyncratic requirements that most analysis functions don't have (the
% 'pad' option and the use of NaNs in <filenums>) that involve a lot of
% extra plumbing and testing.  I have therefore eliminated those features
% from lfp_gatherTrialData.  Such idiosyncracies should really be handled
% in wrappers to lfp_gatherTrialData.
%   lfp_gatherTrialData is considerably older (2011-04-29) than
% lfp_getSamples (2012-05-09).  lfp_getSamples was created to serve the
% needs of lfp_xcov, and lfp_gatherTrialData was created for
% lfp_burstAnalysis.  lfp_gatherTrialData offers no options other than
% 'logical', and does no default processing on <trials> or <filenums>,
% whereas lfp_getSamples offers a variety of frequently used pre-processing
% options and does the standard substitutions on <trials> and <filenums>.
% lfp_getSamples should probably be the standard function to use, and it
% should probably call lfp_gatherTrialData, but at this point (28-Aug-2013)
% there are just these two annoyingly different versions of the same thing.
%   Regression test: test_lfp_gatherTrialData.

%$Rev: 305 $
%$Date: 2013-09-05 18:24:19 -0400 (Thu, 05 Sep 2013) $
%$Author: dgibson $

global lfp_Samples lfp_XLimAll lfp_SamplePeriod

if isempty(trials)
    error('lfp_gatherTrialData:notrials', ...
        '<trials> is empty.');
end
if isempty(filenums)
    error('lfp_gatherTrialData:nofiles', ...
        '<filenums> is empty.');
end

commontimeopts = {};
logicalflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'logical'
            logicalflag = true;
        otherwise
            error('lfp_gatherTrialData:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end


reftime = []; %#ok<NASGU>

% Compute <interval>
if isempty(win)
    win = lfp_XLimAll;
end
% Standard window processing:  if win or lfp_XLimAll is not empty, use
% that, trimmed to recsegs; else use the common time of whole trials,
% trimmed to trial boundaries.
if ~isempty(win)
    commontimeopts{end+1} = 'recseg';
end
[interval, rawtrialinfo] = ...
    lfp_findCommonTime(trials, commontimeopts{:});
if ~isempty(win)
    winpts = round(win/lfp_SamplePeriod);
    interval(1) = max(winpts(1), interval(1));
    interval(2) = min(winpts(2), interval(2));
end
if any(rawtrialinfo(:,3)==0)
    error('lfp_gatherTrialData:noref', ...
        'The reference event is missing from trial(s) %s', ...
        dg_canonicalSeries(trials(rawtrialinfo(:,3)==0)) );
end
% <numpts> is the number of samples from each trial:
numpts = interval(2) - interval(1) + 1;
% Construct <data> based on <dataidx>, which we construct in trials X
% timepoints format relative to each <refsample> (see lfp_findCommonTime).
dataidx = repmat(rawtrialinfo(:,3), 1, numpts) ...
    + repmat(interval(1):interval(2), size(rawtrialinfo,1), 1);
if numel(trials) == 1
    % <data> is in filenums X timepoints
    data = reshape(lfp_Samples{filenums(1)}(dataidx), 1, []);
    for k = 2:length(filenums)
        data(k,:) = reshape(lfp_Samples{filenums(k)}(dataidx), 1, []);
    end
elseif numel(filenums) == 1
    % <data> is in trials X timepoints
    data = lfp_Samples{filenums}(dataidx);
else
    for fileidx = 1:numel(filenums)
        data(fileidx, :, :) = permute( ...
            lfp_Samples{filenums(fileidx)}(dataidx), [3 2 1]);
    end
end
if logicalflag
    data = logical(data);
end
