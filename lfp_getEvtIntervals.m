function [result, trialIDs] = ...
    lfp_getEvtIntervals(trials, segspecstrings, varargin)
%[result, trialIDs] = lfp_getEvtIntervals(trials, segspecstrings)
% Computes inter-event result much like PerTrialStats.  The events that
% define each interval must both belong to the same trial.  See lfp_IEI to
% compute result that may cross trial boundaries.

%INPUTS
% <segspecstrings> is a cell array of strings suitable for passing to
% dg_ParseSegmentSpec, that is, of the form (for example):
%   '31,38 +200 to 18 -15'
% which specifies the time interval starting 200 ms after the first
% instance of either of events 31 or 38 and ending 15 ms before the next
% instance of event 18. If left empty or omitted, then the duration from
% lfp_NominalTrialStart to lfp_NominalTrialEnd is used.  Substitutes values
% of globals and variables that are defined in lfp_getEvtIDs, except for
% any by the following names:
%   segspecs, specidx, specstring, token, tosubstitute, trials
% Also, note that the following variable names must NOT be set by
% lfp_getEvtIDs or else they will interfere with correct operation of this
% function:
%   segspecstrings, trials
%   
% <trials> is as for lfp_disp.
%
%OUTPUTS
%<trialIDs> is a cell column vector containing one element for each
% included trial; it may be shorter in length than <trials> if
% lfp_BadTrials is not empty. 
%<result> contains the inter-event intervals, with one column for each
% element of <segspecstrings> and one row for each trial, i.e. it is
% length(trialIDs) x length(segspecstrings). 
%
%OPTIONS
% 'rate', clustnum - instead of returning inter-event result, returns
%   the firing rate for spike channel <clustnum> computed over each
%   inter-event interval.  To get average over trials, compute the average
%   of each column.

%$Rev: 90 $
%$Date: 2009-10-30 19:11:42 -0400 (Fri, 30 Oct 2009) $
%$Author: dgibson $

lfp_declareGlobals;
clustnum = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'rate'
            argnum = argnum + 1;
            clustnum = varargin{argnum};
        otherwise
            error('lfp_getEvtIntervals:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
if nargin < 2 || isempty(segspecstrings)
    startstr = sprintf('%d,', lfp_NominalTrialStart);
    startstr(end) = [];
    endstr = sprintf('%d,', lfp_NominalTrialEnd);
    endstr(end) = [];
    segspecstrings = sprintf('%s to %s', startstr, endstr);
end
if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
lfp_getEvtIDs;
if strcmp(class(trials), 'char')
    trials = lfp_parseTrialStr(trials, session);
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);


for specidx = 1:numel(segspecstrings)
    tosubstitute = segspecstrings{specidx};
    specstring = '';
    while ~isempty(tosubstitute)
        [token, tosubstitute] = strtok(tosubstitute);
        if exist(token) == 1 && isequal(class(eval(token)), 'double')
            specstring = [ specstring sprintf('%d,', eval(token)) ];
            specstring(end) = [];
        else
            specstring = [ specstring ' ' token ' ' ];
        end
    end
    segspecs(specidx) = dg_ParseSegmentSpec(specstring);
end

trialIDs = cell(length(trials),1);
result = NaN(length(trialIDs), length(segspecstrings));
for trialidx = 1:length(trials)
    trialnum = trials(trialidx);
    trialIDs{trialidx} = lfp_getTrialID(trialnum);
    startevtidx = lfp_TrialIndex(trialnum,1);
    endevtidx = lfp_TrialIndex(trialnum,2);
    trialevents = lfp_Events(startevtidx : endevtidx, :);
    for specidx = 1:length(segspecs)
        startidx = find(ismember( ...
            trialevents(:,2), segspecs(specidx).start.idlist ));
        if ~isempty(startidx)
            starttime = trialevents(startidx(1), 1) ...
                + segspecs(specidx).start.offset * 1e-3;
            endidx = find(ismember( ...
                trialevents(startidx(1):end,2), segspecs(specidx).stop.idlist )) ...
                + startidx(1) - 1;
            if ~isempty(endidx)
                endtime = trialevents(endidx(1), 1) ...
                    + segspecs(specidx).stop.offset * 1e-3;
                if isempty(clustnum)
                    result(trialidx, specidx) = endtime - starttime;
                else
                    numspikes = sum( lfp_Spikes{clustnum} > starttime ...
                        & lfp_Spikes{clustnum} <= endtime );
                    result(trialidx, specidx) = numspikes / ...
                        (endtime - starttime);
                end
            end
        end
    end
end
