function trialnum = lfp_getTrialNum(trialid, session)
%LFP_GETTRIALNUM converts from string trial ID to numeric internal trial
%number.
%trialnum = lfp_getTrialNum(trialid)
%trialnum = lfp_getTrialNum(trialid, session)
%   <trialid> is a string of the form returned by lfp_getTrialID.
%   <trialnum> is the integer index into lfp_SelectedTrials or row number
%   in lfp_TrialIndex.  The session name is normally parsed out from
%   <trialid>.  <session>, if given, is a default session name to use in
%   case <trialid> contains only the original trial number
%   (lfp_OrigTrialNums).  Else, if there is only one session loaded, then
%   that's the session.  Else, error.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 2
    session = '';
end

delimiter = strfind(trialid, '-');
if isempty(delimiter)
    sessionname = session;
    trialstring = trialid;
else
    sessionname = trialid(1:delimiter(end)-1);
    trialstring = trialid(delimiter(end)+1:end);
end
if isempty(sessionname) 
    if (length(lfp_SessionNames) == 1)
        sessionname = lfp_SessionNames{1};
    else
        error('lfp_getTrialNum:nosession', ...
            'Cannot identify session for Trial ID %s', trialid );
    end
end

sessionnum = find(ismember(lfp_SessionNames, sessionname));
if isempty(sessionnum)
    error('lfp_getTrialNum:nosuchsession', ...
        'The session ''%s'' is not loaded', sessionname );
end
firsttrialnum = lfp_SessionFirstTrials(sessionnum);
if sessionnum < length(lfp_SessionFirstTrials)
    lasttrialnum = lfp_SessionFirstTrials(sessionnum + 1) - 1;
else
    lasttrialnum = length(lfp_OrigTrialNums);
end

orgtrialnum = str2num(trialstring);
if isempty(orgtrialnum)
    error('lfp_getTrialNum:badtrialnum', ...
        'Trial number is not numeric: %s', trialstring );
end
relativetrialnum = find(lfp_OrigTrialNums(firsttrialnum:lasttrialnum) == orgtrialnum);
if isempty(relativetrialnum)
    error('lfp_getTrialNum:badtrialnum2', ...
        'Trial number is not loaded: %s', trialid );
end
trialnum = firsttrialnum - 1 + relativetrialnum;

