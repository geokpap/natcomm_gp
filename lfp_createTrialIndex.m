function trialindex = lfp_createTrialIndex(events, numTP)
% Create trial structure: return an index whose columns are:
% start event index, end event index.  Also sets lfp_SelectedTrials to same
% length as lfp_TrialIndex, all trials selected.
% <numTP> is optional, and specifies the number of sets of trial parameters
% that is used for the error check at the end; if not given,
% numel(lfp_TrialParams) is used.

%$Rev: 414 $
%$Date: 2022-01-14 21:08:58 -0500 (Fri, 14 Jan 2022) $
%$Author: dgibson $

lfp_declareGlobals;
if nargin < 2
    numTP = numel(lfp_TrialParams);
end

% First, find good trials and starts and ends of good trials:
switch lfp_SetupType
    case {'monkey' 'naotaka'}
        [trialstartindex, trialendindex] = ...
            lfp_monkeyCreateTrialIndex(events);
    case 'rodent'
    % case 'rodent', but not when using kludged tone-noise test data
        [trialstartindex, trialendindex] = ...
            lfp_rodentCreateTrialIndex(events);
    otherwise
%     case {'theresa' 'wheel'}
        trialstartindex = find(ismember(events(:,2), lfp_NominalTrialStart))';
        trialendindex = find(ismember(events(:,2), lfp_NominalTrialEnd))';
end

[pairs, extraStarts, extraEnds] = ...
    dg_zip(trialstartindex, trialendindex);
% A single trailing extra start should be ignored; delete it:
if ~isempty(extraStarts) && (extraStarts(end) == length(trialstartindex))
    trialstartindex(extraStarts(end)) = [];
    extraStarts(end) = [];
end
% A single leading extra end should be ignored; delete it:
if ~isempty(extraEnds) && (extraEnds(1) == 1)
    trialendindex(1) = [];
    extraEnds(1) = [];
    % adjust indices into trialendindex because we deleted the first
    % element:
    pairs(:,2) = pairs(:,2) - 1;
    extraEnds = extraEnds - 1;
end
if ~isempty(extraStarts)
    warning('lfp_createTrialIndex:extraStarts', ...
        'Skipping extra trial start event(s) at timestamp(s):\n%s', ...
        dg_thing2str(events(trialstartindex(extraStarts), 1)) );
    s = warning('query', 'lfp_createTrialIndex:extraStarts');
    if isequal(s.state, 'on')
        lfp_log(sprintf(...
            'Skipping extra trial start event(s) at timestamp(s):\n%s', ...
            mat2str(events(trialstartindex(extraStarts), 1)) ));
    end
    trialstartindex(extraStarts) = [];
end
if ~isempty(extraEnds)
    warning('lfp_createTrialIndex:extraEnds', ...
        'Skipping extra trial end event(s) at timestamp(s):\n%s', ...
        dg_thing2str(events(trialendindex(extraEnds), 1)) );
    s = warning('query', 'lfp_createTrialIndex:extraEnds');
    if isequal(s.state, 'on')
        lfp_log(sprintf(...
            'Skipping extra trial end event(s) at timestamp(s):\n%s', ...
            mat2str(events(trialendindex(extraEnds), 1)) ));
    end
    trialendindex(extraEnds) = [];
end

% Finally, construct lfp_TrialIndex:
if ( isempty(trialstartindex) ...
        || isempty(trialendindex) ...
        )
    error('lfp_createTrialIndex:emptycol', ...
        'One or more columns of lfp_TrialIndex is empty.');
end
trialindex = [ trialstartindex' trialendindex' ];
if numTP && ...
        (size(trialindex,1) ~= numTP)
    error('lfp_createTrialIndex:TPmismatch', ...
        'length of lfp_TrialParams = %d, number of trials = %d.', ...
        numTP, size(trialindex,1) );
end


