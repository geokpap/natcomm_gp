function [trialstarts, trialends] = lfp_rodentCreateTrialIndex(events)
% Requires non-empty lfp_GoodTrialEvent; value may be scalar or vector.
% Does not change any global values.

%$Rev: 318 $
%$Date: 2014-02-08 22:12:19 -0500 (Sat, 08 Feb 2014) $
%$Author: dgibson $

global lfp_GoodTrialEvent lfp_NoWaitbar lfp_NominalTrialStart ...
    lfp_NominalTrialEnd

% First, find the markers for good trials:
goodtrials = find(ismember(events(:,2), lfp_GoodTrialEvent));
if isempty(goodtrials)
    error('lfp_rodentCreateTrialIndex:nogood', ...
        'There are no good trials marked in the events file.');
end

% Then search backwards from each one to find the preceding
% lfp_NominalTrialStart and lfp_NominalTrialEnd (it is an error if there is
% no lfp_NominalTrialStart between successive goodtrialevent's), which we
% save in trialstarts and trialends respectively:
if ~lfp_NoWaitbar
    hWaitBar = waitbar(0, '', 'Name', ...
        'Locating Trial Beginnings and Endings');
end
WaitBarSteps = length(goodtrials);
CurrentWaitBarStep = 0;
trialstarts = [];
trialends = [];
for goodindex = goodtrials'
    if ~lfp_NoWaitbar
        CurrentWaitBarStep = CurrentWaitBarStep + 1;
        waitbar(CurrentWaitBarStep/WaitBarSteps, hWaitBar);
    end
    index = goodindex;
    while ~ismember(events(index,2), lfp_NominalTrialStart)
        if ismember(events(index,2), lfp_NominalTrialEnd)
            trialends = [ trialends index ]; %#ok<AGROW>
        end
        index = index-1;
        if index < 1
            if ~lfp_NoWaitbar
                close(hWaitBar);
            end
            error('lfp_rodentCreateTrialIndex:nostart', ...
                ['Could not find lfp_NominalTrialStart event' ...
                    ' for first trial.' ]);
        end
    end
    if ~isempty(trialstarts) && index == trialstarts(end)
        if ~lfp_NoWaitbar
            close(hWaitBar);
        end
        error('lfp_rodentCreateTrialIndex:missingTrialStart', ...
            ['There is no trial start event for the trial ending at ' ...
                num2str(events(goodindex,1)) 'seconds.' ]);
    end
    trialstarts = [ trialstarts index ]; %#ok<AGROW>
end
if ~lfp_NoWaitbar
    close(hWaitBar);
end

