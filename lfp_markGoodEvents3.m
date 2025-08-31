%LFP_MARKGOODEVENTS3
% Script for use by lfp_parseEyeTraces3, NOT coded as a function in
% order to save the time and space used by Matlab's compulsion to copy
% parameter values on function entry/exit.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

% Presume events to be good until proven bad:
eventgood = repmat(true, 1, length(rawevents));
% preceding Blink End index:
pBlinkIdx = dg_binsearch(blinkends, rawevents(1)) - 1;
% next Blink Start index:
nBlinkIdx = dg_binsearch(blinkstarts, rawevents(1));
% next Ambiguous event index:
nFailedIdx = dg_binsearch(failedblinks, rawevents(1));
% preceding Ambiguous event index:
if ~isempty(nFailedIdx)
    pFailedIdx = nFailedIdx - 1;
else
    pFailedIdx = [];
end
recidx = find(lfp_RecSegments(:,1) < rawevents(1));
recidx = recidx(end);
for idx = 1:length(rawevents)
    % keep the bracketing blinks and RecSegment current, taking advantage
    % of the fact that all the lists are ordered:
    if (rawevents(idx) > lfp_RecSegments(recidx, 2))
        recidx = recidx + 1;
    end
    if (pBlinkIdx < length(blinkends)) && ...
            (rawevents(idx) > blinkends(pBlinkIdx + 1))
        % We need to stop incrementing at the end of the blink lists, which
        % can have 2 configurations: (1) All integral blinks, nBlinkIdx =
        % end + 1 (i.e. do not test). (2) End w/ blinkstart
        % (blinkstarts(end) > blinkends(end)), nBlinkIdx = end (we soon run
        % out of rawevents in real eye traces). BUT the inclusion of Rec On
        % / Offs implies that there can be any number of unmatched
        % blinkstarts and blinkends, so we can't actually increment but
        % must instead find the proper events each time.
        if nBlinkIdx == length(blinkstarts) ...
                && blinkstarts(end) > blinkends(end)
            % do nothing
        elseif nBlinkIdx > length(blinkstarts)
            % also do nothing
        else
            pBlinkIdx = dg_binsearch(blinkends, rawevents(idx)) - 1;
            nBlinkIdx = dg_binsearch(blinkstarts, rawevents(idx));
        end
    end
    % Ditto for ambiguous events:
    if (~isempty(pFailedIdx) && pFailedIdx < length(failedblinks)) && ...
            (rawevents(idx) > failedblinks(pFailedIdx + 1))
        if ~isempty(nFailedIdx) && nFailedIdx > length(failedblinks)
            % also do nothing
        else
            nFailedIdx = dg_binsearch(failedblinks, rawevents(idx));
            if ~isempty(nFailedIdx)
                pFailedIdx = nFailedIdx - 1;
            else
                pFailedIdx = [];
            end
        end
    end
    % event times must be at least NBlinkTO after the preceding blink
    % end or ambiguous sacc/blink, and NRecTO after the preceding
    % Record On:
    if (pBlinkIdx > 0) && ...
            (rawevents(idx) < blinkends(pBlinkIdx) + NBlinkTO)
        eventgood(idx) = false;
    end
    if ~isempty(pFailedIdx) && (pFailedIdx > 0) && ...
            (rawevents(idx) < failedblinks(pFailedIdx) + NBlinkTO)
        eventgood(idx) = false;
    end
    if (rawevents(idx) < lfp_RecSegments(recidx, 1) + NRecTO)
        eventgood(idx) = false;
    end
    % ... and also at least NBlinkTO before the next blink start
    %  or ambiguous sacc/blink, and NRecTO before the next Rec Off:
    if (nBlinkIdx <= length(blinkstarts)) && ...
            (rawevents(idx) > blinkstarts(nBlinkIdx) - NBlinkTO)
        eventgood(idx) = false;
    end
    if ~isempty(nFailedIdx) && (nFailedIdx <= length(failedblinks)) && ...
            (rawevents(idx) > failedblinks(nFailedIdx) - NBlinkTO)
        eventgood(idx) = false;
    end
    if (rawevents(idx) > lfp_RecSegments(recidx, 2) - NRecTO)
        eventgood(idx) = false;
    end
end
