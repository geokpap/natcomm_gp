function eventgood = lfp_markGoodEvents4Dummies(rawevents, ...
    blinkstarts, blinkends, NBlinkTO, NRecTO)
% Totally stupid and simple function to replace the fast & clever &
% dysfunctional lfp_markGoodEvents4 script.  If there are any unmatched
% blinkstarts or blinkends, complains and ignores them.

%$Rev: 43 $
%$Date: 2009-01-12 17:59:34 -0500 (Mon, 12 Jan 2009) $
%$Author: dgibson $

lfp_declareGlobals;
eventgood = repmat(true, 1, length(rawevents));
[pairs, extraL, extraR] = dg_zip(blinkstarts, blinkends);
if ~isempty(extraL) || ~isempty(extraR)
    msg = sprintf( ...
        'There are %d unpaired blink starts and %d unpaired blink ends', ...
        length(extraL), length(extraR) );
    warning('lfp_markGoodEvents4Dummies:unpaired', ...
         '%s', msg);
     lfp_log(msg);
    blinkstarts = blinkstarts(pairs(:,1));
    blinkends = blinkends(pairs(:,2));
end
for k = 1:length(blinkstarts)
    eventgood( rawevents > blinkstarts(k) - NBlinkTO ...
        & rawevents < blinkends(k) + NBlinkTO ) = false;
end
for k = 1:(size(lfp_RecSegments,1) - 1)
    eventgood( rawevents > lfp_RecSegments(k,2) - NRecTO ...
        & rawevents < lfp_RecSegments(k+1,1) + NRecTO ) = false;
end
eventgood( rawevents < lfp_RecSegments(1,1) + NRecTO ) = false;
eventgood( rawevents > lfp_RecSegments(end,2) - NRecTO ) = false;
if ~any(eventgood)
    error('lfp_markGoodEvents4Dummies:oops', ...
        'oops');
end
