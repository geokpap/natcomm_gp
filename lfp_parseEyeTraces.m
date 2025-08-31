function ts = lfp_parseEyeTraces(eyeX, eyeY, velo, threshes)
%LFP_PARSEEYETRACES for use with lfp_createEvents, marks eye movement
% related events.  
%ts = lfp_parseEyeTraces(eyeX, eyeY, velo, threshes)
% eyeX, eyeY, velo: filenums of X, Y, and velocity traces.
% threshes: a structure containing the thresholds used, whose elements are:
%   blinkv - optional threshold for identifying blinks in velocity trace
%   blinkx - threshold for identifying blinks in X trace  NOTE:  both
%       blinkx and blinkv must be satisfied in order to mark a blink.
%   sacc - positive-going threshold for identifying saccades in velo
%   fix - negative-going threshold for identifying fixations in velo
%   blinksrc - time window around the blink threshold crossings to search
%       for maximum "velocity" (used to pinpoint blink time)
%   blinkTO - time out around each blink during which to discard other
%       events as artifactual
%   recTO - time out around each recording gap during which to discard
%       other events as artifactual
% Returns a 2D array of timestamps with one column for each of:
%   Start Blink
%   End Blink
%   Start Saccade
%   End Saccade
%   Start Fixation
%   End Fixation

% Only the X trace is used to locate blinks.
% This function calls the script lfp_markGoodEvents, which is therefore
% effectively part of the inline code.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

BlinkThreshX = threshes.blinkx;
SaccThresh = threshes.sacc;
FixThresh = threshes.fix;
BlinkSearchWindow = threshes.blinksrc;
BlinkTimeout = threshes.blinkTO;
RecTimeout = threshes.recTO;
if ~ismember('blinkv', fieldnames(threshes))
    threshes.blinkv = 0;
end

% Derived constants:
NBlinkSearchWindow = floor(BlinkSearchWindow / lfp_SamplePeriod);
NBlinkTimeout = round(BlinkTimeout / lfp_SamplePeriod);
NRecTimeout = round(RecTimeout / lfp_SamplePeriod);

% Find blinks:
if BlinkThreshX < 0
    startxings = lfp_findNegThresh(eyeX, BlinkThreshX, 'sample');
    endxings = lfp_findPosThresh(eyeX, BlinkThreshX, 'sample');
else
    startxings = lfp_findPosThresh(eyeX, BlinkThreshX, 'sample');
    endxings = lfp_findNegThresh(eyeX, BlinkThreshX, 'sample');
end
if isempty(startxings)
    error('lfp_parseEyeTraces:noblinkstarts', ...
        'There were no blink starts found' );
    blinkstarts = [];
else
    blinkstarts = reshape( lfp_findlocalmax(velo, ...
        startxings - NBlinkSearchWindow, ...
        startxings + NBlinkSearchWindow), [], 1);
    if threshes.blinkv > 0
        blinkstarts = blinkstarts( ...
            lfp_Samples{velo}(blinkstarts) > threshes.blinkv );
    end
end
if isempty(endxings)
    error('lfp_parseEyeTraces:noblinkstarts', ...
        'There were no blink ends found' );
    blinkends = [];
else
    blinkends = reshape(lfp_findlocalmax(velo, ...
        endxings - NBlinkSearchWindow, ...
        endxings + NBlinkSearchWindow), [], 1);
    if threshes.blinkv > 0
        blinkends = blinkends( ...
            lfp_Samples{velo}(blinkends) > threshes.blinkv );
    end
end

% Find saccade starts:
rawevents = reshape(lfp_findPosThresh(velo, SaccThresh, 'sample'), [], 1);
lfp_markGoodEvents;
saccstarts = rawevents(eventgood);

% Find saccade ends:
rawevents = reshape(lfp_findNegThresh(velo, SaccThresh, 'sample'), [], 1);
lfp_markGoodEvents;
saccends = rawevents(eventgood);

% Find fixation starts:
rawevents = reshape(lfp_findNegThresh(velo, FixThresh, 'sample'), [], 1);
lfp_markGoodEvents;
fixstarts = rawevents(eventgood);

% Find fixation ends:
rawevents = reshape(lfp_findPosThresh(velo, FixThresh, 'sample'), [], 1);
lfp_markGoodEvents;
fixends = rawevents(eventgood);

% Convert from sample indices to timestamps and construct result:
ts = {lfp_index2time(blinkstarts) lfp_index2time(blinkends) ...
        lfp_index2time(saccstarts) lfp_index2time(saccends) ...
        lfp_index2time(fixstarts) lfp_index2time(fixends)};
