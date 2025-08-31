function ts = lfp_parseEyeTraces2(eyeX, eyeY, velo, threshes)
%LFP_PARSEEYETRACES for use with lfp_createEvents, marks eye movement
% related events.  
%ts = lfp_parseEyeTraces(eyeX, eyeY, velo, threshes)
% eyeX, eyeY, velo: filenums of X, Y, and velocity traces.
% threshes: a structure containing the thresholds used.  Thresholds
%   generally come in pairs, the idea being that the region between the
%   lower threshold and the upper threshold is considered ambiguous and
%   probably deserves special attention.  Note that <sacc> and <fix>
%   effectively constitute such a pair, but have kept their old names and
%   meanings from the previous version.
% The fields in the structure are:
%   blinkv1 - lower threshold for identifying blinks in velocity trace
%   blinkv2 - higher threshold for identifying blinks in velocity trace
%   blinkx - threshold for identifying blinks in X trace
%   blinky - threshold for identifying blinks in Y trace
%   blinkdur - the maximum duration of a blink
%   sacc - positive-going threshold for identifying saccades in velo
%   fix - negative-going threshold for identifying fixations in velo
%   blinkTO - time out around each blink during which lfp_EyeTabulation
%       discards events as artifactual
%   recTO - time out around each recording gap during which to discard
%       other events as artifactual
% A 'blink' is now defined as the following sequence of signals:
%   1. a period of blinkTO sec whose median X value is below blinkx
%   2. a "grace period" of blinkTO sec (there are no conditions on this
%   period)
%   3. a v peak that exceeds blinkv2 
%   4. a series of X samples (a) whose median exceeds blinkx and (b) whose
%   duration does not exceed blinkdur
%   5. another v peak that exceeds blinkv2
%   6. another "grace period" of blinkTO sec
%   7. another period of blinkTO sec whose median X value is below blinkx
% The algorithm for finding blinks is to find the first v peak (condition 3
% above), then test each of the other conditions in the following order:
% v peaks that fail condition 1 are considered "Lone End Blinks".
% If cond. 4b fails, it is a "Non Blink".  
% If cond. 4a fails, it is a "Non Blink".  
% If cond. 7 fails, it is a "Non Blink".
% In order to improve the effectiveness of majority vote logic for blink
% identification, I (DG) am adding (6-Apr-2008) a set of Y position
% criteria that echo the X position criteria, using thresh.blinky instead
% of threshes.bliinkx.
% A separate list of v peaks that exceed blinkv1 but that do not exceed
% blinkv2 is returned as "Ambiguous V Peak".
% A separate list of blinkx or blinky threshold crossings that are not
% represented in the above is returned as "Low V blinkx crossing".
% Returns a cell vector of timestamp column vectors, with one column for
% each of: 
%   Start Blink
%   End Blink
%   Start Saccade
%   End Saccade
%   Start Fixation
%   End Fixation
%   Lone End Blink
%   Non Blink
%   Ambiguous V Peak (sacc/blink)
%   Ambiguous V Peak (microsaccade/saccade)
%   Low V blinkx crossing
%
% Note that this function uses lfp_findPosThresh3, lfp_findNegThresh3
% instead of lfp_findPosThresh, lfp_findNegThresh (which are used by the
% original lfp_parseEyeTraces).  Also, it calls the script
% lfp_markGoodEvents2, which is therefore effectively part of the inline
% code.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

SaccThresh = threshes.sacc;
FixThresh = threshes.fix;
BlinkTimeout = threshes.blinkTO;
RecTimeout = threshes.recTO;
if ~ismember('blinkv', fieldnames(threshes))
    threshes.blinkv = 0;
end

% Derived constants:
NBlinkTimeout = round(BlinkTimeout / lfp_SamplePeriod);
NBlinkDur = round(threshes.blinkdur / lfp_SamplePeriod);
NRecTimeout = round(RecTimeout / lfp_SamplePeriod);

% Find ambiguous microsaccade/saccades:
ambigmicrosaccs = lfp_findPeaks(velo, FixThresh, SaccThresh);

% Find blinks.  
vmaxidx = lfp_findPeaks(velo, threshes.blinkv1);
blinkidx = vmaxidx(lfp_Samples{velo}(vmaxidx) > threshes.blinkv2);
ambigsaccblinks = vmaxidx(lfp_Samples{velo}(vmaxidx) <= threshes.blinkv2);
% At this point, blinkidx and ambigsaccblinks are indices into
% lfp_Samples{velo} of v peaks in the blink range and the ambiguous range
% respectively.
loneEndBlinks = [];
nonBlinks = [];
blinkstarts = [];
blinkends = [];
blinkidx2 = 1;  % index into blinkidx ("2nd order" blink index);
while blinkidx2 < length(blinkidx)
    bidx = blinkidx(blinkidx2);
    % condition 1:
    if median(lfp_Samples{eyeX}(bidx-2*NBlinkTimeout ...
            : bidx-NBlinkTimeout )) >= threshes.blinkx
        loneEndBlinks(end+1,1) = bidx;
        blinkidx2 = blinkidx2 + 1;
        continue
    end
    %condition 4b:
    if (blinkidx(blinkidx2+1) > bidx + NBlinkDur)
        nonBlinks(end+1,1) = bidx;
        blinkidx2 = blinkidx2 + 1;
        continue
    end
    % condition 4a:
    if median(lfp_Samples{eyeX}(bidx+1 : blinkidx(blinkidx2+1))) ...
            < threshes.blinkx
        nonBlinks(end+1,1) = bidx;
        blinkidx2 = blinkidx2 + 1;
        continue
    end
    % condition 7:
    if median(lfp_Samples{eyeX}( blinkidx(blinkidx2+1) + NBlinkTimeout ...
            : blinkidx(blinkidx2+1) + 2*NBlinkTimeout )) >= threshes.blinkx
        nonBlinks(end+1,1) = bidx;
        blinkidx2 = blinkidx2 + 1;
        continue
    end
    blinkstarts(end+1,1) = bidx;
    blinkends(end+1,1) = blinkidx(blinkidx2+1);
    blinkidx2 = blinkidx2 + 2;
end
if blinkidx2 == length(blinkidx)
    nonBlinks(end+1,1) = bidx;
end
lowvblinks = dg_tolSetdiff([ 
    lfp_findPosThresh3(eyeX, threshes.blinkx, 'sample')'
    lfp_findNegThresh3(eyeX, threshes.blinkx, 'sample')'
    ], ...
    [blinkstarts; blinkends; ambigsaccblinks; loneEndBlinks; nonBlinks], ...
    NBlinkTimeout );

failedBlinks = ...
    sort([lowvblinks; ambigsaccblinks; loneEndBlinks; nonBlinks]);

% Find saccade starts:
rawevents = reshape(lfp_findPosThresh3(velo, SaccThresh, 'sample'), [], 1);
lfp_markGoodEvents2;
saccstarts = rawevents(eventgood);

% Find saccade ends:
rawevents = reshape(lfp_findNegThresh3(velo, SaccThresh, 'sample'), [], 1);
lfp_markGoodEvents2;
saccends = rawevents(eventgood);

% Find fixation starts:
rawevents = reshape(lfp_findNegThresh3(velo, FixThresh, 'sample'), [], 1);
lfp_markGoodEvents2;
fixstarts = rawevents(eventgood);

% Find fixation ends:
rawevents = reshape(lfp_findPosThresh3(velo, FixThresh, 'sample'), [], 1);
lfp_markGoodEvents2;
fixends = rawevents(eventgood);

% Convert from sample indices to timestamps and construct result:
ts = {lfp_index2time(blinkstarts) lfp_index2time(blinkends) ...
        lfp_index2time(saccstarts) lfp_index2time(saccends) ...
        lfp_index2time(fixstarts) lfp_index2time(fixends) ...
        lfp_index2time(loneEndBlinks) lfp_index2time(nonBlinks) ...
        lfp_index2time(ambigsaccblinks) lfp_index2time(ambigmicrosaccs) ...
        lfp_index2time(lowvblinks) };
