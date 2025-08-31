function ts = lfp_parseEyeTraces3(eyeX, eyeY, velo, threshes)
% NOT TESTED/DEBUGGED!
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
%   blinksrc - time interval within which successive candidate blink events
%   are considered to be all one event
% A prototypical 'blink' is typified by the following sequence of
% conditions:
%   1. a period of blinkTO sec whose median X value is below blinkx and
%   median Y value is below blinky
%   2. a "grace period" of blinkTO sec (there are no conditions on this
%   period)
%   3. a v peak that exceeds blinkv2, or a threshold crossing of blinkx or
%   blinky
%   4. a series of samples whose duration does not exceed blinkdur
%   5. median X of that series of samples exceeds blinkx and median Y
%   exceeds blinky 
%   6. another v peak that exceeds blinkv2
%   7. another "grace period" of blinkTO sec
%   8. another period of blinkTO sec whose median X value is below
%   blinkxand median Y is below blinky
% The algorithm for finding blinks is to find candidate blink events
% consisting of a v peak that exceeds blinkv2, a blinkx crossing, or a
% blink y crossing (condition 3 above); any series of such events occurring
% within blinksrc of each other is considered to be a single event.
% Taking the first candidate blink event as a putative blink start, we
% test each of the following conditions and come to a conclusion by
% majority vote:
%   a. Each v peak > blinkv1 within +/- 2*blinksrc of the event casts one
%   vote
%   b. Each v peak > blinkv2 within +/- 2*blinksrc of the event casts one
%   vote (together with a., these peaks end up cumulatively casting two
%   votes each)
%   c. Condition 1 casts one vote each for X and Y
%   d. The next candidate event satisfying conditions 3 and 4 is now
%   considered as a putative blink end event, and condition 5 casts a vote
%   for each of X and Y.  There must be a majority vote accumulated up to
%   this point in order to continue considering the putative blink end
%   event; if there is not, the tally of conditions met and failed up to
%   this point is logged, and the process is repeated with the next event
%   satisfying conditions 3 and 4.  If there are no such events, then the
%   putative Start Blink event is marked as a Failed Blink, and the
%   the timestamp of the Failed Blink is logged.
%   e. Condition 8 also casts one vote for each of X and Y.  If a good
%   blink end has not yet been found, then a putative blink end that causes
%   majority vote to fail at this point is simply ignored, and the search
%   for a blink end event continues.  If a good blink end has been found,
%   then an attempt is made to extend the blink to the next candidate blink
%   event, and majority vote logic is again applied to the extended blink.
%   When the majority vote logic finally fails, the last winner is marked
%   as End Blink, and its timestamp is logged together with the tally of
%   conditions met and failed.
% If no successful End Blink is found, then the putative Start Blink is
% marked as a Failed Blink, and the timestamp of the Failed Blink
% Start is logged.  In either case, the next candidate blink event that was
% not part of the blink or Failed Blink is used as the next putative
% Start Blink.
% Returns a cell vector of timestamp column vectors, with one column for
% each of: 
%   Start Blink
%   End Blink
%   Start Saccade
%   End Saccade
%   Start Fixation
%   End Fixation
%   Failed Blink
%   Ambiguous V Peak (microsaccade/saccade)
%
% Note that this function uses lfp_findPosThresh3, lfp_findNegThresh3
% instead of lfp_findPosThresh, lfp_findNegThresh (which are used by the
% original lfp_parseEyeTraces).  Also, it calls the script
% lfp_markGoodEvents3, which is therefore effectively part of the inline
% code.
%
% 9-Apr-2008 DG moved all "...end" events one sample earlier so that events
% lists don't become confusing with multiple events in random order at a
% single time point.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~ismember('blinkv', fieldnames(threshes))
    threshes.blinkv = 0;
end

% Derived constants:
NBlinkTO = round(threshes.blinkTO / lfp_SamplePeriod);
NBlinkDur = round(threshes.blinkdur / lfp_SamplePeriod);
NRecTO = round(threshes.recTO / lfp_SamplePeriod);
NBlinkSrc = round(threshes.blinksrc / lfp_SamplePeriod);

% Find ambiguous microsaccade/saccades:
ambigmicrosaccs = lfp_findPeaks(velo, threshes.fix, threshes.sacc);

%
% Find blinks.  
%
NBlinkTimeout2 = round(NBlinkTO/2);

% find sample indices of candidate blink events <blinkidx>:
vmaxidx = lfp_findPeaks(velo, threshes.blinkv1);
blinkidx = ...
    dg_tolUnion(vmaxidx(lfp_Samples{velo}(vmaxidx) > threshes.blinkv2), ...
    [
    lfp_findPosThresh3(eyeX, threshes.blinkx, 'sample')'
    lfp_findNegThresh3(eyeX, threshes.blinkx, 'sample')'
    lfp_findPosThresh3(eyeY, threshes.blinky, 'sample')'
    lfp_findNegThresh3(eyeY, threshes.blinky, 'sample')'
    ], ...
    NBlinkSrc );

% blinkidx and vmaxidx are indices into lfp_Samples.
% blinkidx2, etc., are indices into blinkidx and vmaxidx
% ("2nd order" indices); they are used to obviate the need to
% search all of blinkidx and vmaxidx on each iteration.
blinkstarts = [];
blinkends = [];
failedblinks = [];
blinkidx2 = 1;  % index into blinkidx ("2nd order" blink index);
vmaxidx2 = 1;  % starting index into vmaxidx for votes a & b search
while blinkidx2 < length(blinkidx)
    bidx = blinkidx(blinkidx2); % sample # of current putative blink start
    % find range of vmaxidx to search for votes a and b (aka
    % vmaxidx2a:vmaxidx2b).  We need to start at vmaxidx(vmaxidx2a) =
    % blinkidx(blinkidx2), 
    vmaxidx2a = vmaxidx2;
    while vmaxidx(vmaxidx2a) <= bidx - 2*NBlinkSrc
        vmaxidx2a = vmaxidx2a + 1;
    end
    vmaxidx2b = vmaxidx2a;
    while vmaxidx2b <= length(vmaxidx) && ...
            vmaxidx(vmaxidx2b) < bidx + 2*NBlinkSrc
        vmaxidx2b = vmaxidx2b + 1;
    end
    % so now we've gone one too far to the right in finding vmaxidx2b:
    vmaxidx2b = vmaxidx2b - 1;
    vmaxidx2 = vmaxidx2a;   % update starting point for next time
    vote_a = sum(lfp_Samples{velo}(vmaxidx(vmaxidx2a:vmaxidx2b)) ...
        > threshes.blinkv1 );
    vote_b = sum(lfp_Samples{velo}(vmaxidx(vmaxidx2a:vmaxidx2b)) ...
        > threshes.blinkv2 );
    vote_c = (median( ...
        lfp_Samples{eyeX}(bidx-2*NBlinkTO : bidx-NBlinkTO )) ...
        < threshes.blinkx) + (median( ...
        lfp_Samples{eyeY}(bidx-2*NBlinkTO : bidx-NBlinkTO )) ...
        < threshes.blinky) ;
    bendidx2 = blinkidx2 + 1;
    goodbendidx = 0;
    logmsg = '';
    while bendidx2 <= length(blinkidx) && ...
            blinkidx(bendidx2) <= bidx + NBlinkDur
        vote_d = (median( ...
            lfp_Samples{eyeX}(bidx+1 : blinkidx(bendidx2)-1)) ...
            >= threshes.blinkx) + (median( ...
            lfp_Samples{eyeY}(bidx+1 : blinkidx(bendidx2)-1)) ...
            >= threshes.blinky) ;
        if vote_a + vote_b + vote_c + vote_d > ...
                ~vote_a + ~vote_b + ~vote_c + ~vote_d
            % subtotal through vote d has passed; keep evaluating:
            vote_e = (median(lfp_Samples{eyeX}( ...
                blinkidx(bendidx2) + NBlinkTO ...
                : blinkidx(bendidx2) + 2*NBlinkTO )) ...
                < threshes.blinkx ) + (median(lfp_Samples{eyeY}( ...
                blinkidx(bendidx2) + NBlinkTO ...
                : blinkidx(bendidx2) + 2*NBlinkTO )) ...
                < threshes.blinky) ;
            if goodbendidx
                if ~(vote_a + vote_b + vote_c + vote_d + vote_e > ...
                        ~vote_a + ~vote_b + ~vote_c + ~vote_d + ~vote_e )
                    % existing goodbendidx is the final Blink End
                    break
                else
                    % another good blink end, keep trying to extend
                    goodbendidx = blinkidx(bendidx2);
                    logmsg = [ logmsg sprintf(...
                        'Extended blink to %.6f; a:%d b:%d c:%d d:%d e:%d\n', ...
                        lfp_index2time(goodbendidx), ...
                        vote_a, vote_b, vote_c, vote_d, vote_e )];
                    bendidx2 = bendidx2 + 1;
                end
            else
                % First good blink end for this blink start - but note that
                % although vote_e has been calculated, it did not
                % contribute to this conclusion.  That's probably a bug.
                goodbendidx = blinkidx(bendidx2);
                logmsg = [ logmsg sprintf(...
                    'Found blink end at %.6f; a:%d b:%d c:%d d:%d e:%d\n', ...
                    lfp_index2time(goodbendidx), ...
                    vote_a, vote_b, vote_c, vote_d, vote_e )];
                bendidx2 = bendidx2 + 1;
            end
        else
            logmsg = [ logmsg sprintf(...
                'Skipping evt at %.6f; a:%d b:%d c:%d d:%d\n', ...
                lfp_index2time(blinkidx(bendidx2)), ...
                vote_a, vote_b, vote_c, vote_d ) ];
            bendidx2 = bendidx2 + 1;
        end
    end
    if goodbendidx
        blinkstarts(end+1,1) = bidx;
        blinkends(end+1,1) = goodbendidx;
        % One way or another, bendidx2 is now pointing at the first
        % candidate blink event that is not be part of the blink, so
        % this becomes our next starting point for blink detection:
        blinkidx2 = bendidx2;
        % skip forward past blinkTO:
        while blinkidx2 <= length(blinkidx) ...
                && blinkidx(blinkidx2) < (goodbendidx + NBlinkTO)
            blinkidx2 = blinkidx2 + 1;
        end
    else
        % There was no good blink end event at all.
        failedblinks(end+1, 1) = bidx;
        logmsg = [ logmsg sprintf(...
            'Failed Blink marked at %.6f\n', lfp_index2time(bidx) )];
        logmsg = [ logmsg sprintf(...
            'Candidate blink timestamps: %s', ...
            mat2str(lfp_index2time(blinkidx(blinkidx2:bendidx2))) )];
        % lfp_log opens and closes the log file, so it should not be called
        % often:
        lfp_log(logmsg);
        % Try again with the next candidate blink event:
        blinkidx2 = blinkidx2 + 1;
    end
end

% Find saccade starts:
rawevents = reshape(lfp_findPosThresh3(velo, threshes.sacc, 'sample'), [], 1);
lfp_markGoodEvents3;
saccstarts = rawevents(eventgood);

% Find saccade ends:
rawevents = reshape(lfp_findNegThresh3(velo, threshes.sacc, 'sample'), [], 1);
lfp_markGoodEvents3;
saccends = rawevents(eventgood);

% Find fixation starts:
rawevents = reshape(lfp_findNegThresh3(velo, threshes.fix, 'sample'), [], 1);
lfp_markGoodEvents3;
fixstarts = rawevents(eventgood);

% Find fixation ends:
rawevents = reshape(lfp_findPosThresh3(velo, threshes.fix, 'sample'), [], 1);
lfp_markGoodEvents3;
fixends = rawevents(eventgood);

% Convert from sample indices to timestamps and construct result:
ts = {lfp_index2time(blinkstarts) lfp_index2time(blinkends-1) ...
        lfp_index2time(saccstarts) lfp_index2time(saccends-1) ...
        lfp_index2time(fixstarts) lfp_index2time(fixends-1) ...
        lfp_index2time(failedblinks) lfp_index2time(ambigmicrosaccs) };
