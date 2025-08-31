function [ts, saccwave] = lfp_parseEyeTraces4(eyeX, eyeY, velo, threshes)
%LFP_PARSEEYETRACES for use with lfp_createEvents and lfp_createWave via 
% adapters; marks eye movement related events.  
%ts = lfp_parseEyeTraces(eyeX, eyeY, velo, threshes)
% eyeX, eyeY, velo: filenums of X, Y, and velocity traces.
% threshes: a structure containing the thresholds used.  Thresholds
%   generally come in pairs, the idea being that the region between the
%   lower threshold and the upper threshold is considered ambiguous and
%   probably deserves special attention.  Note that <sacc> and <fix>
%   effectively constitute such a pair, but have kept their old names and
%   meanings from the previous version (they bracket the peak values listed
%   in <ambigmicrosaccs>).
% The fields in the structure are:
%   blinkTO - time out around each blink during which lfp_EyeTabulation
%       discards events as artifactual
%   blinkdur - the maximum duration of a blink
%   blinksrc - time interval within which successive candidate blink events
%       are considered to be all one event
%   blinkv1 - lower threshold for identifying blinks in velocity trace
%   blinkv2 - higher threshold for identifying blinks in velocity trace
%   blinkx - threshold for identifying blinks in X trace
%   blinky - threshold for identifying blinks in Y trace
%   recTO - time out around each recording gap during which to discard
%       other events as artifactual
%   sacc - the threshold below which the v trace is ignored for fitting
%       saccades
%   fix - the threshold above which we call it an <ambigmicrosaccs> and not
%       part of a fixation.
% A prototypical 'blink' is typified by the following sequence of
% conditions:
%   1. a period of blinkTO sec whose median X value is below blinkx and
%   median Y value is below blinky
%   2. a "grace period" of blinkTO sec (there are no conditions on this
%   period)
%   3. a v peak that exceeds blinkv2, or a threshold crossing of blinkx or
%   blinky
%   4. a series of samples whose duration does not exceed blinkdur
%   5. median X of that series of samples exceeds blinkx or median Y
%   exceeds blinky 
%   6. another v peak that exceeds blinkv2
%   7. another "grace period" of blinkTO sec
%   8. another period of blinkTO sec whose median X value is below
%   blinkx and median Y is below blinky
% The algorithm for finding blinks is to find candidate blink events
% consisting of a v peak that exceeds blinkv2, a blinkx crossing, or a
% blink y crossing (condition 3 above). Taking the first candidate blink
% event as a putative blink start, we test the conditions as follows:
%   a. Each v peak > blinkv1 within +/- 2*blinksrc of the event casts one
%   <vote_a>.  Note that extreme glitchiness may produce more than one
%   peak, making it possible to overcome vote deficiencies when examing the
%   putative blink end (see below).
%   b. Each v peak > blinkv2 within +/- 2*blinksrc of the event casts one
%   <vote_b>.  Together with a., EACH of the peaks that meet this condition
%   ends up cumulatively casting two votes.
%   u. Condition 1 is strictly required
%   v. The next candidate event satisfying conditions 3 and 4 is now
%   considered as a putative blink end event, and condition 5 must be met
%   for at least one of X and/or Y.  If the required conditions are not met
%   at this point, the tally of conditions met and failed up to this point
%   is logged, and the process is repeated with the next event satisfying
%   conditions 3 and 4.  If there are no such events, then the putative
%   Start Blink event is marked as a Big Saccade, and the the timestamp of
%   the Big Saccade is logged.
%   c. If the putative blink start has v > blinkv1, <vote_c> is 1.  If the
%   putative blink start has v > blinkv2, <vote_c> is 2. 
%   w. Condition 8 is also strictly required.  The sum of votes is then
%   tallied for a strict majority.  If a good blink end has not yet been
%   found, then a putative blink end that causes majority vote to fail at
%   this point is simply ignored, and the search for a blink end event
%   continues.  If a good blink end has previously been found, then the
%   blink is extended to each candidate blink event that still wins the
%   vote. When the majority vote logic finally fails, the last winner is
%   marked as End Blink, and its timestamp is logged together with the
%   tally of conditions met and failed.
% If no successful End Blink is found, then the putative Start Blink is
% marked as a Big Saccade, and its timestamp
% is logged.  In either case, the next candidate blink event that was
% not part of the blink or Big Saccade is used as the next putative
% Start Blink.
%
% Having found the blinks, the Failed Blinks are assumed to be "big" (i.e.
% not corrective or micro-) saccades.  Each is fitted to an analytic curve
% from which the saccade start and end times are determined.  Only the
% portion of each v peak that is above threshes.sacc is used for fitting
% (i.e. from vmax down to vmax/2).   
% <ts> is a cell vector of timestamp column vectors, with one column for
% each of: 
%   Start Blink
%   End Blink
%   Start Saccade
%   End Saccade
% <saccwave> is an array of sample data the same size as
% lfp_Samples{velo}.  To save execution time, it is produced simply by
% setting the points in the range 3*(-c(3):c(3)) around the v peak equal to
% the fitted curve, where c(3) is the time constant of the fit (i.e. the
% point at which the curve falls to 1/e of its peak value).  Therefore, if
% saccades are spaced more closely than three times the sum of their time
% constants, the later saccade will overwrite the earlier one.  FWIW, the
% savings in execution time is most likely negligible compared to the
% execution time of actually doing the curve fitting, but why change it
% now.
% The saccade start and end are computed as 2*[-c(3) c(3)] around the v
% peak.
%
% This function calls the script lfp_markGoodEvents, which is therefore
% effectively part of the inline code.
%
% 9-Apr-2008 DG moved all "...end" events one sample earlier so that events
% lists don't become confusing with multiple events in random order at a
% single time point.
%
% 15-Apr-2008 DG radically overhauled blink identification logic but didn't
% bother to change the name of the function, well, uh, 'cause the overhaul
% shoulda been there to start with ya know...
%
% 16-Apr-2008 TMD added a waitbar and changed the width at which saccades
% get marked.

%$Rev: 422 $
%$Date: 2023-09-08 14:25:38 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

lfp_declareGlobals;

if ~ismember('blinkv', fieldnames(threshes))
    threshes.blinkv = 0;
end

saccwave = zeros(size(lfp_Samples{velo}));

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
    unique( [
    vmaxidx(lfp_Samples{velo}(vmaxidx) > threshes.blinkv2)
    lfp_findPosThresh3(eyeX, threshes.blinkx, 'sample')'
    lfp_findNegThresh3(eyeX, threshes.blinkx, 'sample')'
    lfp_findPosThresh3(eyeY, threshes.blinky, 'sample')'
    lfp_findNegThresh3(eyeY, threshes.blinky, 'sample')'
    ] );

% blinkidx and vmaxidx are indices into lfp_Samples.
% blinkidx2, etc., are indices into blinkidx and vmaxidx
% ("2nd order" indices); they are used to obviate the need to
% search all of blinkidx and vmaxidx on each iteration.
blinkstarts = [];
blinkends = [];
failedblinks = [];
bigsaccs = [];
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
    vote_c = NaN;
    condition_1 = (median( ...
        lfp_Samples{eyeX}(max(1, bidx-2*NBlinkTO) : max(1, bidx-NBlinkTO) )) ...
        < threshes.blinkx) ...
        && (median( ...
        lfp_Samples{eyeY}(max(1, bidx-2*NBlinkTO) : max(1, bidx-NBlinkTO) )) ...
        < threshes.blinky) ;
    condition_8 = NaN;
    bendidx2 = blinkidx2 + 1;
    goodbendidx = 0;
    logmsg = '';
    while bendidx2 <= length(blinkidx) && ...
            blinkidx(bendidx2) <= bidx + NBlinkDur
        condition_5 = (median( ...
            lfp_Samples{eyeX}(bidx+1 : blinkidx(bendidx2)-1)) ...
            >= threshes.blinkx) ...
            || (median( ...
            lfp_Samples{eyeY}(bidx+1 : blinkidx(bendidx2)-1)) ...
            >= threshes.blinky) ;
        if condition_1 && condition_5
            % conditions met so far; keep evaluating:
            % This case always either increments bendidx2 or breaks the
                % 'while bendidx2' loop:
            vote_c = (lfp_Samples{velo}(blinkidx(bendidx2)) ...
                > threshes.blinkv1 ) ...
                + (lfp_Samples{velo}(blinkidx(bendidx2)) ...
                > threshes.blinkv2 );
            condition_8 = (median(lfp_Samples{eyeX}( ...
                blinkidx(bendidx2) + NBlinkTO ...
                : blinkidx(bendidx2) + 2*NBlinkTO )) ...
                < threshes.blinkx ) ...
                && (median(lfp_Samples{eyeY}( ...
                blinkidx(bendidx2) + NBlinkTO ...
                : blinkidx(bendidx2) + 2*NBlinkTO )) ...
                < threshes.blinky) ;
            if goodbendidx
                % This case will either increment bendidx2 or break the
                % 'while bendidx2' loop:
                if condition_8 && ...
                        (vote_a + vote_b + vote_c > 2)
                    % another good blink end, keep trying to extend
                    goodbendidx = blinkidx(bendidx2);
                    logmsg = [ logmsg sprintf(...
                        'Extended blink to %.6f; a:%d b:%d c:%d\n', ...
                        lfp_index2time(goodbendidx), ...
                        vote_a, vote_b, vote_c )];
                    bendidx2 = bendidx2 + 1;
                else
                    % existing goodbendidx is the final Blink End
                    break
                end
            else
                % This always increments bendidx2:
                if condition_8 && ...
                        (vote_a + vote_b + vote_c > 2)
                    % First good blink end for this blink start - note that
                    % the "probably a bug" has been fixed here.
                    goodbendidx = blinkidx(bendidx2);
                    logmsg = [ logmsg sprintf(...
                        'Found blink end at %.6f; a:%d b:%d c:%d\n', ...
                        lfp_index2time(goodbendidx), ...
                        vote_a, vote_b, vote_c )];
                end
                bendidx2 = bendidx2 + 1;
            end
        else
            % Initial conditions not met; log as Big Saccade.
            % This case always either increments bendidx2 or breaks:
            logmsg = [ logmsg sprintf(...
                'Initial conditions not met;\nputative blink end is Big Saccade at %.6f\n', ...
                lfp_index2time(blinkidx(bendidx2)) ) ];
            logmsg = [ logmsg sprintf(...
                'Candidate blink timestamps: %s\n', ...
                mat2str(lfp_index2time(blinkidx(blinkidx2:bendidx2))) )];
            logmsg = [ logmsg sprintf( ...
                'cond1:%d cond5:%d cond8:%d a:%d b:%d c:%d\n', ...
                condition_1, condition_5, condition_8, ...
                vote_a, vote_b, vote_c )];
            if goodbendidx
                % found bad blink end after good blink end, quit looking:
                break;
            else
                % found nothing but bad blink ends so far, keep looking:
                bendidx2 = bendidx2 + 1;
            end
        end
    end
    % The Search For A Blink End is now over, for better or worse.
    % blinkidx2 gets incremented by at least 1, one way or another:
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
        logmsg = [ logmsg sprintf(...
            'No good blink ends;\nputative blink start is Big Saccade at %.6f\n', ...
            lfp_index2time(bidx) )];
        logmsg = [ logmsg sprintf(...
            'Candidate blink timestamps: %s\n', ...
            mat2str(lfp_index2time(blinkidx(blinkidx2:bendidx2))) )];
        logmsg = [ logmsg sprintf( ...
            'cond8:%d a:%d b:%d c:%d\n', condition_8, ...
            vote_a, vote_b, vote_c )];
        % lfp_log opens and closes the log file, so it should not be called
        % often:
        lfp_log(logmsg);
        % Try again with the next candidate blink event:
        blinkidx2 = blinkidx2 + 1;
    end
end

% Find big saccades and remove those that are too close to blinks:
bigsaccs = dg_tolSetdiff(vmaxidx, [blinkstarts; blinkends], NBlinkTO);
% NOTE: this function no longer calls the script
% lfp_markGoodEvents3 at this point.
% bigsaccs = rawevents(eventgood);
saccstarts = zeros(size(bigsaccs));
saccends = zeros(size(bigsaccs));

% Find saccade starts and ends:
h = waitbar(0,'Finding saccade starts and ends...');
for saccidx = 1:length(bigsaccs)
    vmax = lfp_Samples{velo}(bigsaccs(saccidx));
    % Find the half-velocity points:
    startpt = bigsaccs(saccidx);
    while startpt > 0 && lfp_Samples{velo}(startpt) > threshes.sacc
        startpt = startpt - 1;
    end
    endpt = bigsaccs(saccidx);
    while endpt <= numel(lfp_Samples{velo}) && ...
            lfp_Samples{velo}(endpt) > threshes.sacc
        endpt = endpt + 1;
    end
    % skip it if this "sacc" is not actually the max in the interval
    if any(lfp_Samples{velo}(startpt:endpt) > vmax)
        continue
    end
    % must normalize because 'fit' can't cope with radically different scales
    % along different dimensions of the param search space:
    x = (startpt:endpt) - bigsaccs(saccidx);
    y = lfp_Samples{velo}(startpt:endpt)/vmax;
    cfun = fit(x', y', 'gauss1', ...
        fitoptions( 'Method','NonlinearLeastSquares', ...
        'Upper', [Inf Inf 0.02/lfp_SamplePeriod] ));
    c = coeffvalues(cfun);
    % De-normalize the vertical scale param:
    c(1) = c(1) * vmax;
    smarkingwidth = 1.4; % used to = 2, then 1.8 changed by TMD 4/17/08
    emarkingwidth = 1.7;
    saccstarts(saccidx) = round(c(2) - smarkingwidth*c(3)) + bigsaccs(saccidx) - 1;
    saccends(saccidx) = round(c(2) + emarkingwidth*c(3)) + bigsaccs(saccidx) - 1;
    fitrange = round(-3*c(3)):round(3*c(3));
    saccwaverange = fitrange + bigsaccs(saccidx) - 1;
    if saccwaverange(1) < 1
        fitrange = fitrange(saccwaverange>0);
        saccwaverange = saccwaverange(saccwaverange>0);
    end
    saccwave(saccwaverange) = c(1)*exp(-((fitrange' - c(2)) / c(3)).^2);
    waitbar(saccidx/length(bigsaccs));
end
close(h)

% Filter saccade starts:
rawevents = saccstarts;
lfp_markGoodEvents4;
saccstarts = rawevents(eventgood);

% Filter saccade ends:
rawevents = saccends;
lfp_markGoodEvents4;
saccends = rawevents(eventgood);

% Convert from sample indices to timestamps and construct result:
ts = {lfp_index2time(blinkstarts) lfp_index2time(blinkends-1) ...
    lfp_index2time(saccstarts) lfp_index2time(saccends-1) };


