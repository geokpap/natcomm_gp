function [ts, saccwave, mergethresh, revtrimthresh] = ...
    lfp_parseEyeTraces6(eyeX, eyeY, velo, accel, threshes, vpeakthresh)
%lfp_parseEyeTraces6 for use with lfp_createEvents and lfp_createWave via 
% adapters; marks eye movement related events.  
%[ts, saccwave, mergethresh, revtrimthresh] = ...
%   lfp_parseEyeTraces6(eyeX, eyeY, velo, accel, threshes)

% MAINTENANCE NOTES:  
% 1) The strategic outline of this program is documented in
%	"lfp_parseEyeTraces6.doc".  Any changes made to the code should entail
%	corresponding changes to  "lfp_parseEyeTraces6.doc".
% 2) This code is Theresa-specific because of hard-coded references to
%   events 122 and 123, but could easily be generalized
% 3) VARIABLE NAMES:  variables that contain sample numbers or have nothing
%   whatsoever to do with sample numbers have no special suffixes.
%   However, indices into lists of sample numbers will have the suffix idx;
%   indices into lists of indices into sample numbers will be idx2; etc.
%   Timestamps in seconds (rather than samples) will have the suffix TS.
%   Iterators that index something that has nothing to do with sample
%   numbers have the suffix num.

% eyeX, eyeY, velo, accel: filenums of X, Y, velocity and acceleration
%   traces ("velocity" here really means "speed", i.e. the magnitude of the
%   velocity vector)
% threshes: a structure containing the thresholds used. The fields in the
%   structure are:
%   blinkTO - time out around each blink during which lfp_parseEyeTraces6
%       discards velocity peaks from consideration as saccades; also used
%       in identifying blinks (see below)
%   blinkdur - the maximum duration of a blink
%   blinksrc - time interval within which successive candidate blink events
%       are considered to be all one event
%   blinkv1 - lower threshold to help in identifying "low velocity" blinks
%       in velocity trace and for identifying "big" (i.e. not corrective)
%       saccades
%   blinkv2 - higher threshold for identifying blinks in velocity trace
%   blinkx - threshold for identifying blinks in X trace
%   blinky - threshold for identifying blinks in Y trace
%   minfixdur - duration threshold below which "fixations" are not marked
%   recTO - time out around each recording gap during which to discard
%       other events as artifactual
%   sacc - the threshold below which the v trace is ignored for fitting
%       saccades when no local minima are found above that threshold
%   vpeakratiomergethresh - optional value; if supplied and not empty, then
%       it is used as the value of mergethresh instead of computing it from
%       the v peak ratios histogram.
%   vpeakratiorevtrimthresh - optional value; if supplied and not empty,
%       then it is used as the value of revtrimthresh instead of computing
%       it from the v peak ratios histogram.
% <vpeakthresh>: 'blinkv1' or 'maxfixv': There are two thresholds that can 
%   be used to construct the initial list of velocity peaks (vpeaks) which
%   are candidate saccades and blinks: blinkv1 or maxfixv 
%   (1.75 * threshes.sacc). For Theresa M1 requires maxfixv and M2 requires
%   blinkv1.
%
% A prototypical 'blink' is typified by the following sequence of
% conditions:
%   1. a period of 2*blinkTO sec whose median X value is below blinkx and
%   median Y value is below blinky
%   2. (There is no longer any condition 2.)
%   3. a v peak that exceeds blinkv2, or a threshold crossing of blinkx or
%   blinky
%   4. a series of samples whose duration does not exceed blinkdur
%   5. median X of that series of samples exceeds blinkx or median Y
%   exceeds blinky 
%   6. another v peak that exceeds blinkv2
%   7. (There is no longer any condition 7.)
%   8. another period of 2*blinkTO sec whose median X value is below
%   blinkx and median Y is below blinky
% The algorithm for finding blinks is to find candidate blink events,
% defined as either a v peak that exceeds blinkv2, a blinkx crossing, or a
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
% not corrective or micro-) saccades and are accordingly called <bigsaccs>
% (the elements of which are sample indices).
% Before further processing, the <bigsaccs> list is scanned for runs of v
% peaks that occur within 50 ms (hard-coded value of <minsaccint>) of each
% other, and the runs are handled as follows.  The ratio of each peak v in
% the run to the first peak v in the run is computed, and if it is above
% <mergethresh> (see below), then the peak is considered to be part of a
% "split" saccade together with the first peak, and they are fitted en
% masse; otherwise, the entire tail of the run from the first v peak that
% is below <mergethresh> onwards is ignored, effectively deleting those
% (micro-) saccades from the list of saccades to be fitted. If more than
% three peaks are to be fitted en masse, that by definition makes it a LNM
% (Loch Ness Monster), which means it has to be fit all the way to the tip
% of its boat-smashing, terrorizing tail.
% <mergethresh> is found by constructing a histogram of the ratios of peak
% v's described above for peaks during the targets-on period ({122 123}
% evtbounds), applying enough smoothing to produce no more than 3 minima
% between the global maximum bin and 1, and then finding the first
% (left-most) such minimum.  A "reverse trimming threshold" <revtrimthresh>
% is obtained from the same smoothed histogram by finding the larger of
% either the first bin whose value is below 0.1, or 3 times <mergethresh>.
% If values for <mergethresh> or <revtrimthresh> cannot be extracted from
% the histogram, they default to 0.6 and 2 respectively.
% Each individual or "split" saccade is then fitted by an analytic curve
% from which the saccade start and end times are determined.  Only the
% portion of each v peak that is above threshes.sacc AND is at or above the
% highest local minimum before the first or after the last peak in the v
% trace is used for fitting.  Saccades that are too close to blinks are
% eliminated at the v-peak-finding stage, so this means that sacc starts
% (but not sacc ends) are guaranteed to be more than blinkTO before a blink
% start, and sacc ends (but not starts) are guaranteed to be more than
% blinkTO after a blink end.  recTO works similarly for rec starts and
% ends.  It also means that there are guaranteed to be equal numbers of
% sacc starts and sacc ends.
%
% Finally, fixations are initially considered to be "all the rest", and
% then potentially shortened on both ends to meet position and velocity
% criteria. "All the rest" is actually not as simple as it sounds.
% Practically speaking, a fixation can start:
%   1 sample after sacc end
%   1 blinkTO after blink end
%   1 recTO after rec start
% and a fixation can end:
%   1 sample before sacc start
%   1 blinkTO before blink start
%   1 recTO before rec end
% The first pass of parsing inevitably results in some poorly parsed eye
% movements.  These are discriminated on the basis of the long
% dimension of the bounding rectangle that just fits the eye trace during a
% fixation, and the maximum eye speed during the fixation.
%
% <ts> is a cell vector of timestamp column vectors, with one column for
% each of: 
%   Start Blink
%   End Blink
%   Start Simple Saccade
%   End Simple Saccade
%   Start Complex Saccade
%   End Complex Saccade
%   Post-Blink Fixation Starts
%   Post-Blink Fixation Ends
%   Other Fixation Starts
%   Other Fixation Ends
%   Bad Fixation Starts (when adjustment fails)
%   Bad Fixation Ends (when adjustment fails)
%   Start Merged Saccade Series
%   End Merged Saccade Series
%   Fit Starts
%   Fit Ends
%   Fixations for Adjustment (marked at midpoint of unadjusted fixation)
%   Old Fixation Starts (before adjusting fixations)
%   Old Fixation Ends (before adjusting fixations)
%   Blink-Terminated Sacc Starts
%   Blink-Terminated Sacc Ends
%
% <saccwave> is the waveform of the compound Gaussian curves fitted to the
%   saccades.
% <mergethresh>, <revtrimthresh> - values actually used for parsing.

% NOTES
% Basically a fatter, more conciliatory rendition of lfp_parseEyeTraces5.
% Surreptitiously saves otherfixdiams and blinkfixdiams in eponymous .mat
% files in CWD.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

lfp_log('Starting lfp_parseEyeTraces6');

Nsamples = numel(lfp_Samples{velo});
saccwave = zeros(size(lfp_Samples{velo}));

% Derived constants:
NBlinkTO = round(threshes.blinkTO / lfp_SamplePeriod);
NBlinkDur = round(threshes.blinkdur / lfp_SamplePeriod);
NRecTO = round(threshes.recTO / lfp_SamplePeriod);
NBlinkSrc = round(threshes.blinksrc / lfp_SamplePeriod);
NMinFixDur = round(threshes.minfixdur / lfp_SamplePeriod);
maxfixv = 1.75 * threshes.sacc;

%
% Find blinks.  
%
NBlinkTimeout2 = round(NBlinkTO/2);

% Find sample indicies of candidate velocity peaks
% This is the master list for both blinks and saccades.
h = waitbar(0,'Finding candidate eye events...');
switch vpeakthresh
    case 'blinkv1'
        vpeaks = lfp_findPeaks(velo, threshes.blinkv1);
    case 'maxfixv'
        vpeaks = lfp_findPeaks(velo, maxfixv);
    otherwise
        error('lfp_parseEyeTraces6:badvpeakthresh', ...
            ['The option "' vpeakthresh '" for is not recognized.'] );
end
waitbar(1/5, h);
fPTx = lfp_findPosThresh3(eyeX, threshes.blinkx, 'sample')';
waitbar(2/5, h);
fNTx = lfp_findNegThresh3(eyeX, threshes.blinkx, 'sample')';
waitbar(3/5, h);
fPTy = lfp_findPosThresh3(eyeY, threshes.blinky, 'sample')';
waitbar(4/5, h);
fNTy = lfp_findNegThresh3(eyeY, threshes.blinky, 'sample')';
waitbar(5/5, h);
% find sample indices of candidate blink events <canblinks>:
canblinks = ...
    unique( [
    vpeaks(lfp_Samples{velo}(vpeaks) > threshes.blinkv2)
    fPTx
    fNTx
    fPTy
    fNTy
    ] );
close(h);

% canblinks and vpeaks are indices into lfp_Samples.
% canblinksidx, etc., are indices into canblinks and vpeaks; they are used
% to obviate the need to search all of canblinks and vpeaks on each
% iteration.
blinkstarts = [];
blinkends = [];
canblinksidx = 1;  % index into canblinks
vpeaksidx = 1;  % starting index into vpeaks for votes a & b search
h = waitbar(0,'Finding blinks...');
while canblinksidx < length(canblinks)
    bsamp = canblinks(canblinksidx); % sample # of current putative blink start
    % find range of vpeaks to search for votes a and b (aka
    % vpeaksidx_a:vpeaksidx_b).  We need to start at vpeaks(vpeaksidx_a) =
    % canblinks(canblinksidx), 
    vpeaksidx_a = vpeaksidx;
    while vpeaks(vpeaksidx_a) <= bsamp - 2*NBlinkSrc
        vpeaksidx_a = vpeaksidx_a + 1;
    end
    vpeaksidx_b = vpeaksidx_a;
    while vpeaksidx_b <= length(vpeaks) && ...
            vpeaks(vpeaksidx_b) < bsamp + 2*NBlinkSrc
        vpeaksidx_b = vpeaksidx_b + 1;
    end
    % so now we've gone one too far to the right in finding vpeaksidx_b:
    vpeaksidx_b = vpeaksidx_b - 1;
    vpeaksidx = vpeaksidx_a;   % update starting point for next time
    vote_a = sum(lfp_Samples{velo}(vpeaks(vpeaksidx_a:vpeaksidx_b)) ...
        > threshes.blinkv1 );
    vote_b = sum(lfp_Samples{velo}(vpeaks(vpeaksidx_a:vpeaksidx_b)) ...
        > threshes.blinkv2 );
    vote_c = NaN;
    condition_1 = (median( ...
        abs(lfp_Samples{eyeX}( max(1, bsamp-2*NBlinkTO) : bsamp ))) ...
        < threshes.blinkx ) ...
        && (median( ...
        abs(lfp_Samples{eyeY}( max(1, bsamp-2*NBlinkTO) : bsamp ))) ...
        < threshes.blinky ) ; % time before vpeak is < blinkx & blinky
    condition_8 = NaN;
    bendidx = canblinksidx + 1; %index into blinkindex + 1 to look for end
    goodbend = 0;
    logmsg = sprintf('Cumulative comment:\n');
    while bendidx <= length(canblinks) && ...
            canblinks(bendidx) <= bsamp + NBlinkDur
        condition_5 = (median( ...
            abs(lfp_Samples{eyeX}(bsamp+1 : canblinks(bendidx)-1))) ...
            >= threshes.blinkx) ...
            || (median( ...
            abs(lfp_Samples{eyeY}(bsamp+1 : canblinks(bendidx)-1))) ...
            >= threshes.blinky) ; % time after vpeak is > blinkx | blinky
        if condition_1 && condition_5
            % conditions met so far; keep evaluating:
            % This case always either increments bendidx or breaks the
                % 'while bendidx' loop:
            vote_c = (lfp_Samples{velo}(canblinks(bendidx)) ...
                > threshes.blinkv1 ) ...
                + (lfp_Samples{velo}(canblinks(bendidx)) ...
                > threshes.blinkv2 );
            condition_8 = (median(abs(lfp_Samples{eyeX}( ...
                canblinks(bendidx) ...
                : min( Nsamples, canblinks(bendidx) + 2*NBlinkTO ) )) ) ...
                < threshes.blinkx ) ...
                && (median(abs(lfp_Samples{eyeY}( ...
                canblinks(bendidx) ...
                : min( Nsamples, canblinks(bendidx) + 2*NBlinkTO ) ))) ...
                < threshes.blinky) ; % time after ending vpeak < blinkx & blinky
            if goodbend
                % This case will either increment bendidx or break the
                % 'while bendidx' loop:
                if condition_8 && ...
                        (vote_a + vote_b + vote_c > 2)
                    % another good blink end, keep trying to extend
                    goodbend = canblinks(bendidx);
                    logmsg = [ logmsg sprintf(...
                        'Extended blink to %.6f; a:%d b:%d c:%d\n', ...
                        lfp_index2time(goodbend), ...
                        vote_a, vote_b, vote_c )];
                    bendidx = bendidx + 1;
                else
                    % existing goodbend is the final Blink End
                    break
                end
            else
                % This always increments bendidx:
                if condition_8 && ...
                        (vote_a + vote_b + vote_c > 2)
                    % First good blink end for this blink start - note that
                    % the "probably a bug" has been fixed here.
                    goodbend = canblinks(bendidx);
                    logmsg = [ logmsg sprintf(...
                        'Found blink end at %.6f; a:%d b:%d c:%d\n', ...
                        lfp_index2time(goodbend), ...
                        vote_a, vote_b, vote_c )];
                end
                bendidx = bendidx + 1;
            end
        else
            % Initial conditions not met; log as Big Saccade.
            % This case always either increments bendidx or breaks:
            logmsg = [ logmsg sprintf(...
                'Initial conditions not met;\nputative blink end is Big Saccade at %.6f\n', ...
                lfp_index2time(canblinks(bendidx)) ) ];
            logmsg = [ logmsg sprintf(...
                'Candidate blink timestamps: %s\n', ...
                mat2str(lfp_index2time(canblinks(canblinksidx:bendidx))) )];
            logmsg = [ logmsg sprintf( ...
                'cond1:%d cond5:%d cond8:%d a:%d b:%d c:%d\n', ...
                condition_1, condition_5, condition_8, ...
                vote_a, vote_b, vote_c )];
            if goodbend
                % found bad blink end after good blink end, quit looking:
                break;
            else
                % found nothing but bad blink ends so far, keep looking:
                bendidx = bendidx + 1;
            end
        end
    end
    % The Search For A Blink End is now over, for better or worse.
    % canblinksidx gets incremented by at least 1, one way or another:
    if goodbend
        if length(logmsg) > 20
            lfp_log(logmsg);
        end
        blinkstarts(end+1,1) = bsamp;
        blinkends(end+1,1) = goodbend - 1;
        % One way or another, bendidx is now pointing at the first
        % candidate blink event that is not to be part of the blink, so
        % this becomes our next starting point for blink detection.  It
        % could also be pointing past the end of canblinks, however, in
        % which case the canblinksidx loop will exit.
        canblinksidx = bendidx;
        % skip forward past blinkTO:
        while canblinksidx <= length(canblinks) ...
                && canblinks(canblinksidx) < (goodbend + NBlinkTO)
            canblinksidx = canblinksidx + 1;
        end
    else
        % There was no good blink end event at all.
        logmsg = [ logmsg sprintf(...
            'No good blink ends;\nputative blink start is Big Saccade at %.6f\n', ...
            lfp_index2time(bsamp) )];
        logmsg = [ logmsg sprintf(...
            'Candidate blink timestamps: %s\n', ...
            mat2str(lfp_index2time(canblinks( ...
            canblinksidx:min(bendidx, length(canblinks)) ))) )];
        logmsg = [ logmsg sprintf( ...
            'cond8:%d a:%d b:%d c:%d\n', condition_8, ...
            vote_a, vote_b, vote_c )];
        % lfp_log opens and closes the log file, so it should not be called
        % often:
        if length(logmsg) > 20
            lfp_log(logmsg);
        end
        % Try again with the next candidate blink event:
        canblinksidx = canblinksidx + 1;
    end
    waitbar(canblinksidx/length(canblinks), h);
end
close(h)

if length(blinkstarts) ~= length(blinkends)
    error('lfp_parseEyeTraces6:oops1', ...
        'Unbalanced blink events');
end
if isempty(blinkstarts)
    warning('lfp_parseEyeTraces6:oops2', ...
        'There are no blink events!');
    lfp_log('There are no blink events!');
end
lfp_log('lfp_parseEyeTraces6 finding saccades');
% Find big saccades <bigsaccs> that are not too close to (or identical
% with) blinks:
rawevents = vpeaks;
eventgood = lfp_markGoodEvents4Dummies(rawevents, ...
    blinkstarts, blinkends, NBlinkTO, NRecTO);
bigsaccs = rawevents(eventgood);

% Find the ones that are in the {122 123} evtbounds for histogram analysis:
inbounds = false(size(bigsaccs));
for trial = 1:length(lfp_SelectedTrials)
    trialevts = lfp_Events(lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), :);
    eix122 = find(trialevts(:,2)==122) + lfp_TrialIndex(trial,1) - 1;
    eix123 = find(trialevts(:,2)==123) + lfp_TrialIndex(trial,1) - 1;
    if ~isempty(eix123) && ~isempty(eix123)
        inbounds( ...
            (bigsaccs >= lfp_time2index(lfp_Events(eix122,1))) & ...
            (bigsaccs <= lfp_time2index(lfp_Events(eix123,1))) ) = true;
    else
        warning('lfp_parseEyeTraces6:evtbounds', ...
            'Evt 122 or 123 is missing from trial %d', trial );
    end
end      
boundedbigsaccs = bigsaccs(inbounds);

% Find runs of close-interval saccades and compute intra-run peak v ratios
% for histogram analysis:
maxscale = 9;
minsaccint = .05;   % shortest interval that will be considered separate
Nminsaccint = round(minsaccint/lfp_SamplePeriod);
% The "unique" is probably not necessary, but "better safe":
saccs2plot = unique(boundedbigsaccs);
saccsdif = diff(saccs2plot);
is_short = saccsdif < Nminsaccint;
% First v peak in run:
if isempty(is_short)
    runstartidx = [];
else
    runstartidx = find([ is_short(1)
        is_short(2:end) & ~is_short(1:end-1) ]);
end
if length(runstartidx) < maxscale * 100
    logmsg = 'Not enough data for reliable v-peak ratios histo in {122 123}, trying whole trials';
    lfp_log(['lfp_parseEyeTraces6: ' logmsg]);
    warning('lfp_parseEyeTraces6:dubious3', '%s', logmsg);
    saccs2plot = unique(bigsaccs);
    saccsdif = diff(saccs2plot);
    is_short = saccsdif < Nminsaccint;
    % First v peak in run:
    if isempty(is_short)
        runstartidx = [];
    else
        runstartidx = find([ is_short(1)
            is_short(2:end) & ~is_short(1:end-1) ]);
    end
    if length(runstartidx) < maxscale * 100
        logmsg = 'Not enough data for reliable v-peak ratios in whole trials either';
        lfp_log(['lfp_parseEyeTraces6: ' logmsg]);
        warning('lfp_parseEyeTraces6:dubious4', '%s', logmsg);
    end
end
% <runendidx> is the last v-peak-short-interval-pair in run, i.e. it points to 
% the second-to-last v peak in the run:
runendidx = find([ ~is_short(2:end) & is_short(1:end-1)
    is_short(end) ]);
ratios = [];
ratios2 = [];
for runnum = 1:length(runstartidx)
    ratios(end+1,1) = lfp_Samples{velo}(saccs2plot(runstartidx(runnum)+1)) ...
        / lfp_Samples{velo}(saccs2plot(runstartidx(runnum)));
    if runstartidx(runnum) < runendidx(runnum)
        ratios2(end+1,1) = lfp_Samples{velo}(saccs2plot(runstartidx(runnum)+2)) ...
            / lfp_Samples{velo}(saccs2plot(runstartidx(runnum)));
    end
end

% Defaults when histo fails:
mergethresh = 0.6;
revtrimthresh = 2;

if isempty(runstartidx)
    logmsg = 'No short-interval saccade runs for ratio histogram; using defaults';
    lfp_log(['lfp_parseEyeTraces6: ' logmsg]);
    warning('lfp_parseEyeTraces6:ratiohisto1', '%s', logmsg);
else
    % Histogram with 0.1-wide or narrower bins
    nbins = max(maxscale * 10, round(length(runstartidx)/10));
    binedges = (0:nbins)*maxscale/nbins;
    binctrs = (binedges(1:end-1) + binedges(2:end))/2;
    rawhisto = reshape(histc(ratios, binedges), [], 1);
    if ~isempty(ratios2)
        rawhisto2 = reshape(histc(ratios2, binedges), [], 1);
    else
        rawhisto2 = zeros(length(binedges),1);
    end

    % Titrate to find right smoothing, defined as that which produces 3 or
    % fewer minima between the absolute max of histo and ratio = 1.
    [maxval, maxbin] = max(rawhisto);
    [m, onebin] = min(abs(binctrs-1));
    startbin = maxbin + 1;
    if maxbin > onebin
        logmsg = 'V peak ratio histo has max above 1';
        lfp_log(['lfp_parseEyeTraces6: ' logmsg]);
        warning('lfp_parseEyeTraces6:ratiohisto2', '%s', logmsg);
    end
    smoothfrac = 0.002;
    smoothfracstrs = {};
    k = 1;
    nummin = Inf;
    maxiters = ceil(log(1/smoothfrac)/log(1.1));
    while (nummin>3 && k<maxiters)
        smoothfracstrs{k,1} = sprintf('%.4f', smoothfrac);
        smoothpts = ceil(smoothfrac*nbins);
        smoothfunc = hanning(2*smoothpts+1);
        smoothed = conv(rawhisto, smoothfunc) / sum(smoothfunc);
        smoothedhisto(:,k) = smoothed(1+smoothpts : end-smoothpts);
        if ~isempty(ratios2)
            smoothed = conv(rawhisto2, smoothfunc) / sum(smoothfunc);
            smoothedhisto2(:,k) = smoothed(1+smoothpts : end-smoothpts);
        else
            smoothedhisto2(:,k) = zeros(length(binedges),1);
        end
        % Must use flat-bottommed min finder, not 3-point:
        isvalley = false(onebin - startbin + 1, 1);
        wasdecreasing = false;   % latch to remember prior decreasing behavior
        for binnum = (startbin+1) : (onebin-1)
            if wasdecreasing
                if smoothedhisto(binnum) < smoothedhisto(binnum+1)
                    % we have found a minimum, possibly with a flat bottom
                    isvalley(binnum - startbin + 1) = true;
                    wasdecreasing = false;
                end
            else
                wasdecreasing = smoothedhisto(binnum) < smoothedhisto(binnum-1);
            end
        end
        nummin = sum(isvalley);
        smoothfrac = smoothfrac * 1.1;
        k = k+1;
    end
    if k >= maxiters
        warning('lfp_parseEyeTraces6:maxiters', ...
            'There is no smoothing that produces the required number of minima in ratios histo; using default mergethresh & revtrimthresh.');
    else
        if nummin > 0
            minima = find(isvalley);
            mergethresh = binctrs(minima(1) + startbin - 1);
            almostnone = find(smoothedhisto((startbin+1):end, end) < 0.1);
        else
            warning('lfp_parseEyeTraces6:nomergethresh', ...
                'Could not find value for mergethresh in ratio histogram');
        end
        if nummin == 0 || isempty(almostnone)
            warning('lfp_parseEyeTraces6:norevtrimthresh', ...
                'Could not find value for revtrimthresh in ratio histogram');
        else
            revtrimthresh = max(binctrs(almostnone(1) + startbin - 1), ...
                3 * mergethresh);
        end
    end

    % plot histos
    hF = figure;
    hA = axes;
    plot(hA, binctrs, [rawhisto(1:end-1) smoothedhisto(1:end-1,:)]);
    legend(['raw'; smoothfracstrs]);
    hold on;
    hL = plot([mergethresh mergethresh], get(gca, 'YLim'),'r');
    set(hL, 'ButtonDownFcn', 'disp(''mergethresh'')');
    hL = plot([revtrimthresh revtrimthresh], get(gca, 'YLim'),'g');
    set(hL, 'ButtonDownFcn', 'disp(''revtrimthresh'')');
    ylabel(hA, 'Numbers of ratios');
    xlabel(hA, 'Peak2/Peak1 velocity ratio');
    if max(lfp_OrigTrialNums) == length(lfp_SelectedTrials)
        fragname = lfp_SessionNames{1};
    else
        [pathstr, name] = fileparts(lfp_DataDir);
        fragname = sprintf('%s %s', lfp_SessionNames{1}, name);
    end
    figname = sprintf('%s v-peak ratios', fragname);
    title(hA, figname, 'Interpreter', 'none');

    if ~isempty(ratios2)
        hF2 = figure;
        hA2 = axes;
        plot(hA2, binctrs, [rawhisto2(1:end-1) smoothedhisto2(1:end-1,:)]);
        legend(['raw'; smoothfracstrs]);
        hold on;
        hL = plot([mergethresh mergethresh], get(gca, 'YLim'),'r');
        set(hL, 'ButtonDownFcn', 'disp(''mergethresh'')');
        hL = plot([revtrimthresh revtrimthresh], get(gca, 'YLim'),'g');
        set(hL, 'ButtonDownFcn', 'disp(''revtrimthresh'')');
        ylabel(hA2, 'Numbers of ratios');
        xlabel(hA2, 'Peak3/Peak1 velocity ratio');
        if max(lfp_OrigTrialNums) == length(lfp_SelectedTrials)
            fragname = lfp_SessionNames{1};
        else
            [pathstr, name] = fileparts(lfp_DataDir);
            fragname = sprintf('%s %s', lfp_SessionNames{1}, name);
        end
        figname2 = sprintf('%s v-peak skip-one ratios', fragname);
        title(hA2, figname2, 'Interpreter', 'none');
    end

    saveas(hF, fullfile(lfp_DataDir, [figname '.fig']));
    % close(hF);
    if ~isempty(ratios2)
        saveas(hF2, fullfile(lfp_DataDir, [figname2 '.fig']));
        %     close(hF2);
    end
end

if ismember('vpeakratiomergethresh', fieldnames(threshes)) && ...
        ~isempty(threshes.vpeakratiomergethresh)
    msg = sprintf( ...
        'Overriding computed mergethresh=%d with supplied threshes.vpeakratiomergethresh=%d', ...
        mergethresh, threshes.vpeakratiomergethresh );
    mergethresh = threshes.vpeakratiomergethresh;
    warning('lfp_parseEyeTraces6:override1', '%s', msg);
    lfp_log(msg);
end
if ismember('vpeakratiorevtrimthresh', fieldnames(threshes)) && ...
        ~isempty(threshes.vpeakratiorevtrimthresh)
    msg = sprintf( ...
        'Overriding computed revtrimthresh=%d with supplied threshes.vpeakratiorevtrimthresh=%d', ...
        revtrimthresh, threshes.vpeakratiorevtrimthresh );
    revtrimthresh = threshes.vpeakratiorevtrimthresh;
    warning('lfp_parseEyeTraces6:override2', '%s', msg);
    lfp_log(msg);
end
lfp_log(sprintf( ...
    'lfp_parseEyeTraces6 starting saccade grouping, mergethresh=%d revtrimthresh=%d minsaccint=%d', ...
    mergethresh, revtrimthresh, minsaccint ));

% Find isolated saccades
% <shortintidx>, <runstartidx>, <runendidx>, <isolatedsaccidx> are indices
% into bigsaccs; members of <isolatedsaccidx> are neither the first of a
% pair with a short interval, nor the second (hence the setdiff) - i.e.
% they are truly isolated, not close to any other member of bigsaccs.
saccsdif = diff(bigsaccs);
is_short = saccsdif < Nminsaccint;
shortintidx = find(is_short);
isolatedsaccidx = setdiff(find(~is_short), shortintidx + 1); 

% Group the complex saccades.
% The intention here is that "split" saccades at the beginning of a run
% should be merged as long as their v peak ratio is >= mergethresh, but
% starting with the first pair whose ratio is < mergethresh, all the rest
% of the saccades in the run should be considered corrective and be trimmed
% away.  
%   There may also be pairs that need to be "reverse trimmed", i.e. the
% first saccade in the pair should be deleted when the ratio is >=
% revtrimthresh. For the purpose of doing the subsequent curve fitting, it
% is sufficient to have a list of isolated single saccades to fit, and a
% separate list of pairs of starting and ending saccades of a run that is
% to be merged.
%   Last (but actually handled first in the code), a special category of
% complex saccade is identified by looking at the median velocity in a
% delayed time window after the first peak; if the median velocity is above
% threshes.blinkv1, then the entire run is considered to be rhino/hedgehog.
%
% <run1idx2> indexes <runstartidx>, <runendidx>.  Its members are isolated pairs,
% which are handled the same as runs.  They should all end up in
% <whittledDownToOneidx2>, possibly among others.  None of them should be
% the same as <isolatedsaccidx>.  <run1idx2> is not used computationally,
% but is calculated for debugging purposes.
if isempty(is_short)
    runstartidx = [];
else
    runstartidx = find([ is_short(1)
        is_short(2:end) & ~is_short(1:end-1) ]);
end
runendidx = find([ ~is_short(2:end) & is_short(1:end-1)
    is_short(end) ]);
runlength = runendidx - runstartidx + 1;
run1idx2 = find(runlength==1);
if ~isempty(runstartidx) && any(ismember(runstartidx(run1idx2), isolatedsaccidx))
    error('lfp_parseEyeTraces6:bug2', ...
        'Something impossible happened.  Please despair.');
end
mergeidx = [];  % col 1: starting saccidx; col 2: ending sacc idx
whittledDownToOneidx2 = [];
rhinohedgeTS = [];
rhinosampleoffset = ...
    round(.01/lfp_SamplePeriod) : round(.04/lfp_SamplePeriod);
for runnum = 1:length(runstartidx)
    % Compute ratios of 2nd through "last+1st" peak of the run to the first
    % peak of the run:
    ratios = lfp_Samples{velo}( ...
        bigsaccs( (runstartidx(runnum):runendidx(runnum)) + 1 ) ) ...
        / lfp_Samples{velo}(bigsaccs(runstartidx(runnum)));
    % <relmergestart> is an index into the current run; strictly speaking,
    % it should be called relmergestartidx but that's just getting too
    % long, plus it's not an absolute index but only relative to
    % runstartidx(runnum).
    relmergestart = 1; 
    % The saccades to merge start with the first pair after those that need
    % to be reverse trimmed, and end with the last pair before those that
    % need to be trimmed.  However, since ratios are relative to the first
    % saccade, each saccade that gets reverse trimmed makes it necessary to
    % recalculate <ratios>:
    while relmergestart < (runendidx(runnum) - runstartidx(runnum) + 1) ...
            && ratios(1) >= revtrimthresh
        relmergestart = relmergestart + 1;
        ratios = lfp_Samples{velo}( bigsaccs( ...
            ((runstartidx(runnum)+relmergestart-1):runendidx(runnum)) + 1 )) ...
            / lfp_Samples{velo}(bigsaccs(runstartidx(runnum)+relmergestart-1));
    end
    % <relmergestart> now points at first sacc to merge.
    % <rhinowindow> is 10ms to 40ms after the peak:
    rhinowindow = min( Nsamples, bigsaccs(relmergestart + runstartidx(runnum) - 1) ...
        + rhinosampleoffset );
    if median(lfp_Samples{velo}(rhinowindow)) > threshes.blinkv1
        % Rhinoceros/hedgehog detected; keep it intact, except if it has a
        % tail of one or two peaks that fail the ratio test, AND the tail
        % is preceded by a minimum that falls below maxfixv, then dock 
        % them:
        rhinohedgeTS(end+1,1) = ...
            lfp_index2time(bigsaccs(relmergestart + runstartidx(runnum) - 1));
        mergeendidx = runendidx(runnum) + 1;
        if ratios(end) < mergethresh
            % find previous minimum:
            prevmin = bigsaccs(mergeendidx);
            while prevmin > 0 ...
                    && ( lfp_Samples{velo}(prevmin) >= ...
                    lfp_Samples{velo}(max(1, prevmin-1)) )
                prevmin = prevmin - 1;
            end
            if lfp_Samples{velo}(prevmin) < maxfixv
                % Last peak is a tail; check next to last
                if length(ratios) > 1 && ratios(end-1) < mergethresh
                    prevmin = bigsaccs(mergeendidx);
                    while prevmin > 0 ...
                            && ( lfp_Samples{velo}(prevmin) >= ...
                            lfp_Samples{velo}(max(1, prevmin-1)) )
                        prevmin = prevmin - 1;
                    end
                    if lfp_Samples{velo}(prevmin) < maxfixv
                        % 2nd to last is also part of tail:
                        mergeendidx = runendidx(runnum) - 1;
                    else
                        mergeendidx = runendidx(runnum);
                    end
                else
                    mergeendidx = runendidx(runnum);
                end
            end
        end
        mergestartidx = relmergestart + runstartidx(runnum) - 1;
        if mergestartidx == mergeendidx
            % We're down to a single saccade at <mergestartidx>
            whittledDownToOneidx2(end+1,1) = mergestartidx;
        else
            mergeidx(end+1, :) = [mergestartidx mergeendidx];
        end
    else
        % Apply the ratio criterion.  Note that <pairs2trimidx> points to
        % the FIRST peak of each SHORT-INTERVAL PEAK PAIR whose SECOND peak
        % should get trimmed.  Also note that <ratios> is now relative to
        % relmergestart, not to runstartidx(runnum).
        pairs2trimidx = find(ratios < mergethresh) + relmergestart - 1;
        if isempty(pairs2trimidx)
            % the last sacc in the run:
            relmergeend = runendidx(runnum) - runstartidx(runnum) + 2;
        else
            relmergeend = pairs2trimidx(1);
        end
        % <relmergeend> now points at the last sacc to merge
        if relmergestart > relmergeend
            warning('lfp_parseEyeTraces6:bogusrun', ...
                'BUG! Run from t = %.6f to t = %.6f is impossible.', ...
                lfp_index2time(bigsaccs(runstartidx(runnum))), ...
                lfp_index2time(bigsaccs(runendidx(runnum))) );
        elseif relmergestart == relmergeend
            % We're down to a single saccade at <relmergestart>
            whittledDownToOneidx2(end+1,1) = ...
                relmergestart + runstartidx(runnum) - 1;
        else
            mergeidx(end+1, :) = ...
                [relmergestart relmergeend] + runstartidx(runnum) - 1;
        end
    end
end
if isempty(rhinohedgeTS)
    lfp_log('No rhino/hedges found.');
else
    rhinofilename = fullfile(lfp_DataDir, 'rhinohedge.mat');
    save(rhinofilename, 'rhinohedgeTS');
    lfp_log(sprintf( ...
        '%d rhino/hedges saved to %s', ...
        length(rhinohedgeTS), rhinofilename ));
end

%
% Do the curve fitting to find saccstarts and saccends
%

saccstarts = zeros(size(bigsaccs));
saccends = zeros(size(bigsaccs));
btsaccstarts = zeros(size(bigsaccs));
btsaccends = zeros(size(bigsaccs));
cplxsaccstarts = zeros(size(bigsaccs));
cplxsaccends = zeros(size(bigsaccs));
fitstarts = zeros(size(bigsaccs));
fitends = zeros(size(bigsaccs));
% <singlesaccidx> is an index into bigsaccs (unique should be unnecessary,
% but it's nice to have it sorted): 
singlesaccidx = unique([isolatedsaccidx; whittledDownToOneidx2]);

% First the single saccades:
h = waitbar(0,'Finding saccade starts and ends...');
prevwaitbarval = 0;
totalwaitlength = (length(singlesaccidx) + size(mergeidx,1))*1.01;
% <rmses> cells contain the following.  1: isolated; 2-4: multipeaks 2-4;
% 5: multipeaks > 4; 6: whittled Down To One   
rmses = cell(1,6); 
for saccidx2 = 1:length(singlesaccidx)
    saccidx = singlesaccidx(saccidx2);
    [c, numpeaks, gof, category, startfit, endfit, term] = fitsaccade( ...
        velo, accel, bigsaccs, saccidx, threshes, blinkstarts);
    if isempty(c)
        warning('lfp_parseEyeTraces6:nofit1', ...
            'Fitting for saccade at t=%s was skipped.', ...
            mat2str(lfp_index2time(bigsaccs(saccidx))) );
        continue
    end
    if saccidx2 > length(isolatedsaccidx)
        rmses{6}(end+1) = gof.rmse;
    else
        rmses{1}(end+1) = gof.rmse;
    end
    [fittedwave, saccwaverange, saccstartpt, saccendpt] = ...
        multigauss(c, bigsaccs, saccidx, numel(saccwave), threshes);
    startpt = saccstartpt + saccwaverange(1) - 1;
    rawendpt = saccendpt + saccwaverange(1) - 1;
    adjendpt = adjustsaccend(startpt, ...
            blinkstarts, rawendpt, velo, Inf, ...
            round(0.8 * (rawendpt - startpt + 1)) );
    if term == 'b'
        btsaccstarts(saccidx) = startpt;
        btsaccends(saccidx) = adjendpt;
    elseif category == 1
        saccstarts(saccidx) = startpt;
        saccends(saccidx) = adjendpt;
    else
        cplxsaccstarts(saccidx(1)) = startpt;
        cplxsaccends(saccidx(1)) = adjendpt;
    end
    saccwave(saccwaverange) = saccwave(saccwaverange) + fittedwave;
    fitstarts(saccidx(1)) = startfit;
    fitends(saccidx(1)) = endfit;
    newwaitbarval = saccidx2/totalwaitlength;
    if newwaitbarval >= prevwaitbarval + .01
        waitbar(newwaitbarval, h);
        prevwaitbarval = newwaitbarval;
    end
end

% Then the merged runs.
for mergeidx2 = 1:size(mergeidx,1)
    saccidx = mergeidx(mergeidx2,:);
    c = zeros(12,1);
    c(3:3:length(c)) = 1; % these are divisors, so should not default to 0
    [coeffs, numpeaks, gof, category, startfit, endfit, term] = fitsaccade( ...
        velo, accel, bigsaccs, saccidx, threshes, blinkstarts);
    if numpeaks > 0 && numpeaks <=4
        rmses{numpeaks}(end+1) = gof.rmse;
    else
        rmses{5}(end+1) = gof.rmse;
    end
    if isempty(coeffs)
        warning('lfp_parseEyeTraces6:nofit2', ...
            'Fitting for saccade at t=%s was skipped.', ...
            mat2str(lfp_index2time(bigsaccs(saccidx))) );
        continue
    end
    c(1:length(coeffs)) = coeffs;
    % For multiple peaks, there is no simple analytic formula for saccstart
    % and saccend in this case, so they must be found from the fitted
    % waveform.
    [fittedwave, saccwaverange, saccstartpt, saccendpt] = ...
        multigauss(c, bigsaccs, saccidx, numel(saccwave), threshes);
    saccwave(saccwaverange) = saccwave(saccwaverange) + fittedwave;
    startpt = saccstartpt + saccwaverange(1) - 1;
    rawendpt = saccendpt + saccwaverange(1) - 1;
    if term == 'b'
        btsaccstarts(saccidx(1)) = startpt;
        btsaccends(saccidx(1)) = adjustsaccend(startpt, ...
            blinkstarts, rawendpt, velo, Inf, ...
            round(0.8 * (rawendpt - startpt + 1)) );
    elseif category == 1
        error('lfp_parseEyeTraces6:bug', ...
            'Something impossible has occurred.  Please panic.');
    else
        cplxsaccstarts(saccidx(1)) = startpt;
        cplxsaccends(saccidx(1)) = adjustsaccend(startpt, ...
            blinkstarts, rawendpt, velo, maxfixv, ...
            round(0.8 * (rawendpt - startpt + 1)) );
    end
    fitstarts(saccidx(1)) = startfit;
    fitends(saccidx(1)) = endfit;
    newwaitbarval = (mergeidx2+length(singlesaccidx)) ...
        / totalwaitlength;
    if newwaitbarval >= prevwaitbarval + .01
        waitbar(newwaitbarval, h);
        prevwaitbarval = newwaitbarval;
    end
end
% Remove the pre-allocated empties:
fitstarts(fitstarts==0) = [];
saccstarts(saccstarts==0) = [];
btsaccstarts(btsaccstarts==0) = [];
cplxsaccstarts(cplxsaccstarts==0) = [];
waitbar((mergeidx2 + length(singlesaccidx) + .005*totalwaitlength) ...
        / totalwaitlength );
fitends(fitends==0) = [];
saccends(saccends==0) = [];
btsaccends(btsaccends==0) = [];
cplxsaccends(cplxsaccends==0) = [];
close(h)

if length(saccstarts) ~= length(saccends)
    error('lfp_parseEyeTraces6:oops2', ...
        'Unbalanced saccade events');
end

if any(saccends <= saccstarts)
    isbadsacc = saccends <= saccstarts;
    warning('lfp_parseEyeTraces6:disorderedsacc', ...
        'Saccades have end at or before start at these sample indices:\n%s', ...
        dg_thing2str(saccstarts(isbadsacc)) );
end

% Remove any simple saccs that span < 30 pixels in space.
saccspans = sqrt( ...
    (lfp_Samples{eyeX}(saccends) - lfp_Samples{eyeX}(saccstarts)).^2 ...
    + (lfp_Samples{eyeY}(saccends) - lfp_Samples{eyeY}(saccstarts)).^2 );
isbadsacc = saccspans < 30;
saccstarts(isbadsacc) = [];
saccends(isbadsacc) = [];
    
% Search for complex saccades that overlap simple saccades.  We assume that
% the event lists are sorted in temporally increasing order.  Assuming that
% an "end" is always strictly greater than its matching "start", "overlap" 
% implies one of the following cases:
%   1. saccstart < cplxsaccstart && saccend >= cplxsaccstart
%   2. saccstart >= cplxsaccstart && saccstart <= cplxsaccend
% In the same loop, remove any cplxsaccs that span < 30 pixels in space.
saccidx = 1;
isbadsacc = false(size(saccstarts));
for cplxsaccidx = 1:length(cplxsaccstarts)
    while saccends(saccidx) < cplxsaccstarts(cplxsaccidx) ...
        && saccidx < length(saccends)
        saccidx = saccidx + 1;
    end
    if saccends(saccidx) < cplxsaccstarts(cplxsaccidx)
        % there are no more saccs to check
        break
    end
    % we now have saccends(saccidx) >= cplxsaccstarts(cplxsaccidx)
    while saccidx < length(saccends) && ...
            saccstarts(saccidx) < cplxsaccstarts(cplxsaccidx) ...
            && saccends(saccidx) >= cplxsaccstarts(cplxsaccidx)
        % overlap condition 1
        isbadsacc(saccidx) = true;
        saccidx = saccidx + 1;
    end
    while saccidx < length(saccends) && ...
            saccstarts(saccidx) >= cplxsaccstarts(cplxsaccidx) ...
            && saccstarts(saccidx) <= cplxsaccends(cplxsaccidx)
        % overlap condition 2
        isbadsacc(saccidx) = true;
        saccidx = saccidx + 1;
    end
end
% Remove overlapping simple saccades
if any(isbadsacc)
    warning('lfp_parseEyeTraces6:overlap', ...
        'Removing %d overlapped simple saccades', ...
        sum((isbadsacc)) );
    saccstarts(isbadsacc) = [];
    saccends(isbadsacc) = [];
end

% Construct fixations by maintaining three state variables: inblink,
% inrecTO, insacc, while running through a time-sorted list of just those
% events.  A fixation starts when ALL variables become false, and ends when
% ANY variable becomes true.  The event codes in this strictly local event
% list are > 0 for event starts, < 0 for event ends, and magnitude 1, 2, 3
% for inblink, inrecTO, insacc respectively.  Keep identifiable post-blink
% fixations separate from other fixations, using state variable wasblink
% for this purpose (wasblink remains true until the end of the fixation).
% There is also an event 4 that marks the actual start and end of the
% recording breaks themselves, for the purpose of clearing the state
% variables.

% First we must find functional blink boundaries, which for the purpose of
% finding fixations means the last point moving out from the blink where
% velocity is above threshes.sacc; those boundaries are denoted
% <blinkstarts2>, <blinkends2>.
blinkstarts2 = blinkstarts;
for k = 1:length(blinkstarts2)
    while blinkstarts2(k) > 0 && lfp_Samples{velo}(blinkstarts2(k)) > threshes.sacc
        blinkstarts2(k) = blinkstarts2(k) - 1;
    end
    blinkstarts2(k) = blinkstarts2(k) + 1;
end
blinkends2 = blinkends;
for k = 1:length(blinkends2)
    while blinkends2(k) <= numel(lfp_Samples{velo}) ...
            && lfp_Samples{velo}(blinkends2(k)) > threshes.sacc
        blinkends2(k) = blinkends2(k) + 1;
    end
    blinkends2(k) = blinkends2(k) - 1;
end

% Ditto for functional saccade boundaries
saccstarts2 = [btsaccstarts; saccstarts; cplxsaccstarts];
for k = 1:length(saccstarts2)
    rawsample = saccstarts2(k);
    while saccstarts2(k) > 0 && lfp_Samples{velo}(saccstarts2(k)) > threshes.sacc
        saccstarts2(k) = saccstarts2(k) - 1;
    end
    if saccstarts2(k) ~= rawsample
        % If while loop ran at all, then we have gone one sample too far:
        saccstarts2(k) = saccstarts2(k) + 1;
    end
end
saccends2 = [btsaccends; saccends; cplxsaccends];
for k = 1:length(saccends2)
    rawsample = saccends2(k);
    while saccends2(k) <= numel(lfp_Samples{velo}) ...
            && lfp_Samples{velo}(saccends2(k)) > threshes.sacc
        saccends2(k) = saccends2(k) + 1;
    end
    if saccends2(k) ~= rawsample
        % If while loop ran at all, then we have gone one sample too far:
        saccends2(k) = saccends2(k) - 1;
    end
end

events = [
    blinkstarts2, repmat(1, size(blinkstarts))
    blinkends2, repmat(-1, size(blinkends))
    lfp_RecSegments(:,2) - NRecTO + 1, repmat(2, size(lfp_RecSegments(:,1)))
    lfp_RecSegments(:,1) + NRecTO - 1, repmat(-2, size(lfp_RecSegments(:,1)))
    lfp_RecSegments(:,2), repmat(4, size(lfp_RecSegments(:,1)))
    lfp_RecSegments(:,1), repmat(-4, size(lfp_RecSegments(:,1)))
    saccstarts2, repmat(3, size(saccstarts2))
    saccends2, repmat(-3, size(saccends2))
    ];
events = sortrows(events, 1);
isrepeat = events(1:end-1, 2) == events(2:end, 2);
if any(isrepeat)
    warning('lfp_parseEyeTraces6:repeatevts', ...
        'Deleting %d repeated events in "events"', ...
        sum(isrepeat) );
    repnum = find(isrepeat);
    repIDs = unique(abs(events(repnum, 2)));
    isevt2delete = false(size(events,1), 1);
    for id = reshape(repIDs, 1, [])
        isrepstart = events(repnum, 2) == id;
        isrepend = events(repnum, 2) == -id;
        if sum(isrepstart) ~= sum(isrepend)
            warning('lfp_parseEyeTraces6:unbalrep', ...
                'There are unequal numbers of repeated starts and stops of class %d', ...
                id );
        end
        % merge the repeats by eliminating the second starts and the first
        % ends
        isevt2delete(repnum(isrepstart) + 1) = true;
        isevt2delete(repnum(isrepend)) = true;
    end
    % prep for log
    evt2deletelist = events(isevt2delete, :);
    evt2deletelist(:,1) = lfp_index2time(evt2deletelist(:,1));
    lfp_log( sprintf('Deleting %d repeated events in "events":\nState codes: >0 start, <0 end, 1 inblink, 2 inrecTO, 3 insacc\n%s', ...
        sum(isrepeat), mat2str(evt2deletelist)) );
    events(isevt2delete, :) = [];
end
inblink = false;
inrecTO = true;
insacc = false;
infix = false;
blinkfixstarts = zeros(size(events,1),1);
blinkfixends = zeros(size(events,1),1);
blinkfixdiams = NaN(size(events,1),1);
blinkfixpeakv = NaN(size(events,1),1);
otherfixstarts = zeros(size(events,1),1);
otherfixends = zeros(size(events,1),1);
otherfixdiams = NaN(size(events,1),1);
otherfixpeakv = NaN(size(events,1),1);
wasblink = false;
fixnum = 1;
for k = 1:size(events,1)
    % Update state variables, but only the leading edge of wasblink - it
    % must remain true until after fixation processing is entirely done.
    switch abs(events(k,2))
        case 1
            inblink = events(k,2) > 0;
            if ~inblink
                % This is where wasblink goes true, only if there is no
                % saccade in progress: 
                if insacc
                    wasblink = false;
                else
                    wasblink = true;
                    % infix can be already set if there was a sacc end
                    % after a rec break, but now we know we it wasn't
                    % really a legitimate fix.  We have faith that when the
                    % true fix start is found, it will overwrite the
                    % previously written value:
                    infix = false;
                end
            end
        case 2
            inrecTO = events(k,2) > 0;
        case 3
            insacc = events(k,2) > 0;
    end
    if infix
        if any([inblink inrecTO insacc])
            % found fix end; compute the "diameter" of the fixation:
            if wasblink
                blinkfixends(fixnum) = events(k,1) - 1;
                if blinkfixstarts(fixnum) == 0 ...
                        && otherfixstarts(fixnum) ~= 0
                    % Inconsistency handling (just to CMA):
                    msg = sprintf( ...
                        'Bug! Post-blink fix start %.6f in wrong list', ...
                        lfp_index2time(otherfixstarts(fixnum)) );
                    warning('lfp_parseEyeTraces6:missedblinkfix', ...
                        '%s', msg );
                    lfp_log(msg);
                    blinkfixstarts(fixnum) = otherfixstarts(fixnum);
                    otherfixstarts(fixnum) = 0;
                end
                if blinkfixends(fixnum) < blinkfixstarts(fixnum) + NMinFixDur
                    % Too short, back it out; note that it IS logically
                    % possible for blinkfixstarts(fixnum) to be greater than
                    % blinkfixends(fixnum):
                    blinkfixstarts(fixnum) = 0;
                    blinkfixends(fixnum) = 0;
                else
                    blinkfixdiams(k) = max(abs( ...
                        max(lfp_Samples{eyeX}( ...
                        blinkfixstarts(fixnum):blinkfixends(fixnum) )) ...
                        - min(lfp_Samples{eyeX}( ...
                        blinkfixstarts(fixnum):blinkfixends(fixnum) )) ...
                        ), abs( ...
                        max(lfp_Samples{eyeY}( ...
                        blinkfixstarts(fixnum):blinkfixends(fixnum) )) ...
                        - min(lfp_Samples{eyeY}( ...
                        blinkfixstarts(fixnum):blinkfixends(fixnum) )) ...
                        ));
                    blinkfixpeakv(k) = max(lfp_Samples{velo}( ...
                        blinkfixstarts(fixnum):blinkfixends(fixnum) ));
                    fixnum = fixnum + 1;
                end
            else
                if events(k,1) > 1
                    otherfixends(fixnum) = events(k,1) - 1;
                end
                if otherfixstarts(fixnum) == 0 ...
                        && blinkfixstarts(fixnum) ~= 0
                    % Inconsistency handling (just to CMA):
                    msg = sprintf( ...
                        'Bug! "Other" fix start %.6f in wrong list', ...
                        lfp_index2time(blinkfixstarts(fixnum)) );
                    warning('lfp_parseEyeTraces6:missedotherfix', ...
                        '%s', msg );
                    lfp_log(msg);
                    otherfixstarts(fixnum) = blinkfixstarts(fixnum);
                    blinkfixstarts(fixnum) = 0;
                end
                if otherfixends(fixnum) < otherfixstarts(fixnum) + NMinFixDur
                    % Too short, back it out; note that it IS logically
                    % possible for otherfixstarts(fixnum) to be greater than
                    % otherfixends(fixnum):
                    otherfixstarts(fixnum) = 0;
                    otherfixends(fixnum) = 0;
                else
                    otherfixdiams(k) = max(abs( ...
                        max(lfp_Samples{eyeX}( ...
                        otherfixstarts(fixnum):otherfixends(fixnum) )) ...
                        - min(lfp_Samples{eyeX}( ...
                        otherfixstarts(fixnum):otherfixends(fixnum) )) ...
                        ), ...
                        abs( ...
                        max(lfp_Samples{eyeY}( ...
                        otherfixstarts(fixnum):otherfixends(fixnum) )) ...
                        - min(lfp_Samples{eyeY}( ...
                        otherfixstarts(fixnum):otherfixends(fixnum) )) ...
                        ));
                    otherfixpeakv(k) = max(lfp_Samples{velo}( ...
                        otherfixstarts(fixnum):otherfixends(fixnum) ));
                    fixnum = fixnum + 1;
                end
            end
            infix = false;
            % update trailing edge of wasblink now that fixation processing
            % is done:
            wasblink = false;
        end
    else
        if ~any([inblink inrecTO insacc])
            % Found fix start.  After a rec break, it's possible to already
            % have a "false" fix start that has no fix end.  That's no
            % issue if ~wasblink, but if wasblink, then the false fix start
            % is recorded in the wrong list (i.e. it's in otherfixstarts)
            % and will NOT be overwritten here. Consequently, it's also
            % necessary to clear the entry in otherfixstarts 
            infix = true;
            if events(k,1) < numel(lfp_Samples{velo})
                if wasblink
                    blinkfixstarts(fixnum) = events(k,1) + 1;
                    otherfixstarts(fixnum) = 0;
                else
                    otherfixstarts(fixnum) = events(k,1) + 1;
                    % "It vouldn't hoit":
                    blinkfixstarts(fixnum) = 0;
                end
            end
        end
    end
    % update trailing edge of wasblink now that fixation processing is
    % done:
    if events(k,2) == 3
        wasblink = false;
    elseif abs(events(k,2)) == 4
        % Recording breaks should reset everything, including backing out
        % any in-progress fixation.  This "should normally be" redundant,
        % but both starts and ends of record breaks invoke the reset.
        if wasblink
            blinkfixstarts(fixnum) = 0;
        else
            otherfixstarts(fixnum) = 0;
        end
        wasblink = false;
        inblink = false;
        inrecTO = false;
        insacc = false;
        infix = false;
        % recTO ends should also start a new fixation.  It "should be"
        % impossible for wasblink to be true at this point, but it's a
        % crash if it turns out that there's weird junk inserted into the
        % rec break that set wasblink, so we have to allow for the
        % possibility. 
        if events(k,2) == -2
            infix = true;
            if events(k,1) < numel(lfp_Samples{velo})
                if wasblink
                    blinkfixstarts(fixnum) = events(k,1) + 1;
                else
                    otherfixstarts(fixnum) = events(k,1) + 1;
                end
            end
        end
    end
end
blinkfixstarts(blinkfixstarts==0) = [];
blinkfixends(blinkfixends==0) = [];
blinkfixdiams(isnan(blinkfixdiams)) = [];
blinkfixpeakv(isnan(blinkfixpeakv)) = [];
otherfixstarts(otherfixstarts==0) = [];
otherfixends(otherfixends==0) = [];
otherfixdiams(isnan(otherfixdiams)) = [];
otherfixpeakv(isnan(otherfixpeakv)) = [];
% Secret Easter egg feature:
save('blinkfixdiams', 'blinkfixdiams');
save('otherfixdiams', 'otherfixdiams');

% Find "normal" size for fixations if possible; if not default to 100
% pixels.
if length(otherfixdiams) < 50
    maxdiam = 100;
    lfp_log(sprintf( ...
        'Not enough data for fixdiam histo; using default maxdiam = %.2f', ...
        maxdiam));
else
    % Use 1-pixel bins to make life easy:
    binedges = floor(min(otherfixdiams)):ceil(max(otherfixdiams)*(1+eps));
    binctrs = binedges(1:end-1) + .5;
    rawhisto = histc(otherfixdiams, binedges);
    smoothpts = 10;
    smoothfunc = hanning(2*smoothpts+1);
    smoothed = conv(rawhisto, smoothfunc) / sum(smoothfunc);
    smoothedhisto = smoothed(1+smoothpts : end-smoothpts);
    [C,I] = max(smoothedhisto);
    maxbin = I(1);
    fitlimit = min(10*binctrs(maxbin), length(binctrs));
    p = 99;
    maxdiam = prctile(otherfixdiams, p);
    lfp_log(sprintf('maxdiam at %d percentile = %.2f', p, maxdiam));
    hF = figure;
    bins2plot = 1:min(fitlimit*10,length(binctrs));
    plot(binctrs(bins2plot), smoothedhisto(bins2plot));
    hold on;
    plot([maxdiam maxdiam], get(gca, 'YLim'), 'm');
    legend({'otherfixdiams' 'maxdiam'});
    figname = sprintf('%s fixation diameters', fragname);
    title(figname, 'Interpreter', 'none');
    saveas(hF, fullfile(lfp_DataDir, [figname '.fig']));
    close(hF);
end

[blinkfixstarts, blinkfixends, oldblinkfixstarts, oldblinkfixends, ...
    badblinkfixstarts, badblinkfixends, midblinkfix] = ...
    adjustfixations(blinkfixstarts, blinkfixends, blinkfixdiams, ...
    blinkfixpeakv, maxdiam, threshes.blinkv1, eyeX, eyeY, velo, threshes, NMinFixDur);
[otherfixstarts, otherfixends, oldotherfixstarts, oldotherfixends, ...
    badotherfixstarts, badotherfixends, midotherfix] = ...
    adjustfixations(otherfixstarts, otherfixends, otherfixdiams, ...
    otherfixpeakv, maxdiam, maxfixv, eyeX, eyeY, velo, threshes, NMinFixDur);
oldfixstarts = sort([oldblinkfixstarts; oldotherfixstarts]);
badfixstarts = sort([badblinkfixstarts; badotherfixstarts]);
oldfixends = sort([oldblinkfixends; oldotherfixends]);
badfixends = sort([badblinkfixends; badotherfixends]);
midfix = sort([midblinkfix; midotherfix]);

% Convert from sample indices to timestamps and construct result:
ts = {lfp_index2time(blinkstarts) lfp_index2time(blinkends) ...
    lfp_index2time(saccstarts) lfp_index2time(saccends) ...
    lfp_index2time(cplxsaccstarts) lfp_index2time(cplxsaccends) ...
    lfp_index2time(blinkfixstarts) lfp_index2time(blinkfixends) ...
    lfp_index2time(otherfixstarts) lfp_index2time(otherfixends) ...
    lfp_index2time(badfixstarts) lfp_index2time(badfixends) ...
    lfp_index2time(bigsaccs(mergeidx(:,1))) ...
    lfp_index2time(bigsaccs(mergeidx(:,2))) ...
    lfp_index2time(fitstarts) lfp_index2time(fitends) ...
    lfp_index2time(midfix) ...
    lfp_index2time(oldfixstarts) lfp_index2time(oldfixends) ...
    lfp_index2time(btsaccstarts) lfp_index2time(btsaccends) ...
    };
end


function [c, numpeaks, gof, category, startfit, endfit, term] = ...
    fitsaccade(velo, accel, bigsaccs, saccidx, threshes, blinkstarts)
%   If saccidx is a scalar, then one velo peak is fitted.
%   If saccidx is a two-element vector, then the entire interval extending
% from the first peak to the second (with any intervening peaks) is fitted
% unconditionally.
%   In either case, the start of the interval to be fitted is found by
% starting at the peaks given in saccidx, and moving leftwards until either
% the velocity falls below the fixed threshold threshes.sacc, or the
% velocity begins to increase again (i.e. a strict minimum is found).  The
% end is found similarly, except that minima only count if they are below
% blinkv1.  The additional restriction is applied that no data that are
% within blinkTO of the next blink start can be included.  If there's more
% than one velo peak internal to the range being merged, that by definition
% makes it a LNM (Loch Ness Monster), which means it has to be fitted all
% the way to the tip of its boat-smashing, terrorizing tail (i.e. minima
% don't count at all). 
%   <term> is the type of termination, 'b' for blink / end-of-data, 'f' for
% fixation 

% <category>:
%   1. simple ballistic saccade
%   2. single-peaked complex saccade
%   3. 2 to (length(funcnames)-1)-peaked complex saccade
%   4. length(funcnames)-peaked complex saccades
% <startfit>, <endfit>: sample numbers of time interval fitted
global lfp_Samples lfp_SamplePeriod
c = [];
term = 'f';
category = 0;
diffmaxch = 1;
diffminch = 1e-4;
tolfun = 1e-2;
tolx = 1e-3;
maxfunevals = 600;
maxgauss1width = .02;    % seconds
maxmultigausswidth = .01;    % seconds
NBlinkTO = round(threshes.blinkTO / lfp_SamplePeriod);
vmax = max(lfp_Samples{velo}(bigsaccs(saccidx(1)):bigsaccs(saccidx(end))));
startfit = bigsaccs(saccidx(1));
endfit = bigsaccs(saccidx(end));
nextblinkidx = find(blinkstarts>endfit);
if isempty(nextblinkidx)
    lastsample2fit = numel(lfp_Samples{velo}) - 1;
else
    lastsample2fit = blinkstarts(nextblinkidx(1)) - 1;
end
y = lfp_Samples{velo}(startfit:endfit);
numinternalpeaks = sum( ...
    (y(2:end-1) > y(1:end-2) & y(2:end-1) > y(3:end)) );
isLNM = numinternalpeaks > 1;
while startfit > 2 && lfp_Samples{velo}(startfit) > threshes.sacc ...
        && (lfp_Samples{velo}(startfit) >= lfp_Samples{velo}(startfit-1))
    startfit = startfit - 1;
end
% For LNM, search forward for <endfit> until v <= threshes.sacc; for
% non-LNM, search until either v <= threshes.sacc, or v is 
% increasing and <= threshes.blinkv1.
while endfit < lastsample2fit ...
        && lfp_Samples{velo}(endfit) > threshes.sacc ...
        && ( isLNM || ...
        (lfp_Samples{velo}(endfit) >= lfp_Samples{velo}(endfit+1)) ...
        || (lfp_Samples{velo}(endfit) > threshes.blinkv1) )
    endfit = endfit + 1;
end
if endfit == lastsample2fit
    % No suitable endpoint was found; this means we have a blink-terminated
    % or end-of-data-terminated saccade.  Look for the last v minimum
    % before lastsample2fit.
    lastminidx = lastsample2fit - 1;
    while lastminidx > bigsaccs(saccidx(1)) && ...
            lfp_Samples{velo}(lastminidx) >= lfp_Samples{velo}(lastminidx+1)
        % slope of v is negative or flat
        lastminidx = lastminidx - 1;
    end
    if lastminidx == bigsaccs(saccidx(1))
        % v is monotonic nonincreasing over the whole search interval;
        % leave <endfit> pointing at <lastsample2fit>, i.e. do nothing.
    else
        % found an increasing region; continue searching until the
        % preceding decreasing region is found to define the minimum (note:
        % this strategy can find minima with flat bottoms)
        while lastminidx > bigsaccs(saccidx(1)) && ...
                lfp_Samples{velo}(lastminidx) <= lfp_Samples{velo}(lastminidx+1)
            % slope of v is positive or flat
            lastminidx = lastminidx - 1;
        end
        if lastminidx == bigsaccs(saccidx(1))
            error('lfp_paresEyeTraces6:bug', ...
                'The supposed v peak at %.6f is not actually a peak', ...
                lfp_index2time(bigsaccs(saccidx(1))) );
        end
        newendfit = lastminidx + 1;
        msg = sprintf( ...
            ['End of interval to fit as saccade was truncated by blink or ' ...
            'end of data\nat %.6f s; using last v min at %.6f instead'], ...
            lfp_index2time(endfit), lfp_index2time(newendfit) );
        warning('lfp_parseEyeTraces6:endfit', '%s', msg);
        lfp_log(msg);
        endfit = newendfit(end);
        term = 'b';
    end
end
% bigsaccs(saccidx(1)) is arbitrarily used as a reference because it's
% conveniently available:
x = (startfit:endfit) - bigsaccs(saccidx(1));
% must normalize because 'fit' can't cope with radically different scales
% along different dimensions of the param search space:
v = lfp_Samples{velo}(startfit:endfit)/vmax;
% Find v peaks to use as starting values for the fit:
if length(saccidx) == 1
    [peakvals, peakposns] = max(v);
    numpeaks = 1;
else
    ispeak = [false ...
        (v(2:end-1) > v(1:end-2) & v(2:end-1) > v(3:end)) ...
        false];
    peakvals = v(ispeak);
    peakposns = x(ispeak);
    numpeaks = sum(ispeak);
end
funcnames = {
    'gauss1'
    'gauss2'
    'gauss3'
    'gauss4'
    'gauss5'
    'gauss6'
    'gauss7'
    'gauss8'
    };
switch numpeaks
    case 1
        a = lfp_Samples{accel}(startfit:endfit);
        ispeak = [false ...
            (a(2:end-1) > a(1:end-2) & a(2:end-1) > a(3:end)) ...
            false];
        accelpeakvals = v(ispeak);
        accelpeakposns = x(ispeak);
        numaccelpeaks = sum(ispeak);
        if numaccelpeaks == 1
            category = 1;
            fitfunc = 'gauss1';
        else
            category = 2;
            fitfunc = funcnames{max( ...
                min(numaccelpeaks, length(funcnames)), 1 )};
            % accel peaks make better starting values than random, although
            % they will generally shift to the right:
            if numaccelpeaks > length(funcnames)
                [accelpeakvals, ix] = sort(accelpeakvals, 'descend');
                peakvals = accelpeakvals(1:length(funcnames));
                peakposns = accelpeakposns(ix(1:length(funcnames)));
            else
                peakvals = accelpeakvals;
                peakposns = accelpeakposns;
            end
        end
    case mat2cell(2:(length(funcnames)-1),1,ones(size(2:(length(funcnames)-1))))
        category = 3;
        fitfunc = funcnames{numpeaks};
    otherwise
        % well, you've got to stop adding peaks sometime!
        category = 4;
        fitfunc = funcnames{end};
        % Use the highest-v peaks as starting values:
        [peakvals, ix] = sort(peakvals, 'descend');
        peakvals = peakvals(1:length(funcnames));
        peakposns = peakposns(ix(1:length(funcnames)));
end
% 'fit' gets upset sometimes when upper bounds are Inf, so we compute
% finite upper bounds instead:
almostinfv = max(v);
almostinfx = max(x);
if isequal(fitfunc, 'gauss1')
    [cfun,gof,output] = fit(x', v', 'gauss1', ...
        fitoptions( 'Method','NonlinearLeastSquares', ...
        'DiffMaxChange', diffmaxch, 'DiffMinChange', diffminch, ...
        'TolFun', tolfun, 'TolX', tolx', ...
        'MaxFunEvals', maxfunevals, ...
        'Lower', repmat([0 x(1) 1], 1, length(peakvals)), ...
        'Upper', repmat([almostinfv almostinfx maxgauss1width/lfp_SamplePeriod], ...
        1, length(peakvals)) ));
else
    % Alternative minimum number of peaks determined by duration to be
    % fitted; peaks after the first are allowed to be up to 2 *
    % maxmultigausswidth wide:
    minnumpeaks = round( ...
        (length(x)*lfp_SamplePeriod - maxmultigausswidth) ...
        / (2 * maxmultigausswidth) ) + 1;
    if minnumpeaks > numpeaks
        % But numpeaks must also not exceed length(funcnames):
        numpeaks = min(minnumpeaks, length(funcnames));
        % first calculate peak position indices (i.e. indices into x):
        firstposn = round(length(x)/(2*numpeaks));
        lastposn = round(length(x) * (1 - 1/(2*numpeaks)));
        peakposns = round( ...
            linspace(firstposn, lastposn, numpeaks) );
        peakvals = v(peakposns);
        % convert peakposns to x scale:
        peakposns = peakposns + x(1) - 1;
        if numpeaks < length(funcnames)
            fitfunc = funcnames{numpeaks};
        else
            category = 4;
            fitfunc = funcnames{end};
        end
    end
    [cfun,gof,output] = fit(x', v', fitfunc, ...
        fitoptions( 'Method','NonlinearLeastSquares', ...
        'DiffMaxChange', diffmaxch, 'DiffMinChange', diffminch, ...
        'TolFun', tolfun, 'TolX', tolx', ...
        'MaxFunEvals', maxfunevals, ...
        'StartPoint', reshape([ peakvals
        peakposns
        ones(size(peakvals)) ], 1, []), ...
        'Lower', repmat([0 x(1) 1], 1, length(peakvals)), ...
        'Upper', [ almostinfv almostinfx maxmultigausswidth/lfp_SamplePeriod ...
        repmat([almostinfv almostinfx 2*maxmultigausswidth/lfp_SamplePeriod ], ...
        1, length(peakvals) - 2 ) ...
        almostinfv almostinfx maxmultigausswidth/lfp_SamplePeriod ] ));
end
if gof.rmse > .12
    lfp_log(sprintf( ...
        'High RMS error %.3f for %d-peak fit at %.6f to %.6f s', ...
        gof.rmse, length(peakvals), ...
        lfp_index2time(startfit), lfp_index2time(endfit) ));
end
if output.funcCount >= maxfunevals
    lfp_log(sprintf( ...
        'Gave up after %d func evals in fit at %.6f to %.6f s', ...
        output.funcCount, ...
        lfp_index2time(startfit), lfp_index2time(endfit) ));
end
c = coeffvalues(cfun);
% De-normalize the vertical scale params:
for k = 1:length(c)/3
    c((k-1)*3 + 1) = c((k-1)*3 + 1) * vmax;
end

end


function [fittedwave, saccwaverange, saccstartpt, saccendpt] = ...
    multigauss(c, bigsaccs, saccidx, maxsample, threshes)
% Create the fitted waveform as specified by c, and find the start and end
% points of the saccade.  smarkingwidth_merged/emarkingwidth_merged work
% the opposite way from smarkingwidth/emarkingwidth: they are multipliers
% for threshes.sacc to nudge that threshold up or down for marking start
% and end of saccade respectively, so increasing the values results in
% temporally tighter saccade markers.
smarkingwidth_merged = 1;
emarkingwidth_merged = 1;
if mod(length(c), 3)
    error('multigauss:badc', ...
        '<c> must contain a multiple of 3 elements');
end
numpeaks = length(c)/3;
for k = 1:numpeaks
    offset = (k - 1) * 3;
    possiblestarts(k,1) = round(c(2+offset) - 3*c(3+offset));
end
fitstartpt = min(possiblestarts);
for k = 1:numpeaks
    offset = (k - 1) * 3;
    possibleends(k,1) = round(c(2+offset) + 3*c(3+offset));
end
fitendpt = max(possibleends);
fitrange = fitstartpt:fitendpt;
saccwaverange = fitrange + bigsaccs(saccidx(1));
% Adjust fitrange and saccwaverange if necessary:
if saccwaverange(1) < 1
    fitrange = fitrange(saccwaverange>0);
    saccwaverange = saccwaverange(saccwaverange>0);
end
if saccwaverange(end) > maxsample
    fitrange = fitrange(saccwaverange<=maxsample);
    saccwaverange = saccwaverange(saccwaverange<=maxsample);
end
fittedwave = zeros(size(fitrange));
for k = 1:numpeaks
    offset = (k - 1) * 3;
    fittedwave = fittedwave + c(1+offset) * ...
        exp(-((fitrange - c(2+offset)) / c(3+offset)).^2);
end
saccstartpt = 1;
while(fittedwave(saccstartpt) < threshes.sacc*smarkingwidth_merged ...
        && saccstartpt < length(fittedwave))
    saccstartpt =saccstartpt + 1;
end
saccendpt = length(fittedwave);
while(fittedwave(saccendpt) < threshes.sacc*emarkingwidth_merged ...
        && saccendpt > 1)
    saccendpt = saccendpt - 1;
end
end


function saccend = adjustsaccend(startpt, ...
    blinkstarts, saccend, velo, maxendv, NmaxT)
% Move saccade end to the nearest v minimum that is at or below <maxendv>
% (moving later wins ties).  If there is no v min <= <maxendv> within
% <NmaxT> samples of <saccend>, use the v min closest to <saccend>.  If
% there is no v min at all within <NmaxT> samples, do nothing but complain.
% The search for minima at later times only is limited to the last sample
% before the entry in <blinkstarts> that follows <startpt>, in the same
% manner as <NmaxT>.
global lfp_Samples;
highminvpt = [];
blinkstartidx = find(blinkstarts > startpt);
if isempty(blinkstartidx)
    Nblinkstart = Inf;
else
    % Note that Nblinkstart may be negative:
    Nblinkstart = blinkstarts(blinkstartidx(1)) - saccend;
end
k = 0;
while true
    if k < Nblinkstart
        pt1 = max( ...
            min(saccend + k - 1, numel(lfp_Samples{velo})), 1 );
        pt2 = min(saccend + k, numel(lfp_Samples{velo}));
        pt3 = min(saccend + k + 1, numel(lfp_Samples{velo}));
        if lfp_Samples{velo}(pt2) == 0 || ...
                ( lfp_Samples{velo}(pt2) < lfp_Samples{velo}(pt1) ...
                && lfp_Samples{velo}(pt2) < lfp_Samples{velo}(pt3) )
            if lfp_Samples{velo}(pt2) <= maxendv
                saccend = pt2;
                break
            else
                if isempty(highminvpt)
                    highminvpt = pt2;
                end
            end
        end
    end
    if -k < Nblinkstart
        pt1 = max(saccend - k - 1, 1);
        pt2 = min(saccend - k, numel(lfp_Samples{velo}));
        pt3 = max( ...
            min(saccend - k + 1, numel(lfp_Samples{velo})), 1 );
        if lfp_Samples{velo}(pt2) == 0 || ...
                ( lfp_Samples{velo}(pt2) < lfp_Samples{velo}(pt1) ...
                && lfp_Samples{velo}(pt2) < lfp_Samples{velo}(pt3) )
            if lfp_Samples{velo}(pt2) <= maxendv
                saccend = pt2;
                break
            else
                if isempty(highminvpt)
                    highminvpt = pt2;
                end
            end
        end
    end
    k = k + 1;
    if k > NmaxT
        if ~isempty(highminvpt)
            saccend = highminvpt;
            break
        else
            msg = sprintf('Could not adjust sacc end at %.6f s', ...
                lfp_index2time(saccend) );
            warning('lfp_parseEyeTraces6:badsaccend', ...
                '%s', msg);
            lfp_log(msg);
            break
        end
    end
end
end

function [fixstarts, fixends, oldfixstarts, oldfixends, ...
    badfixstarts, badfixends, midfix] = ...
    adjustfixations(fixstarts, fixends, fixdiams, fixpeakv, ...
    bigfixthresh, maxfixv, eyeX, eyeY, velo, threshes, NMinFixDur)
% Find the fixations that are abnormally large or high velocity.  Mark them
% right in the middle as <midfix>.  <fixstarts>, <fixends> are sample
% numbers.  <fixdiams>, <fixpeakv> are vectors of same length as
% <fixstarts>, <fixends>.

lfp_declareGlobals;

% <bigfixidx> indexes <fixdiams>, <fixpeakv>, <fixstarts>, <fixends>.  
bigfixidx = find(fixdiams>bigfixthresh | fixpeakv>maxfixv);
midfix = round( ...
    (fixstarts(bigfixidx) + fixends(bigfixidx)) / 2 );

% Adjust big fixations until they're not big anymore.
fNTv = lfp_findNegThresh3(velo, threshes.sacc, 'sample');
badfixidx = []; % indexes <fixdiams>, <fixpeakv>, <fixstarts>, <fixends>
oldfixstarts = [];
oldfixends = [];
% <midfixidx2> DOES index <midfix> and <bigfixidx>
logmsg = sprintf('Cumulative comment:\n');
for midfixidx2 = 1:length(bigfixidx)
    fix2fix = bigfixidx(midfixidx2);
    if fixends(fix2fix) < fixstarts(fix2fix) + NMinFixDur
        logmsg = sprintf( ...
            '%sfixation at %.6f - %.6f s too short to adjust\n', ...
            logmsg, ...
            lfp_index2time(fixstarts(fix2fix)), ...
            lfp_index2time(fixends(fix2fix)) );
        badfixidx(end+1) = fix2fix;
        continue
    end
    % <samplerange> is the range to search for positions and velocities
    samplerange = fixstarts(fix2fix) : fixends(fix2fix);
    % Start by trimming off any weird high v peaks from the end.  The Inf
    % is just a way to pre-allocate an extra element.
    myvmaxidx = [ dg_findPeaks( lfp_Samples{velo}(samplerange), ...
        maxfixv ) + samplerange(1) - 1
        Inf ];
    myvmaxidx = myvmaxidx(myvmaxidx>midfix(midfixidx2));
    if lfp_Samples{velo}(samplerange(end)) > threshes.sacc
        % In this case, the v peak causing the problem might actually be
        % off the end of samplerange, so we must include the last point as
        % an honorary v peak:
        myvmaxidx(end) = samplerange(end);
    else
        myvmaxidx(end) = [];
    end
    if isempty(myvmaxidx)
        logmsg = sprintf( ...
            '%sFixation end at %.6f s cannot be adjusted\n', ...
            logmsg, ...
            lfp_index2time(fixends(fix2fix)) );
    else
        newfixend = myvmaxidx(1);
        while (newfixend > fixstarts(fix2fix) + NMinFixDur - 1) ...
                && lfp_Samples{velo}(newfixend) > threshes.sacc
            newfixend = newfixend - 1;
        end
        if newfixend >= fixstarts(fix2fix) + NMinFixDur
            % Good fix end adjustment
            oldfixends(end+1,1) = fixends(fix2fix);
            fixends(fix2fix) = newfixend;
            samplerange = fixstarts(fix2fix) : newfixend;
        else
            logmsg = sprintf( ...
                '%sEnd-adjusted fix at %.6f - %.6f s is too short\n', ...
                logmsg, ...
                lfp_index2time(fixstarts(fix2fix)), ...
                lfp_index2time(fixends(fix2fix)) );
            badfixidx(end+1) = fix2fix;
            continue
        end
    end
    % Now adjust the start.  Try every threshes.sacc negative crossing
    % that occurs between the big fix start and minfixdur before the big
    % fix end.  <bigfixidx> is an index into <fixstarts>, <fixends>.
    % <midfix> is a sample number.  The existing fixation start,
    % fixstarts(fix2fix), is a candidate because the end adjustment might
    % have fixed the fix already.
    candidates = [ fixstarts(fix2fix) ...
        fNTv( fNTv > samplerange(1) & ...
        fNTv < (samplerange(end) - NMinFixDur) ) ];
    adjusted = false;
    for candidx = 1:length(candidates)
        candiam = max(abs( ...
            max(lfp_Samples{eyeX}(candidates(candidx):fixends(fix2fix))) ...
            - min(lfp_Samples{eyeX}(candidates(candidx):fixends(fix2fix))) ...
            ), ...
            abs( ...
            max(lfp_Samples{eyeY}(candidates(candidx):fixends(fix2fix))) ...
            - min(lfp_Samples{eyeY}(candidates(candidx):fixends(fix2fix))) ...
            ));
        candmaxv = max(lfp_Samples{velo}( ...
            candidates(candidx) : fixends(fix2fix) ));
        if candiam <= bigfixthresh && candmaxv <= maxfixv
            % adjustment found; save new values
            oldfixstarts(end+1,1) = fixstarts(fix2fix);
            fixstarts(fix2fix) = candidates(candidx) + 1;
            adjusted = true;
            break
        end
    end
    if ~adjusted
        % Explain why not
        if candiam > bigfixthresh && candmaxv > maxfixv
            reasonstr = 'diameter too big and max v too high';
        elseif candiam > bigfixthresh
            reasonstr = 'diameter too big';
        elseif candmaxv > maxfixv
            reasonstr = 'max v too high';
        else
            reasonstr = 'no good reason; shoot the programmer';
        end
        logmsg = sprintf( ...
            '%sFixation at %.6f - %.6f s could not be adjusted, marked bad because %s\n', ...
            logmsg, ...
            lfp_index2time(fixstarts(fix2fix)), ...
            lfp_index2time(fixends(fix2fix)), ...
            reasonstr );
        badfixidx(end+1) = fix2fix;
    end
end
if length(logmsg) > 20
    lfp_log(logmsg);
end

badfixstarts = fixstarts(badfixidx);
fixstarts(badfixidx) = [];
badfixends = fixends(badfixidx);
fixends(badfixidx) = [];
end
