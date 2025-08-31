function [ts, saccwave] = lfp_parseEyeTraces5(eyeX, eyeY, velo, threshes)
%lfp_parseEyeTraces5 for use with lfp_createEvents and lfp_createWave via 
% adapters; marks eye movement related events.  
%ts = lfp_parseEyeTraces5(eyeX, eyeY, velo, threshes)
% eyeX, eyeY, velo: filenums of X, Y, and velocity traces.
% threshes: a structure containing the thresholds used. The fields in the
% structure are:
%   blinkTO - time out around each blink during which lfp_EyeTabulation
%       discards events as artifactual
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
% not corrective or micro-) saccades and are accordingly called <bigsaccs>.
% Before further processing, the <bigsaccs> list is scanned for runs of v
% peaks that occur within 50 ms (hard-coded value of <minsaccint>) of each
% other, and the runs are handled as follows.  The ratio of each peak v in
% the run to the first peak v in the run is computed, and if it is above
% <mergethresh>, then the peak is considered to be part of a "split"
% saccade together with the first peak, and they are fitted en masse;
% otherwise, the entire tail of the run from the first v peak that is below
% <mergethresh> onwards is ignored, effectively deleting those (micro-)
% saccades from the list of saccades to be fitted. If more than three peaks
% are to be fitted en masse, that by definition makes it a LNM (Loch Ness
% Monster), which means it has to be fit all the way to the tip of its
% boat-smashing, terrorizing tail.
% <mergethresh> is found by constructing a histogram of the ratios of peak
% v's described above, applying enough smoothing to produce no more than 3
% minima between the global maximum bin and 1, and then finding the first
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
% blinkTO after a blink end; similarly for rec starts and ends.  It also
% means that there are guaranteed to be equal numbers of sacc starts and
% sacc ends.
%
% Finally, fixations are considered to be "all the rest", which is actually
% not as simple as it sounds.  Practically speaking, a fixation can start:
%   1 sample after sacc end
%   1 blinkTO after blink end
%   1 recTO after rec start
% and a fixation can end:
%   1 sample before sacc start
%   1 blinkTO before blink start
%   1 recTO before rec end
%
% <ts> is a cell vector of timestamp column vectors, with one column for
% each of: 
%   Start Blink
%   End Blink
%   Start Saccade
%   End Saccade
%   Start Merged Saccade Series
%   End Merged Saccade Series
%   Fixation Starts
%   Fixation Ends
%   Big Fixations (marked at midpoint of fixation)
% <saccwave> is an array of sample data the same size as
% lfp_Samples{velo}.  
%

% NOTES
% Basically a leaner, meaner rendition of lfp_parseEyeTraces4 and
% cleanEyeEvents (the concept, not any special version) in one function,
% with a sprinkle of lfp_vhistos and td_findshortsacc for generally
% improved functioning.
% 10:37PM 30-Apr-2008 DG replaced single Gaussian with multi-Gaussian for
% merged saccade fitting.
% 5:49PM 07-May-2008 DG added detection of "big fixations"

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
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
NMinFixDur = round(threshes.minfixdur / lfp_SamplePeriod);

%
% Find blinks.  
%
NBlinkTimeout2 = round(NBlinkTO/2);

% find sample indices of candidate blink events <blinkidx>:
h = waitbar(0,'Finding candidate eye events...');
vmaxidx = lfp_findPeaks(velo, threshes.blinkv1);
waitbar(1/5, h);
fPTx = lfp_findPosThresh3(eyeX, threshes.blinkx, 'sample')';
waitbar(2/5, h);
fNTx = lfp_findNegThresh3(eyeX, threshes.blinkx, 'sample')';
waitbar(3/5, h);
fPTy = lfp_findPosThresh3(eyeY, threshes.blinky, 'sample')';
waitbar(4/5, h);
fNTy = lfp_findNegThresh3(eyeY, threshes.blinky, 'sample')';
waitbar(5/5, h);
blinkidx = ...
    unique( [
    vmaxidx(lfp_Samples{velo}(vmaxidx) > threshes.blinkv2)
    fPTx
    fNTx
    fPTy
    fNTy
    ] );
close(h);

% blinkidx and vmaxidx are indices into lfp_Samples.
% blinkidx2, etc., are indices into blinkidx and vmaxidx
% ("2nd order" indices); they are used to obviate the need to
% search all of blinkidx and vmaxidx on each iteration.
blinkstarts = [];
blinkends = [];
blinkidx2 = 1;  % index into blinkidx ("2nd order" blink index);
vmaxidx2 = 1;  % starting index into vmaxidx for votes a & b search
h = waitbar(0,'Finding blinks...');
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
    waitbar(blinkidx2/length(blinkidx), h);
end
close(h)

% Find big saccades that are not too close to (or identical with) blinks:
rawevents = vmaxidx;
lfp_markGoodEvents4;
bigsaccs = rawevents(eventgood);

% Find the ones that are in the {122 123} evtbounds.
inbounds = false(size(bigsaccs));
for trial = 1:length(lfp_SelectedTrials)
    trialevts = lfp_Events(lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), :);
    eix122 = find(trialevts(:,2)==122) + lfp_TrialIndex(trial,1) - 1;
    eix123 = find(trialevts(:,2)==123) + lfp_TrialIndex(trial,1) - 1;
    inbounds( ...
        (bigsaccs >= lfp_time2index(lfp_Events(eix122,1))) & ...
        (bigsaccs <= lfp_time2index(lfp_Events(eix123,1))) ) = true;
end      
boundedbigsaccs = bigsaccs(inbounds);

% Find runs of close-interval saccades and compute intra-run peak v ratios.
minsaccint = .05;   % shortest interval that will be considered separate
Nminsaccint = round(minsaccint/lfp_SamplePeriod);
% The "unique" is probably not necessary, but "better safe":
saccs2plot = unique(boundedbigsaccs);
saccsdif = diff(saccs2plot);
is_short = saccsdif < Nminsaccint;
runstart = find([ is_short(1)
    is_short(2:end) & ~is_short(1:end-1) ]);
runend = find([ ~is_short(2:end) & is_short(1:end-1)
    is_short(end) ]);
ratios = [];
ratios2 = [];
for runidx = 1:length(runstart)
    ratios(end+1,1) = lfp_Samples{velo}(saccs2plot(runstart(runidx)+1)) ...
        / lfp_Samples{velo}(saccs2plot(runstart(runidx)));
    if runstart(runidx) < runend(runidx)
        ratios2(end+1,1) = lfp_Samples{velo}(saccs2plot(runstart(runidx)+2)) ...
            / lfp_Samples{velo}(saccs2plot(runstart(runidx)));
    end
end

% Defaults when histo fails:
mergethresh = 0.6;
revtrimthresh = 2;

% Histogram with 0.1-wide or narrower bins
maxscale = 9;
if length(runstart) < maxscale * 100
    warning('lfp_parseEyeTraces5:dubious3', ...
        'Not enough data for reliable v-peak ratios histo');
end
nbins = max(maxscale * 10, round(length(runstart)/10));
binedges = (0:nbins)*maxscale/nbins;
binctrs = (binedges(1:end-1) + binedges(2:end))/2;
rawhisto = histc(ratios, binedges);
if ~isempty(ratios2)
    rawhisto2 = histc(ratios2, binedges);
else
    rawhisto2 = zeros(length(binedges),1);
end

% Titrate to find right smoothing, defined as that which produces 3 or
% fewer minima between the absolute max of histo and ratio = 1.
[maxval, maxbin] = max(rawhisto);
[m, onebin] = min(abs(binctrs-1));
startbin = maxbin + 1;
if maxbin > onebin
    warning('lfp_parseEyeTraces5:histo', ...
        'V peak ratio histo has max above 1' );
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
    isvalley = [false
        smoothedhisto((startbin+1):onebin-1,k) < smoothedhisto(startbin:onebin-2,k) ...
        & smoothedhisto((startbin+1):onebin-1,k) < smoothedhisto((startbin+2):onebin,k)
        false];
    nummin = sum(isvalley);
    smoothfrac = smoothfrac * 1.1;
    k = k+1;
end
if k >= maxiters
    warning('lfp_parseEyeTraces5:maxiters', ...
        'There is no smoothing that produces the required number of minima');
else
    if nummin > 0
        minima = find(isvalley);
        mergethresh = binctrs(minima(1) + startbin - 1);
        almostnone = find(smoothedhisto((startbin+1):end, end) < 0.1);
    end
    if nummin == 0 || isempty(almostnone)
        warning('lfp_parseEyeTraces5:norevtrimthresh', ...
            'Could not find value for revtrimthresh in histogram');
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
plot([mergethresh mergethresh], get(gca, 'YLim'));
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
    plot([mergethresh mergethresh], get(gca, 'YLim'));
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
close(hF);
if ~isempty(ratios2)
    saveas(hF2, fullfile(lfp_DataDir, [figname2 '.fig']));
    close(hF2);
end



% Edit the saccades.
% The intention here is that "split" saccades at the beginning of a run
% should be merged as long as their v peak ratio is >= mergethresh, but
% starting with the first pair whose ratio is < mergethresh, all the rest
% of the saccades in the run should be considered corrective and be trimmed
% away.  
% There may also be pairs that need to be "reverse trimmed", i.e. the first
% saccade in the pair should be deleted when the ratio is >= revtrimthresh.
% For the purpose of doing the subsequent curve fitting, it is sufficient
% to have a list of isolated single saccades to fit, and a separate list of
% pairs of starting and ending saccades of a run that is to be merged.
% Each isolated single saccade may have always been isolated, or it may be
% the first of a short-interval pair that failed the mergethresh criterion.
saccsdif = diff(bigsaccs);
is_short = saccsdif < Nminsaccint;
runstart = find([ is_short(1)
    is_short(2:end) & ~is_short(1:end-1) ]);
runend = find([ ~is_short(2:end) & is_short(1:end-1)
    is_short(end) ]);
shortintidx = find(is_short);
% <shortintidx>, <runstart>, <runend>, <isolatedsaccidx> are indices into
% bigsaccs; members of <isolatedsaccidx> are neither the first of a pair
% with a short interval, nor the second (hence the setdiff).
isolatedsaccidx = setdiff(find(~is_short), shortintidx + 1); 
mergeidx = [];  % col 1: starting saccidx; col 2: ending sacc idx
whittledDownToOneidx = [];
for runidx = 1:length(runstart)
    ratios = lfp_Samples{velo}( ...
        bigsaccs((runstart(runidx):runend(runidx))+1) ) ...
        / lfp_Samples{velo}(bigsaccs(runstart(runidx)));
    mergestart = 1;
    % The saccades to merge start with the first pair after those that need
    % to be reverse trimmed, and end with the last pair before those that
    % need to be trimmed.  However, since ratios are relative to the first
    % saccade, each saccade that gets reverse trimmed makes it necessary to
    % recalculate <ratios>:
    while mergestart < (runend(runidx) - runstart(runidx) + 1) ...
            && ratios(1) >= revtrimthresh
        mergestart = mergestart + 1;
        ratios = lfp_Samples{velo}( bigsaccs( ...
            ((runstart(runidx)+mergestart-1):runend(runidx))+1) ) ...
            / lfp_Samples{velo}(bigsaccs(runstart(runidx)+mergestart-1));
    end
    % <mergestart> now points at first sacc to merge.
    pairs2trimidx = find(ratios < mergethresh) + mergestart - 1;
    if isempty(pairs2trimidx)
        mergeend = length(ratios) + 1;
    else
        mergeend = pairs2trimidx(1);
    end
    % <mergeend> now points at the last sacc to merge
    if mergestart > mergeend
        warning('lfp_parseEyeTraces5:bogusrun', ...
            'BUG! Run from t = %.6f to t = %.6f is impossible.', ...
            lfp_index2time(bigsaccs(runstart(runidx))), ...
            lfp_index2time(bigsaccs(runend(runidx))) );
    elseif mergestart == mergeend
        % We're down to a single saccade at <mergestart>
        whittledDownToOneidx(end+1,1) = mergestart + runstart(runidx) - 1;
    else
        mergeidx(end+1, :) = [mergestart mergeend] + runstart(runidx) - 1;
    end
end

%
% Do the curve fitting to find saccstarts and saccends
%

saccstarts = zeros(size(bigsaccs));
saccends = zeros(size(bigsaccs));
singlesaccidx = [isolatedsaccidx; whittledDownToOneidx];

% First the isolated saccades:
h = waitbar(0,'Finding saccade starts and ends...');
totalwaitlength = (length(singlesaccidx) + size(mergeidx,1))*1.01;
smarkingwidth = 1.4; % used to = 2, then 1.8 changed by TMD 4/17/08
emarkingwidth = 1.7;
% <rmses> cells contain the following.  1: isolated; 2-4: multipeaks 2-4;
% 5: multipeaks > 4; 6: whittled Down To One   
rmses = cell(1,6); 
for saccidx2 = 1:length(singlesaccidx)
    saccidx = singlesaccidx(saccidx2);
    [c, numpeaks, gof] = fitsaccade( ...
        velo, bigsaccs, saccidx, threshes);
    if isempty(c)
        continue
    end
    if saccidx2 > length(isolatedsaccidx)
        rmses{6}(end+1) = gof.rmse;
    else
        rmses{1}(end+1) = gof.rmse;
    end
    saccstarts(saccidx) = round(c(2) - smarkingwidth*c(3)) ...
        + bigsaccs(saccidx) - 1;
    saccends(saccidx) = round(c(2) + emarkingwidth*c(3)) ...
        + bigsaccs(saccidx) - 1;
    fitrange = round(-3*c(3)):round(3*c(3));
    saccwaverange = fitrange + bigsaccs(saccidx) - 1;
    if saccwaverange(1) < 1
        fitrange = fitrange(saccwaverange>0);
        saccwaverange = saccwaverange(saccwaverange>0);
    end
    saccwave(saccwaverange) = saccwave(saccwaverange) + ...
        c(1)*exp(-((fitrange - c(2)) / c(3)).^2);
    waitbar(saccidx2/totalwaitlength, h);
end
% Then the merged runs.  smarkingwidth_merged/emarkingwidth_merged work
% the opposite way from smarkingwidth/emarkingwidth: they are multipliers
% for threshes.sacc to nudge that threshold up or down for marking start
% and end of saccade respectively, so increasing the values results in
% temporally tighter saccade markers.
smarkingwidth_merged = 1;
emarkingwidth_merged = 1;
for mergeidx2 = 1:size(mergeidx,1)
    saccidx = mergeidx(mergeidx2,:);
    c = zeros(12,1);
    c(3:3:length(c)) = 1; % these are divisors, so should not default to 0
    [coeffs, numpeaks, gof] = fitsaccade( ...
        velo, bigsaccs, saccidx, threshes);
    if numpeaks > 0 && numpeaks <=4
        rmses{numpeaks}(end+1) = gof.rmse;
    else
        rmses{5}(end+1) = gof.rmse;
    end
    if isempty(coeffs)
        warning('lfp_parseEyeTraces5:nofit', ...
            'Fitting for saccade at t=%s was skipped.', ...
            mat2str(lfp_index2time(bigsaccs(saccidx))) );
        continue
    end
    c(1:length(coeffs)) = coeffs;
    % For multiple peaks, there is no simple analytic formula for saccstart
    % and saccend in this case, so they must be found from the fitted
    % waveform.
    fitstartpt = min([
        round(c(2) - 3*c(3))
        round(c(5) - 3*c(6))
        round(c(8) - 3*c(9))
        round(c(10) - 3*c(11))
        ]);
    fitendpt = max([
        round(c(2) + 3*c(3))
        round(c(5) + 3*c(6))
        round(c(8) + 3*c(9))
        round(c(10) + 3*c(11))
        ]);
    fitrange = fitstartpt:fitendpt;
    saccwaverange = fitrange + bigsaccs(saccidx(1));
    % Adjust fitrange and saccwaverange if necessary:
    if saccwaverange(1) < 1
        fitrange = fitrange(saccwaverange>0);
        saccwaverange = saccwaverange(saccwaverange>0);
    end
    if saccwaverange(end) > numel(saccwave)
        fitrange = fitrange(saccwaverange<=numel(saccwave));
        saccwaverange = saccwaverange(saccwaverange<=numel(saccwave));
    end
    fittedwave = c(1)*exp(-((fitrange - c(2)) / c(3)).^2) + ...
        c(4)*exp(-((fitrange - c(5)) / c(6)).^2) + ...
        c(7)*exp(-((fitrange - c(8)) / c(9)).^2) + ...
        c(10)*exp(-((fitrange - c(11)) / c(12)).^2);
    saccwave(saccwaverange) = saccwave(saccwaverange) + fittedwave;
    startpt = 1;
    while(fittedwave(startpt) < threshes.sacc*smarkingwidth_merged ...
            && startpt < length(fittedwave))
        startpt =startpt + 1;
    end
    saccstarts(saccidx(1)) = startpt + saccwaverange(1) - 1;
    saccendpt = length(fittedwave);
    while(fittedwave(saccendpt) < threshes.sacc*emarkingwidth_merged ...
            && saccendpt > 1)
        saccendpt = saccendpt - 1;
    end
    saccends(saccidx(1)) = saccendpt + saccwaverange(1) - 1;
    waitbar((mergeidx2+length(singlesaccidx)) ...
        / totalwaitlength, h );
end
saccstarts(saccstarts==0) = [];
waitbar((mergeidx2 + length(singlesaccidx) + .005*totalwaitlength) ...
        / totalwaitlength );
saccends(saccends==0) = [];
close(h)

if length(blinkstarts) ~= length(blinkends)
    error('lfp_parseEyeTraces5:oops1', ...
        'Unbalanced blink events');
end
if length(saccstarts) ~= length(saccends)
    error('lfp_parseEyeTraces5:oops2', ...
        'Unbalanced saccade events');
end

% Construct fixations by maintaining three state variables: inblinkTO,
% inrecTO, insacc, while running through a time-sorted list of just those
% events.  A fixation starts when ALL variables become false, and ends when
% ANY variable becomes true.  The event codes in this strictly local event
% list are > 0 for event starts, < 0 for event ends, and magnitude 1, 2, 3
% for inblinkTO, inrecTO, insacc respectively.
events = [
    blinkstarts - NBlinkTO + 1, repmat(1, size(blinkstarts))
    blinkends - 1 + NBlinkTO - 1, repmat(-1, size(blinkends))
    lfp_RecSegments(:,2) - NRecTO + 1, repmat(2, size(lfp_RecSegments(:,1)))
    lfp_RecSegments(:,1) + NRecTO - 1, repmat(-2, size(lfp_RecSegments(:,1)))
    saccstarts, repmat(3, size(saccstarts))
    saccends - 1, repmat(-3, size(saccends))
    ];
events = sortrows(events, 1);
inblinkTO = false;
inrecTO = true;
insacc = false;
infix = false;
fixstarts = zeros(size(events,1),1);
fixends = zeros(size(events,1),1);
fixdiams = NaN(size(events,1),1);
fixnum = 1;
for k = 1:size(events,1)
    switch abs(events(k,2))
        case 1
            inblinkTO = events(k,2) > 0;
        case 2
            inrecTO = events(k,2) > 0;
        case 3
            insacc = events(k,2) > 0;
    end
    if infix
        if any([inblinkTO inrecTO insacc])
            infix = false;
            fixends(fixnum) = events(k,1) - 1;
            if fixends(fixnum) < fixstarts(fixnum) + NMinFixDur
                % Too short, back it out; note that it IS logically
                % possible for fixstarts(fixnum) to be greater than
                % fixends(fixnum):  
                fixstarts(fixnum) = 0;
                fixends(fixnum) = 0;
            else
                % save the "diameter" of the fixation:
                fixdiams(k) = max(abs( ...
                    max(lfp_Samples{eyeX}(fixstarts(fixnum):fixends(fixnum))) ...
                    - min(lfp_Samples{eyeX}(fixstarts(fixnum):fixends(fixnum))) ...
                    ), ...
                    abs( ...
                    max(lfp_Samples{eyeY}(fixstarts(fixnum):fixends(fixnum))) ...
                    - min(lfp_Samples{eyeY}(fixstarts(fixnum):fixends(fixnum))) ...
                    ));
%                 fixdiams(k) = sqrt( ...
%                     (lfp_Samples{eyeX}(fixends(fixnum)) ...
%                     - lfp_Samples{eyeX}(fixstarts(fixnum)))^2 ...
%                     + (lfp_Samples{eyeY}(fixends(fixnum)) ...
%                     - lfp_Samples{eyeY}(fixstarts(fixnum)))^2 );
                
                fixnum = fixnum + 1;
            end
        end
    else
        if ~any([inblinkTO inrecTO insacc])
            infix = true;
            fixstarts(fixnum) = events(k,1) + 1;
        end
    end
end
fixstarts(fixstarts==0) = [];
fixends(fixends==0) = [];
fixdiams(isnan(fixdiams)) = [];

% Find "normal" size for fixations if possible; if not default to 100
% pixels.
maxdiam = 100;
if length(fixdiams) > 50
    binedges = floor(min(fixdiams)):ceil(max(fixdiams)*(1+eps));
    binctrs = binedges(1:end-1) + .5;
    rawhisto = histc(fixdiams, binedges);
    smoothpts = 5;
    smoothfunc = hanning(2*smoothpts+1);
    smoothed = conv(rawhisto, smoothfunc) / sum(smoothfunc);
    smoothedhisto = smoothed(1+smoothpts : end-smoothpts);
    [C,I] = max(smoothedhisto);
    maxbin = I(1);
    diffmaxch = 1;
    diffminch = 1e-4;
    tolfun = 1e-2;
    tolx = 1e-3;
    maxfunevals = 600;
    fitlimit = min(10*binctrs(maxbin), length(binctrs));
    x = reshape(binctrs(binctrs<fitlimit), [], 1);
    yscale = max(smoothedhisto(binctrs<fitlimit));
    y = reshape( ...
        smoothedhisto(binctrs<fitlimit) / yscale, [], 1 );
    % This is a theoretically groundless formula, but the problem is
    % actually rather complicated because (to start with) simple isotropic
    % IID diffusion is NOT a good model for microsaccades, which at the
    % very least tend to revert to the initial position.
    [cfun,gof,output] = fit(x, y, 'gauss2', ...
        fitoptions( 'Method','NonlinearLeastSquares', ...
        'DiffMaxChange', diffmaxch, 'DiffMinChange', diffminch, ...
        'TolFun', tolfun, 'TolX', tolx', ...
        'MaxFunEvals', maxfunevals, ...
        'StartPoint', [1 50 20 1 100 20], ...
        'Lower', [0 0 0 0 0 0] ));
    c = coeffvalues(cfun);
    c([1 4]) = yscale * c([1 4]);
    % single peak fit:
    fitted = c(1)*exp(-((binctrs' - c(2)) / c(3)).^2);
    ratio = smoothedhisto(1:end-1) ./ fitted;
    middiambin = find(ratio(maxbin:end) > 2) + maxbin - 1;
    middiambin = middiambin(1);
    % double peak fit:
    fitted2 = c(1)*exp(-((binctrs' - c(2)) / c(3)).^2) + ...
        c(4)*exp(-((binctrs' - c(5)) / c(6)).^2);
    ratio = smoothedhisto(1:end-1) ./ fitted2;
    maxdiambin = find(ratio(maxbin:end) > 2) + maxbin - 1;
    maxdiambin = maxdiambin(1);
    hF = figure;
    bins2plot = 1:min(fitlimit*10,length(binctrs));
    maxdiam = binctrs(maxdiambin);
    middiam = binctrs(middiambin);
    plot(binctrs(bins2plot), [ smoothedhisto(bins2plot) ...
        fitted(bins2plot) fitted2(bins2plot) ]);
    hold on;
    plot([maxdiam maxdiam], get(gca, 'YLim'), 'k');
    plot([middiam middiam], get(gca, 'YLim'), 'm');
    legend({'fixdiams' 'fitted' 'fitted2' 'maxdiam' 'middiam'});
    figname = sprintf('%s fixation diameters', fragname);
    title(figname, 'Interpreter', 'none');
    saveas(hF, fullfile(lfp_DataDir, [figname '.fig']));
    close(hF);
end

% Find the fixations that are abnormally large.  Mark them right in the
% middle:
bigfixidx = find(fixdiams>maxdiam);
bigfix = round( ...
    (fixstarts(bigfixidx) + fixends(bigfixidx)) / 2 );
midfixidx = find(fixdiams>middiam);
midfix = round( ...
    (fixstarts(midfixidx) + fixends(midfixidx)) / 2 ) - 5;

% Adjust saccade ends that precede big fixations until they're not big
% anymore.  Try every valley and threshes.sacc negative crossing that
% occurs between the preceding sacc end and minfixdur before the big fix
% end.
fNTv = lfp_findNegThresh3(velo, threshes.sacc, 'sample');
bigfixthresh = middiam;
for saccfix2fix = midfixidx'
    prevsaccnum = find(saccends <= (fixstarts(saccfix2fix)));
    if isempty(prevsaccnum)
        % There is nothing to adjust.  
        warning('lfp_parseEyeTraces5:oops', 'What are we doing here?');
        continue
    end
    prevsaccnum = prevsaccnum(end);
    samplerange = (saccends(prevsaccnum) + 1) ...
        : (fixends(saccfix2fix) - NMinFixDur);
    if isempty(samplerange)
        % There is not chance to adjust.
        warning('lfp_parseEyeTraces5:badfix2', ...
            'The saccade-fixation pair at %.6f - %.6f s cannot be adjusted.', ...
            lfp_index2time(saccstarts(prevsaccnum)), ...
            lfp_index2time(fixends(saccfix2fix)) );
        continue
    end
    candidates = fNTv(fNTv > samplerange(1) & fNTv < samplerange(end))';
    adjusted = false;
    for candidx = 1:length(candidates)
        candiam = max(abs( ...
            max(lfp_Samples{eyeX}(candidates(candidx):fixends(saccfix2fix))) ...
            - min(lfp_Samples{eyeX}(candidates(candidx):fixends(saccfix2fix))) ...
            ), ...
            abs( ...
            max(lfp_Samples{eyeY}(candidates(candidx):fixends(saccfix2fix))) ...
            - min(lfp_Samples{eyeY}(candidates(candidx):fixends(saccfix2fix))) ...
            ));
        if candiam <= bigfixthresh
            % adjustment found; save new values
            saccends(prevsaccnum) = candidates(candidx);
            fixstarts(saccfix2fix) = candidates(candidx) + 1;
            adjusted = true;
            break
        end
    end
    if ~adjusted
        warning('lfp_parseEyeTraces5:badfix', ...
            'The saccade-fixation pair at %.6f - %.6f s cannot be adjusted.', ...
            lfp_index2time(saccstarts(prevsaccnum)), ...
            lfp_index2time(fixends(saccfix2fix)) );
    end
end

% Convert from sample indices to timestamps and construct result:
ts = {lfp_index2time(blinkstarts) lfp_index2time(blinkends-1) ...
    lfp_index2time(saccstarts) lfp_index2time(saccends-1) ...
    lfp_index2time(bigsaccs(mergeidx(:,1))) ...
    lfp_index2time(bigsaccs(mergeidx(:,2))) ...
    lfp_index2time(fixstarts) lfp_index2time(fixends) ...
    lfp_index2time(bigfix) lfp_index2time(midfix) ...
    };
end


function [c, numpeaks, gof] = fitsaccade( ...
    velo, bigsaccs, saccidx, threshes)
% If saccidx is a scalar, then one velo peak is fitted.
% If saccidx is a two-element vector, then the entire interval extending
% from the first peak to the second (with any intervening peaks) is fitted
% unconditionally.
% In either case, the interval to be fitted is found by starting at the
% peaks given in saccidx, and moving outwards until either the velocity
% falls below the fixed threshold threshes.sacc, or the velocity begins to
% increase again (i.e. a strict minimum is found).
global lfp_Samples lfp_SamplePeriod
c = [];
diffmaxch = 1;
diffminch = 1e-4;
tolfun = 1e-2;
tolx = 1e-3;
maxfunevals = 600;
vmax = max(lfp_Samples{velo}(bigsaccs(saccidx)));
startpt = bigsaccs(saccidx(1));
endpt = bigsaccs(saccidx(end));
% If there's more than one peak internal to the range being merged, that by
% definition makes it a LNM (Loch Ness Monster), which means it has to be
% fit all the way to the tip of its boat-smashing, terrorizing tail:
y = lfp_Samples{velo}(startpt:endpt);
numinternalpeaks = sum( ...
    (y(2:end-1) > y(1:end-2) & y(2:end-1) > y(3:end)) );
LNM = numinternalpeaks > 1;
while startpt > 0 && lfp_Samples{velo}(startpt) > threshes.sacc ...
        && (lfp_Samples{velo}(startpt) >= lfp_Samples{velo}(startpt-1))
    startpt = startpt - 1;
end
while endpt < numel(lfp_Samples{velo}) && ...
        lfp_Samples{velo}(endpt) > threshes.sacc ...
        && (LNM || ...
        (lfp_Samples{velo}(endpt) >= lfp_Samples{velo}(endpt+1)) )
    endpt = endpt + 1;
end
% bigsaccs(saccidx(1)) is arbitrarily used as a reference because it's
% conveniently available:
x = (startpt:endpt) - bigsaccs(saccidx(1));
% must normalize because 'fit' can't cope with radically different scales
% along different dimensions of the param search space:
y = lfp_Samples{velo}(startpt:endpt)/vmax;
if length(saccidx) == 1
    [peakvals, peakposns] = max(y);
    numpeaks = 1;
else
    ispeak = [false ...
        (y(2:end-1) > y(1:end-2) & y(2:end-1) > y(3:end)) ...
        false];
    peakvals = y(ispeak);
    peakposns = x(ispeak);
    numpeaks = sum(ispeak);
end
switch numpeaks
    case 1
        [cfun,gof,output] = fit(x', y', 'gauss1', ...
            fitoptions( 'Method','NonlinearLeastSquares', ...
            'DiffMaxChange', diffmaxch, 'DiffMinChange', diffminch, ...
            'TolFun', tolfun, 'TolX', tolx', ...
            'MaxFunEvals', maxfunevals, ...
            'Lower', repmat([0 x(1) 1], 1, length(peakvals)), ...
            'Upper', repmat([Inf Inf 0.02/lfp_SamplePeriod], ...
            1, length(peakvals)) ));
    case 2
        [cfun,gof,output] = fit(x', y', 'gauss2', ...
            fitoptions( 'Method','NonlinearLeastSquares', ...
            'DiffMaxChange', diffmaxch, 'DiffMinChange', diffminch, ...
            'TolFun', tolfun, 'TolX', tolx', ...
            'MaxFunEvals', maxfunevals, ...
            'StartPoint', reshape([ peakvals
            peakposns
            ones(size(peakvals)) ], 1, []), ...
            'Lower', repmat([0 x(1) 1], 1, length(peakvals)), ...
            'Upper', repmat([Inf Inf 0.02/lfp_SamplePeriod], ...
            1, length(peakvals)) ));
    case 3
        [cfun,gof,output] = fit(x', y', 'gauss3', ...
            fitoptions( 'Method','NonlinearLeastSquares', ...
            'DiffMaxChange', diffmaxch, 'DiffMinChange', diffminch, ...
            'TolFun', tolfun, 'TolX', tolx', ...
            'MaxFunEvals', maxfunevals, ...
            'StartPoint', reshape([ peakvals
            peakposns
            ones(size(peakvals)) ], 1, []), ...
            'Lower', repmat([0 x(1) 1], 1, length(peakvals)), ...
            'Upper', repmat([Inf Inf 0.02/lfp_SamplePeriod], ...
            1, length(peakvals)) ));
    otherwise
        % 3 might actually be enough, but surely 4 can fit anything
        [peakvals, ix] = sort(peakvals, 'descend');
        peakvals = peakvals(1:4);
        peakposns = peakposns(ix(1:4));
        [cfun,gof,output] = fit(x', y', 'gauss4', ...
            fitoptions( 'Method','NonlinearLeastSquares', ...
            'DiffMaxChange', diffmaxch, 'DiffMinChange', diffminch, ...
            'TolFun', tolfun, 'TolX', tolx', ...
            'MaxFunEvals', maxfunevals, ...
            'StartPoint', reshape([ peakvals
            peakposns
            ones(size(peakvals)) ], 1, []), ...
            'Lower', repmat([0 x(1) 1], 1, length(peakvals)), ...
            'Upper', repmat([Inf Inf 0.02/lfp_SamplePeriod], ...
            1, length(peakvals)) ));
end
% Warn on fits where rmse exceeds the 99.8th percentile of
% Y081007\Sc9correct or where the fitting algorithm gave up:
tails998 = [0.16 0.11 0.10 0.066 0.16];
if gof.rmse > tails998(length(peakvals))
    warning('lfp_parseEyeTraces5:dubious', ...
        'High RMS error %.3f for %d-peak fit at %.6f to %.6f s', ...
        gof.rmse, length(peakvals), lfp_index2time(startpt), lfp_index2time(endpt) );
end
if output.funcCount >= maxfunevals
    warning('lfp_parseEyeTraces5:dubious2', ...
        'Gave up after %d func evals in fit at %.6f to %.6f s', ...
        output.funcCount, ...
        lfp_index2time(startpt), lfp_index2time(endpt) );
end
c = coeffvalues(cfun);
% De-normalize the vertical scale params:
for k = 1:length(c)/3
    c((k-1)*3 + 1) = c((k-1)*3 + 1) * vmax;
end

end
