function [blinkstarts, blinkends] = ...
    lfp_findBlinks(eyeX, eyeY, velo, threshes)
%INPUTS
% eyeX, eyeY, velo: filenums of X, Y, velocity
%   traces ("velocity" here really means "speed", i.e. the magnitude of the
%   velocity vector)
% threshes: a structure containing the thresholds used. The fields in the
%   structure are:
%   blinkTO - time out around each blink during which lfp_parseEyeTraces7
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
% <blinkstarts>, <blinkends> are sample indices.
%OUTPUTS
% <blinkstarts>, <blinkends>: sample indices.
%NOTES
%   A prototypical 'blink' is defined in terms of a chronological
% sequence of conditions being met as task time proceeds.  These conditions
% are represented in the code as <condition_1> through <condition_8> (more
% details below).
%   Depending on recording and calibration conventions, eye position
% data can potentially be inverted, so positive blink thresholds are taken
% to mean that blinks are represented as maximal positive values, and
% negative blink thresholds are taken to mean that blinks go negative. Note
% that "velocity" values cannot be negative, whereas position values can
% be, and so can position thresholds. The words "subthreshold" and
% "suprathreshold" are used below to indicate whether or not the event the
% threshold is supposed to detect has occurred.  E.g. blink thresholds that
% are negative detect blinks when the position value is more negative, and
% blink thresholds that are positive detect blinks when the position is
% more positive.
%   Real blinks often fail to be prototypical, so voting logic is applied
% to score the candidate blink events. The meanings of <condition_1>
% through <condition_8> are as follows:
%   1. a period of 2*blinkTO sec whose median X value is subthreshold for
%   blinkx and median Y value is subthreshold for blinky
%   2. (There is no longer any condition 2.)
%   3. a v peak that exceeds blinkv2, or a threshold crossing away from
%   zero of blinkx or blinky
%   4. a series of samples whose duration does not exceed blinkdur
%   5. median X of that series of samples is suprathreshold to blinkx or
%   median Y is suprathreshold to blinky
%   6. another v peak that exceeds blinkv2
%   7. (There is no longer any condition 7.)
%   8. another period of 2*blinkTO sec as in condition 1
% The algorithm for finding blinks is to find candidate blink events,
% defined as either a v peak that exceeds blinkv2, a blinkx crossing, or a
% blink y crossing (condition 3 above). Taking the first candidate blink
% event as a putative blink start, we test the conditions as follows,
% accumulating votes in categories a, b, and c (the other letters below
% denote aditional comments, not voting categories):
%   a. Each v peak > blinkv1 within +/- 2*blinksrc of the event casts one
%   <vote_a>.  Note that extreme glitchiness may produce more than one
%   peak, making it possible to overcome vote deficiencies when examining
%   the putative blink end (see below).
%   b. Each v peak > blinkv2 within +/- 2*blinksrc of the event casts one
%   <vote_b>.  Together with a., EACH of the peaks that meet this condition
%   ends up cumulatively casting two votes.
%   u. Condition 1 is strictly required.
%   v. The next candidate event satisfying conditions 3 and 4 is now
%   considered as a putative blink end event, and condition 5 must be met
%   for at least one of X and/or Y.  If the required conditions are not met
%   at this point, the tally of conditions met and failed up to this point
%   is logged, and the process is repeated with the next event satisfying
%   conditions 3 and 4.  If there are no such events, then the putative
%   Start Blink event is marked as a Big Saccade, and the timestamp of the
%   Big Saccade is logged.
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
% marked as a Big Saccade, and its timestamp is logged.  In either case,
% the next candidate blink event that was not part of the blink or Big
% Saccade is used as the next putative Start Blink.
%
% MAINTENANCE NOTES:  
% 3) VARIABLE NAMES:  variables that contain sample numbers or have nothing
%   whatsoever to do with sample numbers have no special suffixes.
%   However, indices into lists of sample numbers will have the suffix idx;
%   indices into lists of indices into sample numbers will be idx2; etc.
%   Timestamps in seconds (rather than samples) will have the suffix TS.
%   Iterators that index something that has nothing to do with sample
%   numbers have the suffix num.

%$Rev: 422 $
%$Date: 2023-09-08 14:25:38 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

global lfp_Samples lfp_SamplePeriod lfp_RecSegments %#ok<GVMIS> 
lfp_log('Starting lfp_findBlinks');

Nsamples = numel(lfp_Samples{velo});

% Derived constants:
NBlinkTO = round(threshes.blinkTO / lfp_SamplePeriod);
NBlinkDur = round(threshes.blinkdur / lfp_SamplePeriod);
NBlinkSrc = round(threshes.blinksrc / lfp_SamplePeriod);

% Find sample indicies of candidate velocity peaks
% This is the master list for both blinks and saccades.
vpeaks = lfp_findPeaks(velo, threshes.blinkv1);
if threshes.blinkx > 0
    fPTx = lfp_findPosThresh3(eyeX, threshes.blinkx, 'sample')';
    fNTx = lfp_findNegThresh3(eyeX, threshes.blinkx, 'sample')';
else
    fPTx = lfp_findNegThresh3(eyeX, threshes.blinkx, 'sample')';
    fNTx = lfp_findPosThresh3(eyeX, threshes.blinkx, 'sample')';
end
if threshes.blinky > 0
    fPTy = lfp_findPosThresh3(eyeY, threshes.blinky, 'sample')';
    fNTy = lfp_findNegThresh3(eyeY, threshes.blinky, 'sample')';
else
    fPTy = lfp_findNegThresh3(eyeY, threshes.blinky, 'sample')';
    fNTy = lfp_findPosThresh3(eyeY, threshes.blinky, 'sample')';
end
% find sample indices of candidate blink events <canblinks>:
canblinks = ...
    unique( [
    vpeaks(lfp_Samples{velo}(vpeaks) > threshes.blinkv2)
    fPTx
    fNTx
    fPTy
    fNTy
    lfp_RecSegments(:,2)
    ] );

% canblinks and vpeaks are indices into lfp_Samples.
% canblinksidx, etc., are indices into canblinks and vpeaks; they are used
% to obviate the need to search all of canblinks and vpeaks on each
% iteration.
blinkstarts = [];
blinkends = [];
canblinksidx = 1;  % index into canblinks
vpeaksidx = 1;  % starting index into vpeaks for votes a & b search
fprintf('Vetting candidate blinks...\n');
prevcanblinksidx = 0;

% The whole loop that follows has gradually grown to Gargantuan propertions
% and should really re-factored.  The goal here is to find the longest
% blink start-end pair where the start is the current candidate in
% <canblinks> and the end is subject to some hideous combination of
% conditions 1-8 and votes a-c.  See note "w." and following for more
% details.
while canblinksidx < length(canblinks)
    % <bsamp> is sample # of current putative blink start:
    bsamp = canblinks(canblinksidx);

    % Find range of <vpeaks> <vpeaksidx_a:vpeaksidx_b> to search for votes
    % a and b; <vpeaksidx_a> and <vpeaksidx_b> do NOT correspond to
    % <vote_a> and <vote_b>!  <vpeaksidx> gets updated by the end of each
    % iteration to contain the proper value to use for <vpeaksidx_a> on the
    % following iteration.
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

    % Compute <vote_a> and <vote_b>
    vote_a = sum(lfp_Samples{velo}(vpeaks(vpeaksidx_a:vpeaksidx_b)) ...
        > threshes.blinkv1 );
    vote_b = sum(lfp_Samples{velo}(vpeaks(vpeaksidx_a:vpeaksidx_b)) ...
        > threshes.blinkv2 );
    vote_c = NaN;

    % Compute <condition_1> (does not depend on <vpeaksidx_a>,
    % <vpeaksidx_b>, <vote_a>, or <vote_b>):
    if threshes.blinkx > 0
        cond1x = median( ...
            lfp_Samples{eyeX}( max(1, bsamp-2*NBlinkTO) : bsamp )) ...
            < threshes.blinkx;
    else
        cond1x = median( ...
            lfp_Samples{eyeX}( max(1, bsamp-2*NBlinkTO) : bsamp )) ...
            > threshes.blinkx;
    end
    if threshes.blinky > 0
        cond1y = median( ...
            lfp_Samples{eyeY}( max(1, bsamp-2*NBlinkTO) : bsamp )) ...
            < threshes.blinky;
    else
        cond1y = median( ...
            lfp_Samples{eyeY}( max(1, bsamp-2*NBlinkTO) : bsamp )) ...
            > threshes.blinky;
    end
    condition_1 = cond1x && cond1y; % time before vpeak is < blinkx & blinky
    condition_8 = NaN;

    % Search <canblinks> for a matching blink end starting with the next
    % candidate and running until we exceed <NBlinkDur> (which may happen
    % immediately).  <goodbend == 0> implies that none was found.  Start
    % searching at the next candidate blink after <canblinksidx>.
    goodbend = 0;
    logmsg = sprintf('***** Re: canblinksidx = %d *****\n', canblinksidx);
    bendidx = canblinksidx + 1;
    while bendidx <= length(canblinks) && ...
            canblinks(bendidx) <= bsamp + NBlinkDur
        % Compute <condition_5>:
        if threshes.blinkx > 0
            cond5x = median( ...
                lfp_Samples{eyeX}(bsamp+1 : canblinks(bendidx)-1)) ...
                >= threshes.blinkx;
        else
            cond5x = median( ...
                lfp_Samples{eyeX}(bsamp+1 : canblinks(bendidx)-1)) ...
                <= threshes.blinkx;
        end
        if threshes.blinky > 0
            cond5y = median( ...
                lfp_Samples{eyeY}(bsamp+1 : canblinks(bendidx)-1)) ...
                >= threshes.blinky;
        else
            cond5y = median( ...
                lfp_Samples{eyeY}(bsamp+1 : canblinks(bendidx)-1)) ...
                <= threshes.blinky;
        end
        condition_5 = cond5x || cond5y;

        if condition_1 && condition_5
            % Conditions met so far; keep evaluating. This case always
            % either increments bendidx or breaks out of the 'while
            % bendidx' loop.

            % Compute <vote_c>, for which recording breaks are a special
            % case considered equivalent to "> threshes.blinkv2".
            isrecbreak = ismember(bendidx, lfp_RecSegments(:,2));
            if isrecbreak
                vote_c = 2;
            else
                vote_c = (lfp_Samples{velo}(canblinks(bendidx)) ...
                    > threshes.blinkv1 ) ...
                    + (lfp_Samples{velo}(canblinks(bendidx)) ...
                    > threshes.blinkv2 );
            end

            % Compute <condition_8>:
            if threshes.blinkx > 0
                cond8x = median(lfp_Samples{eyeX}( ...
                    canblinks(bendidx) ...
                    : min( Nsamples, canblinks(bendidx) + 2*NBlinkTO ) )) ...
                    < threshes.blinkx;
            else
                cond8x = median(lfp_Samples{eyeX}( ...
                    canblinks(bendidx) ...
                    : min( Nsamples, canblinks(bendidx) + 2*NBlinkTO ) )) ...
                    > threshes.blinkx;
            end
            if threshes.blinky > 0
                cond8y = median(lfp_Samples{eyeY}( ...
                    canblinks(bendidx) ...
                    : min(Nsamples, canblinks(bendidx) + 2*NBlinkTO) )) ...
                    < threshes.blinky;
            else
                cond8y = median(lfp_Samples{eyeY}( ...
                    canblinks(bendidx) ...
                    : min(Nsamples, canblinks(bendidx) + 2*NBlinkTO) )) ...
                    > threshes.blinky;
            end
            condition_8 = cond8x && cond8y;

            % Combine <condition_8> with the votes and <goodbend> to decide
            % whether to keep looking for a later <goodbend>, or declare
            % victory and break from 'while bendix' loop.
            if goodbend
                % This case will either increment bendidx or break the
                % 'while bendidx' loop:
                if (condition_8 || isrecbreak)...
                        && (vote_a + vote_b + vote_c > 2)
                    % another good blink end, keep trying to extend
                    goodbend = canblinks(bendidx);
                    logmsg = [ logmsg sprintf(...
                        'Extended blink to %.6f; a:%d b:%d c:%d\n', ...
                        lfp_index2time(goodbend), ...
                        vote_a, vote_b, vote_c )]; %#ok<*AGROW> 
                    bendidx = bendidx + 1;
                else
                    % existing <goodbend> is the final Blink End
                    break; % from 'while bendix' loop
                end
            else
                % This case always increments bendidx:
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
            % Initial conditions (i.e. 1 and 5) not met; log as Big
            % Saccade. This case always either increments bendidx ("keep
            % looking") or breaks out of the 'while bendidx' loop ("quit
            % looking":
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
                % <goodbend> is valid
                break; % from'while bendix' loop
            else
                % no <goodbend> yet (this is always the case on the first
                % iteration)
                bendidx = bendidx + 1;
            end
        end
    end
    % The Search For A Blink End is now over, for better or worse.
    % canblinksidx gets incremented by at least 1, one way or another:
    if goodbend
        % 'lfp_log' opens and closes the log file, so it should not be
        % called too often:
        if length(logmsg) > 1e5
            lfp_log(logmsg);
            logmsg = '';
        end
        blinkstarts(end+1,1) = bsamp;
        blinkends(end+1,1) = goodbend - 1;
        % One way or another, bendidx is now pointing at the first
        % candidate blink event that is not to be part of the blink, so
        % this becomes our next starting point for blink detection.  It
        % could also be pointing past the end of canblinks, however, in
        % which case the 'while bendidx' loop will exit after this
        % iteration.
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
        % 'lfp_log' opens and closes the log file, so it should not be
        % called too often:
        if length(logmsg) > 1e5
            lfp_log(logmsg);
            logmsg = '';
        end
        % Try again with the next candidate blink event:
        canblinksidx = canblinksidx + 1;
    end
    % Equivalent to a progress bar, but faster:
    if canblinksidx - prevcanblinksidx > 1000
        fprintf( '%d/%d candidate blinks processed.\n', ...
            canblinksidx, length(canblinks) );
        prevcanblinksidx = canblinksidx;
    end
end
fprintf('Done vetting candidate blinks.\n');
if ~isempty(logmsg)
    lfp_log(logmsg);
end

if length(blinkstarts) ~= length(blinkends)
    error('lfp_findBlinks:oops1', ...
        'Unbalanced blink events');
end
if isempty(blinkstarts)
    warning('lfp_findBlinks:oops2', ...
        'There are no blink events!');
    lfp_log('There are no blink events!');
end

