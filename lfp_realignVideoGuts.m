function [new_VT_timestamps, reportstr, scorecard] = ...
    lfp_realignVideoGuts(testevts, testcol, VT_timestamps, ...
    maxpixerr, ms2, XYchan, frametol)
% Computes an offset separately for each trial such that when video data is
% shifted in time, the video positions match the event positions given in
% the mazespec2 <ms2>.  Note that if <ms2> is bad, the entire process is
% hopeless.
%INPUTS
% testevts - a cell vector of lists of alternative events to use as
%   alignment references.
% testcol - vector of same length as testevts containing column numbers of
%   the relevant coordinate for each testevt, i.e. 1 for X and 2 for Y.
% VT_timestamps - in same units as lfp_Events.
% maxpixerr - threshold used to determine whether each event is
%   consistently aligned between video and events file.
% ms2 - a mazespec2 (see lfp_measureTMaze2).  May contain event IDs that
%   are not in <testevts>, but must contain every element of <testevts>, as
%   tested with isequal.
% XYchan - a two-element vector containing channel numbers for X trace and
%   Y trace, respectively.
% frametol - specifies the maximum number of frames of slippage that is
%   considered to be tolerable.  As we iterate through the trials, the time
%   slippage is assumed to zero until a trial is encountered that has a
%   time offset whose magnitude is more than <frametol> frames.  Once an
%   offset has been applied, the same offset will be used for all
%   subsequent trials until a trial is encountered whose offset differs
%   from the previous trial's by more than <frametol> frames.

%$Rev: 91 $
%$Date: 2009-11-02 18:05:07 -0500 (Mon, 02 Nov 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if size(VT_timestamps,1) ~= 1
    error('lfp_realignVideoGuts:VT_timestamps', ...
        'VT_timestamps must be a row vector');
end
reportstr = '';
dt = diff(VT_timestamps);
framedur = median(dt);  % first approximation
framedur = median(dt(dt<(1.1*framedur)));   % real deal
TStol = framedur * frametol;
distance = NaN(size(testevts));
new_VT_timestamps = VT_timestamps;
scorecard = cell(length(lfp_SelectedTrials), 1);
xings = cell(size(testevts));

for e = 1:length(testevts)
    ismazerow = cellfun(@isequal, ...
        ms2.evtIDs, repmat(testevts(e), size(ms2.evtIDs)) );
    xings{e} = lfp_findPosThresh3(XYchan(testcol(e)), ...
        ms2.medians(ismazerow, testcol(e)), 'sample');
end

prevoffset = 0;
isVTgap = [(dt > 1.5 * framedur) false];
mintrialdur = 2;    % minimum plausible trial duration in seconds
for t = 1:length(lfp_SelectedTrials)
    ts = NaN(size(testevts));   % working storage for evt timestamps
    % Use the realigned timestamps to identify the frames that belong to
    % the trial, since cumulative can add up to more than a trial duration
    trialframes = ...
        (new_VT_timestamps >= lfp_Events(lfp_TrialIndex(t,1),1)) ...
        & (new_VT_timestamps <= lfp_Events(lfp_TrialIndex(t,2),1));
    trialframeidx = find(trialframes);
    if length(trialframeidx) < mintrialdur/framedur
        % There could be such a huge increment in slippage during this
        % particular ITI that even the <new_VT_timestamps> don't work.
        % Look for the next sequence of timestamps that could possibly be a
        % trial, starting with the first timestamp gap after the end of the
        % previous trial.  "Possibly a trial" is defined as a series of
        % of frames without gaps of >= than <mintrialdur> duration.
        if t == 1
            slippedPrevTrialEnd = VT_timestamps(1);
        else
            slippedPrevTrialEnd = lfp_Events(lfp_TrialIndex(t-1, 2), 1) ...
                - prevoffset;
        end
        firstframeidx = find(VT_timestamps > slippedPrevTrialEnd);
        if isempty(firstframeidx)
            warning('lfp_realignVideoGuts:firstframeidx', ...
                'There are no more VT timestamps to match to trial %d', t);
        else
            if ~any(isVTgap)
                error('This is too hard, I give up!');
            else
                nxgapidx = find(VT_timestamps > slippedPrevTrialEnd ...
                    & isVTgap);
                firstframeidx = nxgapidx(1) + 1;
                while (length(nxgapidx) > 1) && ...
                        ( VT_timestamps(nxgapidx(2)) - ...
                        VT_timestamps(firstframeidx) < mintrialdur )
                    nxgapidx(1) = [];
                    firstframeidx = nxgapidx(1) + 1;
                end
            end
            % Guess that the first frame is 2 s before trial start:
            offsetguess = lfp_Events(lfp_TrialIndex(t,1),1) - 2 -...
                VT_timestamps(firstframeidx);
            recstartTS = lfp_Events(lfp_TrialIndex(t,1),1) - 2;
            trialframes = ...
                ( (VT_timestamps + offsetguess) >= recstartTS ) ...
                & ( (VT_timestamps + offsetguess) <= ...
                lfp_Events(lfp_TrialIndex(t,2),1) );
            trialframeidx = find(trialframes);
            if isempty(trialframeidx)
                warning('lfp_realignVideoGuts:badtrial2', ...
                    'Could not find corresponding video timestamps for trial %d', t);
                continue
            end
            offsetguess_samples = round(offsetguess/lfp_SamplePeriod);
        end
    else
        offsetguess_samples = round(prevoffset/lfp_SamplePeriod);
    end
    for e = 1:length(testevts)
        % find the event:
        evtidx = find(ismember( ...
            lfp_Events(lfp_TrialIndex(t,1):lfp_TrialIndex(t,2), 2), ...
            testevts{e} )) ...
            + lfp_TrialIndex(t,1) - 1;
        if length(evtidx) > 1
            warning('lfp_realignVideoGuts:multievt', ...
                'Trial %d contains more than one event %s; using first', ...
                t, testevts{e} );
            evtidx = evtidx(1);
        end
        % compute distance from coordinate specified in <ms2>:
        if ~isempty(evtidx)
            ts(e) = lfp_Events(evtidx, 1);
            evtID(e) = lfp_Events(evtidx, 2);
            sampidx = lfp_time2index(ts(e));
            ismazerow = cellfun(@ismember, ...
                repmat({evtID(e)}, size(ms2.evtIDs)), ms2.evtIDs );
            distance(e) = abs( lfp_Samples{XYchan(testcol(e))}(sampidx) ...
                - ms2.medians(ismazerow, testcol(e)) );
        end
    end
    if all(distance <= maxpixerr)
        scorecard{t} = 'good';
    else
        % To find offsets for all non-NaN ts, get timestamps of
        % matching video threshold crossing.  <offsets> are defined as
        % numbers that can be added to the video timestamps to make the
        % video threshold crossings agree with the recorded event times.
        offsets = cell(size(testevts));
        for e = 1:length(testevts)
            if isnan(ts(e))
                offsets{e} = NaN;
            else
                myxings = xings{e}( ...
                    xings{e} >= (lfp_TrialIndex(t,3) - offsetguess_samples) ...
                    & xings{e} <= (lfp_TrialIndex(t,4) - offsetguess_samples) );
                if isempty(myxings)
                    offsets{e} = NaN;
                else
                    offsets{e} = ts(e) - lfp_index2time(myxings);
                end
            end
        end
        % Find median offset for all events, throwing out the most extreme
        % value if there is an even number of values.  Blindly just use all
        % the values for events with multiple thresh xings. 
        metamedian = nanmedian(cell2mat(offsets));
        if isnan(metamedian) || isinf(metamedian)
            warning('lfp_realignVideoGuts:badoffset3', ...
                'Could not calculate median offset for trial %d', t);
            continue
        end
        % Choose the closest value to <metamedian> for events that had
        % multiple crossings:
        xinglens = cellfun(@length, offsets);
        for e = find(xinglens>1)
            [C,I] = min(abs(offsets{e} - metamedian));
            offsets{e} = offsets{e}(I);
        end
        % All cells in <offsets> are now scalars.  <offsetvals> is a vector
        % of offsets, one for each event in same order as <testevts>.
        offsetvals = cell2mat(offsets);
        offsets2use = true(size(offsets));
        % If there is any even number of offsets, throw out the most
        % extreme so that the median does not get interpolated.
        if mod(sum(~isnan(offsetvals)),2) == 0
            [C,I] = max(abs(offsetvals - metamedian));
            offsets2use(I) = false;
        end
        medoffset = nanmedian(offsetvals(offsets2use));
        % For any offset that differs by > 2 frames from the median,
        % determine if it could be a premature photobeam trigger, i.e. if
        % the recorded event time is excessively early, making the offset
        % less than it should be.
        isearly = (offsetvals - medoffset) < 0;
        isbad = ~isearly & (abs(offsetvals - medoffset) > 2 * framedur);
        reportstr = sprintf('%s\nAdjusting trial %d: ', reportstr, t);
        % Compute realigned new_VT_timestamps, using "trialframeidx(1):end"
        % indices to propagate each trial's offset through the following
        % ITI.
        if any(isbad)
            % determine whether there are gaps in the video frame
            % ts sequence:
            trialframeTS = VT_timestamps(trialframes);
            trialframedur = diff(trialframeTS);
            if any( trialframedur > 1.5 * framedur ...
                    | trialframedur < 0.5 * framedur )
                % Got gaps, i.e. intra-trial slippage.
                gapstring = '';
                for f = find( trialframedur > 1.5 * framedur ...
                        | trialframedur < 0.5 * framedur )
                    gapstring = sprintf( '%s\n%.6f - %.6f', ...
                        gapstring, trialframeTS(f), trialframeTS(f+1) );
                end
                warning('lfp_realignVideoGuts:gaps', ...
                    'Trial %d has intra-trial slippage between frames timestamped:%s', ...
                    t, gapstring );
                % If necessary, a mind-numbingly complicated algorithm for
                % fitting each inter-gap segment separately can be inserted
                % here.  If not necessary (which hopefully will turn out to
                % be the case), then we stick with calling the exact same
                % routine as for gapless trials:
                betteroffset = findBetterOffset(offsetvals);
            else
                % No gaps, no reason to use multiple offsets.
                betteroffset = findBetterOffset(offsetvals);
            end
            if isnan(betteroffset)
                scorecard{t} = 'bad';
                msg = sprintf( ...
                    'Trial %d is hopelessly fudd; evtID & offset sec:', ...
                    t );
                for e = 1:length(testevts)
                    msg = sprintf( '%s\n%10s %.6f', msg, ...
                        dg_thing2str(testevts(e)), offsetvals(e) );
                end
                lfp_log(msg);
                warning('lfp_realignVideoGuts:badtrial', '%s', msg);
                reportstr = sprintf('%s failed', reportstr);
            else
                if abs(betteroffset - prevoffset) <= TStol
                    betteroffset = prevoffset;
                else
                    if isinf(betteroffset) || isnan(betteroffset)
                        error('lfp_realignVideoGuts:badoffset', ...
                            'Bad value calculated for offset on trial %d', ...
                            t );
                    end
                    prevoffset = betteroffset;
                end
                scorecard{t} = 'adjusted with difficulty';
                new_VT_timestamps(trialframeidx(1):end) = ...
                    VT_timestamps(trialframeidx(1):end) + betteroffset;
                reportstr = sprintf('%s %.6f', reportstr, betteroffset);
            end
        else
            % ~any(isbad), so "It's all good!"
            if abs(medoffset - prevoffset) <= TStol
                medoffset = prevoffset;
            else
                if isinf(medoffset) || isnan(medoffset)
                    error('lfp_realignVideoGuts:badoffset2', ...
                        'Bad value calculated for offset on trial %d', ...
                        t );
                end
                prevoffset = medoffset;
            end
            scorecard{t} = 'well-adjusted';
            new_VT_timestamps(trialframeidx(1):end) = ...
                VT_timestamps(trialframeidx(1):end) + medoffset;
            reportstr = sprintf('%s %.6f', reportstr, medoffset);
        end
    end
end
end


function betteroffset = findBetterOffset(offsetvals)
% Attempt to find an offset that makes the majority of events fit within
% maxdiff and the remaining events trigger early.  Returns NaN for
% failure.  Keep liberalizing maxdiff in 10 ms increments starting at 10 ms
% and giving up at 100 ms.  Do not count NaN-valued offsetvals!
betteroffset = NaN;
numneeded = sum(~isnan(offsetvals))/2;
maxoff = max(offsetvals);
for maxdiff = (1:10) * .01
    % 2*maxdiff because we will take the average of maxoff and the closest
    % possible other offset:
    possible = offsetvals > (maxoff - 2*maxdiff);
    if sum(possible) >= numneeded
        possiblevals = sort(offsetvals(possible));
        betteroffset = (possiblevals(1) + possiblevals(end)) / 2;
        return
    end
end
end
