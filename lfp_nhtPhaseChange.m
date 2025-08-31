function [badsamplenum, bf, relphase] = lfp_nhtPhaseChange(freqfilenum, ...
    startsample, direction, bfpts, varargin)
%[badsamplenum, bf] = lfp_nhtPhaseChange(freqfilenum, ...
%     startsample, direction, bfpts)
% TESTING NOTE: Only the combination of options 'disjoint2', 'adj2freqband'
%   has been carefully tested in the current version (Rev 217).
% Measures phase shift by integrating a raw instantaneous frequency trace
% <freq>. "Zero phase shift" means "constant frequency".  We consider a
% "significant" phase shift to be 0.25 cycle or more within a half cycle of
% oscillation.  We therefore compute a simple moving average baseline
% frequency <bf> by averaging <freq> over 2*<bfpts>+1 initially centered
% on <startsample>.  After converting bf to cycles per sample, the phase
% shift in cycles over any fraction <frac> of a cycle is then equal to 
%       phaseshift = sum(freq(1:reltestpt) - reltestpt * bf
% where
%       reltestpt = round((1/bf)*frac)
% and therefore is defined at the point <reltestpt> relative to
% <startsample>.  (Note that if <bf> changes sufficiently fast, then there
% could be points where <phaseshift> is undefined or multiply defined.
% However, none of this concerns us because we are simply looking for the
% first point where phaseshift exceeds a threshold.) This whole strategy is
% based on the assumption that <bfpts> is sufficiently large to ensure
% that at least one whole cycle of <bf> fits within 2*<bfpts>+1, so there
% is a warning if that assumption fails.  The analysis repeats shifted one
% sample in <direction> until a significant phase shift is found.  Positive
% relative phase means that the actual waveform is oscillating faster than
% baseline, negative relative phase indicates slower than baseline.
%INPUTS
% freqfilenum: filenum containing raw frequency trace in Hz, e.g. from
%   lfp_nhtfreq.
% startsample: the sample number at which to start the search for phase
%   shift.
% direction: +1 to search forward, -1 to search backward.
% bfpts: sets the number of samples over which the baseline frequency is
%   computed (equal-weighted simple moving average of 2*<bfpts>+1
%   samples).  If <bfpts> is 0, then <bfpts> is individually calculated
%   at each offset to be equal to one cycle of the oscillation moving
%   opposite to <direction>, rounded up (i.e. <bfpts> is the smallest
%   number that covers at least one full cycle of oscillation).  If a NaN
%   is encountered during the search for a full cycle, then <bfpts> is the
%   number of points available before hitting the NaN.
%OUTPUTS
% badsamplenum: the sample number of the first sample in <direction> from
%   <startsample> where the phase shift is considered significant.  NaN
%   if there is no such sample, which can happen when we run off the end of
%   lfp_Samples{freqfilenum}.  If we encounter a NaN value in
%   lfp_Samples{freqfilenum} while searching for phase shift, the sample
%   number of that NaN is returned as <badsamplenum>.  Otherwise,
%   <badsamplenum> is the far endpoint of the phase measurement interval
%   (i.e., startsample + offset + direction * reltestpt).
% bf: baseline frequency at the offset where phase shift was detected.  NaN
%   if <bf> was never calculated successfully (i.e. to yield a value that
%   does not trigger warning 'lfp_nhtPhaseChange:lowbf').
% relphase: relative phase (i.e. phase shift) at sample where it exceeded
%   maxshift, in cycles.
%OPTIONS
% 'adj2freqband', freqband - <freqband> is a 3-element vector giving the
%   lower and upper limits for the values in <freqfilenum> that are allowed
%   to go into the <bf> calculation, followed by the filenum to use as the
%   frequency measure.  This could be the same as <freqfilenum>, but
%   usually will be a smoothed version of the same trace.  If the interval
%   options ('disjoint' etc) specify an interval that contains points at
%   which the frequency in channel <freqband(3)> is outside of
%   <freqband(1:2)>, then the bf measurement interval is moved in
%   <direction> until there are no such points.  If there is no set of
%   <bfpts> consecutive points in <freqband(3)> that are in <freqband(1:2)>
%   including <startsample>, then the number of points averaged is reduced
%   as needed so that all points are in <freqband(1:2)>, and a warning is
%   raised.  (This probably works with <bfpts> = 0, but I haven't tested it
%   and there might be undesirable interactions. -DG)
% 'centerbf' - the baseline frequency interval is centered on
%   <startsample> + <offset>, the phase measurement interval extends in
%   <direction> from <startsample> + <offset>.  If <bfpts> is even, then
%   the larger number of points is on the side towards <direction>.
% 'disjoint' - sets the start of the phase measurement interval at the
%   first point beyond the end of the baseline frequency interval. Note
%   that the lack of "surrounding context" may make it more prone to false
%   alarms on frequency sweeps.
% 'disjoint2' - achieves disjointness between baseline frequency interval
%   and phase measurement interval by eliminating the half of the bf period
%   that overlaps with the phase period.  That is, the bf period extends
%   over <bfpts>+1 samples opposite to <direction> starting at
%   <startsample> + <offset>, and the phase measurement interval extends in
%   <direction> from <startsample> + <offset> + <direction>.  Note that the
%   lack of "surrounding context" may make it more prone to false alarms on
%   frequency sweeps.
% 'fixedbf' - does not recalculate <bf> for each moving window (i.e. each
%   value of <offset>), just keeps using the value from the first window
%   (i.e. where <offset> = 0). In the case where bf turns out to be low
%   enough to trigger warning 'lfp_nhtPhaseChange:lowbf', that value is
%   ignored and the first good value (i.e. for the lowest absolute value of
%   offset) calculated is used.
% 'frac', frac - fraction of cycle over which to compute phase 'maxshift',
% 'maxshift', maxshift - maximum phase shift in cycles that is tolerated in
%   the phase measurement interval
% 'offset', initval - uses <initval> as the initial value for <offset>.
%   The default is zero.
% 'searchbound', endsample - Stop searching for phase shift after examining
%   <endsample>.  This means that searching is stopped in any circumstance
%   where a sample beyond endsample would be used for any purpose, e.g.
%   calculating <bf>, in which case it returns badsamplenum = endsample +
%   direction.
% 'verbose' - :-P

%$Rev: 220 $
%$Date: 2011-04-21 20:18:01 -0400 (Thu, 21 Apr 2011) $
%$Author: dgibson $

global lfp_Samples lfp_SamplePeriod

badsamplenum = NaN;
bf = NaN;
relphase = NaN;

centerbfflag = false;
disjointflag = false;
disjoint2flag = false;
endsample = [];
fixedbfflag = false;
frac = 0.5; % fraction of cycle over which to compute phase
freqband = [];
initval = 0;
maxsample = numel(lfp_Samples{freqfilenum});
maxshift = 0.25;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'adj2freqband'
            argnum = argnum + 1;
            freqband = varargin{argnum};
        case 'centerbf'
            centerbfflag = true;
        case 'disjoint'
            disjointflag = true;
        case 'disjoint2'
            disjoint2flag = true;
        case 'fixedbf'
            fixedbfflag = true;
        case 'frac'
            argnum = argnum + 1;
            frac = varargin{argnum};
        case 'maxshift'
            argnum = argnum + 1;
            maxshift = varargin{argnum};
        case 'offset'
            argnum = argnum + 1;
            initval = varargin{argnum};
        case 'searchbound'
            argnum = argnum + 1;
            endsample = varargin{argnum};
        case 'verbose'
            verboseflag = true;
        otherwise
            error('lfp_nhtPhaseChange:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end
if disjointflag + disjoint2flag + centerbfflag > 1
    error('lfp_nhtPhaseChange:disjointtoo', ...
        'Like, dude, make up your mind about ''disjoint''/''disjoint2''/''centerbf''');
end
if abs(direction) ~= 1
    warning('lfp_nhtPhaseChange:direction', ...
        '<direction> should be +1 or -1; behavior for %d is undefined', ...
        direction);
end
if sign(endsample - startsample) ~= sign(direction)
    error('lfp_nhtPhaseChange:searchbound', ...
        '''searchbound'' requires an argument that is on the same side of startsample as <direction>.');
end
offset = initval;
done = false;
while ~done
    % Find <numpts>, the number of points to include in the <bf> average,
    % not counting startsample+offset (for consistency with <bfpts>).
    if bfpts == 0
        nextpt = startsample + offset;
        bfphase = 0;
        numpts = 0; % number of points not counting (startsample + offset)
        % find full cycle
        while nextpt > 1 && nextpt < maxsample && abs(bfphase) < 1
            if isnan(lfp_Samples{freqfilenum}(nextpt))
                break % "while nextpt..." (find full cycle) loop
            else
                bfphase = bfphase + ...
                    lfp_Samples{freqfilenum}(nextpt) * lfp_SamplePeriod;
                numpts = numpts + 1;
                nextpt = nextpt - direction;
            end
        end
        numpts = numpts - 1;  % not to count (startsample + offset)
    else
        numpts = bfpts;
    end
    if startsample + offset - numpts < 1 || ...
            startsample + offset + numpts > maxsample
        % sample index out of bounds
        break
    end
    if (fixedbfflag && isnan(bf)) || ~fixedbfflag
        % compute bf = avg freq in window, in cycles per sample:
        if disjoint2flag && bfpts == 0 && isempty(freqband)
            % bfsamprange is the same as the range we used to calculate
            % numpts, and so this is faster than computing with 'mean':
            bf = bfphase/numpts;
            numptsadj = numpts; %#ok<*NASGU>
        else
            if centerbfflag
                % Find the number of points to use on the +direction side
                % and the -direction side of the center of the bf interval:
                if mod(numpts,2)
                    % numpts is odd, center the interval:
                    numptsplus = (numpts - 1) / 2;
                    numptsminus = (numpts - 1) / 2;
                else
                    % numpts is even, put longer interval on plus side:
                    numptsplus = numpts / 2;
                    numptsminus = numpts / 2 - 1;
                end
            end
            if isempty(freqband)
                % nothing fancy
                numptsadj = numpts;
                offsetadj = offset;
            else
                % First of all, if any of the samples we need to use are
                % beyond <endsample>, we simply bail.
                if ~isempty(endsample) && ( direction > 0 ...
                        && startsample + offset + numpts > endsample ...
                        || direction < 0 ...
                        && startsample + offset - numpts < endsample )
                    badsamplenum = endsample + direction;
                    done = true;
                    relphase = NaN;
                    break    % while ~done
                end
                % Find the number of points on each side of startsample +
                % offset that are within freqband.  The relative index of
                % startsample + offset itself is 1, so the number of points
                % not counting startsample + offset is one less than the
                % index of the last good point, which in turn is one less
                % than the index of the first bad point.  
                ptsonR = startsample + offset : ...
                    startsample + offset + numpts;
                firstBadPtR = find( ...
                    lfp_Samples{freqband(3)}(ptsonR) < freqband(1) | ...
                    lfp_Samples{freqband(3)}(ptsonR) > freqband(2), ...
                    1 );
                % We do not count (startsample + offset) itself (for
                % consistency with bfpts), so it's minus 2:               
                if isempty(firstBadPtR)
                    numgoodptsR = numpts;
                else
                    numgoodptsR = firstBadPtR - 2;
                end
                if numgoodptsR < 0
                    % startsample + offset itself is bad, move on:
                    warning('lfp_nhtPhaseChange:badstart', ...
                        'Skipping sample %d + %d due to freqband', ...
                        startsample, offset);
                    offset = offset + direction;
                    continue % "while ~done" loop
                end
                % we process points on the L in reverse order, so
                % <lastBadPtL> is the number of points to the left of
                % (startsample + offset) where the last bad sample is:
                ptsonL = startsample + offset : -1 : ...
                    startsample + offset - numpts;
                lastBadPtL = find( ...
                    lfp_Samples{freqband(3)}(ptsonL) < freqband(1) | ...
                    lfp_Samples{freqband(3)}(ptsonL) > freqband(2), 1 );
                if isempty(lastBadPtL)
                    numgoodptsL = numpts;
                else
                    numgoodptsL = lastBadPtL - 2;
                end
                if numgoodptsL + numgoodptsR < numpts
                    % We don't have as many good points as we want. Compute
                    % <numptsadj>, the number of good points we actually
                    % have, and give a warning.
                    numptsadj = numgoodptsL + numgoodptsR;
                    warning('lfp_nhtPhaseChange:numpts', ...
                        'Only %d samples are in freqband at sample %d + %d', ...
                        numptsadj, startsample, offset);
                else
                    numptsadj = numpts;
                end
                % Fudge offset as needed to make all points in the <bf>
                % interval be in freqband.  <offsetadj> is the "adjusted"
                % version of <offset>.
                if centerbfflag
                    if direction > 0
                        offsetadj = offset + ...
                            max(0, numptsminus - numgoodptsL);
                    else
                        offsetadj = offset - ...
                            max(0, numptsplus - numgoodptsR);
                    end
                else
                    % All other interval options require numptsadj on the
                    % side of (startsample + offsetadj) opposite to
                    % direction:
                    if direction > 0 && numgoodptsL < numptsadj
                        offsetadj = offset + (numptsadj - numgoodptsL);
                    elseif direction < 0 && numgoodptsR < numptsadj
                        offsetadj = offset - (numptsadj - numgoodptsR);
                    else
                        offsetadj = offset;
                    end
                end
            end
            if offsetadj ~= offset
                warning('lfp_nhtPhaseChange:offsetadj', ...
                    'Had to adjust offset to %d at sample %d + %d to measure bf', ...
                    offsetadj, startsample, offset);
            end
            if disjoint2flag
                bfsamprange = startsample + offsetadj ...
                    - direction * (0:numptsadj);
            elseif centerbfflag
                endptminus = startsample + offsetadj ...
                    - direction * numptsminus;
                endptplus = startsample + offsetadj ...
                    + direction * numptsplus;
                if endptplus > endptminus
                    bfsamprange = endptminus:endptplus;
                else
                    bfsamprange = endptplus:endptminus;
                end
            else
                % measurement intervals are either default or 'disjoint',
                % which both use the same <bf> measurement interval.
                bfsamprange = startsample + offsetadj ...
                    + (-numptsadj:numptsadj);
            end
            % At this point, all points in <bfsamprange> should be within
            % freqband on channel freqband(3), and so <bf> should come out
            % reasonable if freqband(3) is not grossly oversmoothed.  But
            % first, the idiot check: are all points really in freqband,
            % and are none of the samples in <bfsamprange> beyond
            % <endsample>?
            if ~isempty(endsample) && ( any(lfp_Samples{freqband(3) ...
                    }(bfsamprange) < freqband(1) | ...
                    lfp_Samples{freqband(3)}(bfsamprange) > freqband(2)) ...
                    || direction > 0 && any(bfsamprange > endsample) ...
                    || direction < 0 && any(bfsamprange < endsample) )
                warning('lfp_nhtPhaseChange:oops', ...
                    'This should not have happened.  But I''m dealing with it anyway.');
                badsamplenum = endsample + direction;
                done = true;
                relphase = NaN;
                break    % while ~done
            end
            bf = mean(lfp_Samples{freqfilenum}(bfsamprange)) ...
                * lfp_SamplePeriod;
        end
    end
    % compute relative phase at reltestpt, which is frac cycle of bf away
    % from (startsample + offset):
    reltestpt = round(frac/bf);
    if reltestpt == 0
        error('Burma!');
    end
    samplerange = startsample + offset + direction * (1:reltestpt);
    if direction > 0
        firstsamp = samplerange(1);
        lastsamp = samplerange(end);
    else
        firstsamp = samplerange(end);
        lastsamp = samplerange(1);
    end
    if disjointflag
        % shift samplerange
        samplerange = samplerange + direction * numpts;
    end
    if direction > 0
        if samplerange(1) < 1 || samplerange(end) > maxsample
            % sample index out of bounds
            break
        end
    else
        if samplerange(end) < 1 || samplerange(1) > maxsample
            % sample index out of bounds
            break
        end
    end
    % check for NaNs, check 'searchbound'; if neither of those then
    % calculate relphase:
    if any(isnan(samplerange)) || ...
            any(isnan(lfp_Samples{freqfilenum}(samplerange)))
        % Find the NaN closest to <startsample> in <direction>
        badsamplenum = startsample;
        while ~isnan(lfp_Samples{freqfilenum}(badsamplenum))
            badsamplenum = badsamplenum + direction;
        end
        done = true;
        relphase = NaN;
    elseif ~isempty(endsample) && ...
            endsample + direction >= firstsamp && ...
            endsample + direction <= lastsamp
        % endsample is within samplerange
        badsamplenum = endsample + direction;
        done = true;
        relphase = NaN;
    else
        % relphase is the difference between actual and extrapolated:
        relphase = sum( lfp_Samples{freqfilenum}(samplerange) ...
            ) * lfp_SamplePeriod - reltestpt * bf;
        if abs(relphase) > maxshift
            done = true;
            badsamplenum = ...
                startsample + offset + direction * reltestpt;
        end
    end
    if verboseflag
        fprintf('offset:%02d; bf=%.2fHz; rel samplerange %d:%d; phase %8.3f\n', ...
            offset, bf/lfp_SamplePeriod, samplerange(1)-startsample, ...
            samplerange(end)-startsample, relphase)
    end
    offset = offset + direction;
end
% At this point, badsamplenum either was found on the last iteration, or
% still has its initial value of NaN.  There is nothing left to do.

        