function burstbounds = lfp_findBurstBounds( ...
    filenums, bfpts, burstTS, maxerr, phase_opts, verboseflag)
%burstbounds = lfp_findBurstBounds( ...
%    [IMFfilenum errfilenum freqfilenum], bfpts, burstTS, maxerr, ...
%    phase_opts, verboseflag ) 
% Processes a list of candidate bursts to determine the endpoints of each
% candidate and decide which candidates to merge.
%INPUTS
% filenums:  a vector of three lfp_lib filenums specifying the CSC channels
%   to use for:
%       1: composite IMF
%       2: smoothed absolute error trace
%       3: instantaneous frequency trace to be used by lfp_nhtPhaseChange
% bfpts: used by lfp_nhtPhaseChange
% burstTS:  vector containing nominal timestamp of each burst (not
%       necessarily a peak, but must be a point inside the burst)
% maxerr: the maximum value that will be allowed in errfilenum as part of a
%   "good burst"; note that this is evaluated between pairs of peaks, so it
%   is possible that there could be an excessively large error that
%   terminates a burst but is nonetheless included in the last
%   quarter-cycle of the burst.
% phase_opts: a cell array that gets passed directly to lfp_nhtPhaseChange
%   for the purpose of specifying options.  Use {} to invoke
%   lfp_nhtPhaseChange defaults.  Do not specify 'searchbound', as this is
%   hard-coded in this function.
% verboseflag: not all THAT verbose, really, just tells you which burst
%   it's working on and how its endpoints get triggered.
%OUTPUTS
% burstbounds: two columns containing timestamps of the beginning and end
%   of each final "good burst" after merging.  Contains one row for each
%   final "good burst".  "Beginning" and "end" refer to the last sample
%   that DOES belong to the burst.
%
% A "good burst" is considered to extend leftwards and rightwards from a
% peak (extremum) in the composite IMF until any ONE of the following occurs:
%   4. An abrupt phase shift (> 90 degrees over one half cycle) is
%   detected.
%   5. A "zeroxing" failure occurs (i.e., there is a minimum in the
%   composite IMF that has a value >= 0 or a maximum <= 0).
% There are also time-based requirements:  the duration of the burst must
% be positive (not zero or negative); the found bounds must encompass the
% nominal burst peak; and the duration must be at least a fraction of one
% full cycle of the burst frequency as measured by integrating the wave in
% <freqfilenum>, where the minimal fraction is (# samples - 1 / # samples).
%NOTES
% Some data may legitimately contain NaNs.  In a merciful world, there
% would not be any bursts that extend into the regions of NaNs.  However,
% in reality, if there is such a burst, we simply mark its boundary as
% being the NaN that terminated the boundary finding.  For some exact
% statistical purposes, it is thus necessary to look at the actual sample
% value at both boundaries of a burst before including that burst in the
% statistics.


%$Rev: 275 $
%$Date: 2012-05-30 17:54:52 -0400 (Wed, 30 May 2012) $
%$Author: dgibson $

IMFfilenum = filenums(1);
errfilenum = filenums(2);
freqfilenum = filenums(3);

zerodurcount = 0;
shortdurcount = 0;
badburstcount = 0;
closeboundcount = 0;
global lfp_Samples lfp_SamplePeriod
burstbounds = zeros(0, 2);
for burstidx = 1:length(burstTS)
    burstsamp = lfp_time2index(burstTS(burstidx));
    if (burstsamp > 0 && burstsamp < numel(lfp_Samples{IMFfilenum})) ...
            && (lfp_Samples{IMFfilenum}(burstsamp) < 0 && ...
            lfp_Samples{IMFfilenum}(burstsamp) > lfp_Samples{IMFfilenum}(burstsamp-1) ...
            && lfp_Samples{IMFfilenum}(burstsamp) > lfp_Samples{IMFfilenum}(burstsamp+1) ...
            || lfp_Samples{IMFfilenum}(burstsamp) > 0 && ...
            lfp_Samples{IMFfilenum}(burstsamp) < lfp_Samples{IMFfilenum}(burstsamp-1) ...
            && lfp_Samples{IMFfilenum}(burstsamp) < lfp_Samples{IMFfilenum}(burstsamp+1))
        warning('lfp_findBurstBounds:badburstsamp', ...
            'Nominal burstsample %d is zeroxing failure', burstsamp);
        badburstcount = badburstcount + 1;
        continue
    end
    if verboseflag
        fprintf('burstidx=%d of %d; burstsamp=%d\n', ...
            burstidx, length(burstTS), burstsamp);
    end
    errtoobigsamp = NaN(1,2);
    for inc = [-1 1] % search backwards and forwards for burst bounds.
        % Find bounds from fit error:
        errtoobigsamp(1+(inc>0)) = burstsamp + inc;
        while errtoobigsamp(1+(inc>0)) > 1 && ...
                errtoobigsamp(1+(inc>0)) < numel(lfp_Samples{errfilenum})
            if lfp_Samples{errfilenum}(errtoobigsamp(1+(inc>0))) > maxerr || ...
                    isnan(lfp_Samples{errfilenum}(errtoobigsamp(1+(inc>0))))
                break
            else
                errtoobigsamp(1+(inc>0)) = errtoobigsamp(1+(inc>0)) + inc;
            end
        end
        % Note that at this point, <errtoobigsamp> could be off the end of
        % the array by one, so it must not be used directly as an index.
        if inc < 1
            startsamp = errtoobigsamp(1+(inc>0)) - inc;
            endsamp = burstsamp;
        else
            startsamp = burstsamp;
            endsamp = errtoobigsamp(1+(inc>0)) - inc;
        end
        % Tighten bounds to exclude NaNs in IMFfilenum:
        nansample = NaN;
        if any(isnan(lfp_Samples{IMFfilenum}(startsamp:endsamp)))
            if inc < 1
                lastnan = find(isnan(lfp_Samples{IMFfilenum}( ...
                    startsamp:endsamp )), 1, 'last') ...
                    + startsamp - 1;
                thisburst(1) = lastnan;
                nansample = lastnan;
                startsamp = lastnan + 1;
            else
                firstnan = find(isnan(lfp_Samples{IMFfilenum}( ...
                    startsamp:endsamp )), 1, 'first') ...
                    + startsamp - 1;
                thisburst(2) = firstnan;
                nansample = firstnan;
                endsamp = firstnan - 1;
            end
        end
        % Check zeroxings:
        [badzeroxing, firstoffender, lastoffender] = ...
            lfp_zeroxingFailure(IMFfilenum, startsamp, endsamp);
        % boundsamp is the last sample that could be part of the burst:
        if badzeroxing
            if inc > 0
                badzeroxingsamp = firstoffender + startsamp - 1;
            else
                badzeroxingsamp = lastoffender + startsamp - 1;
            end
            boundsamp = badzeroxingsamp - inc;
        else
            if inc > 0
                boundsamp = endsamp;
            else
                boundsamp = startsamp;
            end
        end
        if boundsamp ==  burstsamp
            warning('lfp_findBurstBounds:closebound', ...
                'Burst bound is only one sample from nominal burst sample, skipping lfp_nhtPhasechange');
            closeboundcount = closeboundcount + 1;
            phasesamp(1+(inc>0)) = NaN;
            bf(1+(inc>0)) = NaN;
            relphase(1+(inc>0)) = NaN;
        else
            % Check phase shift:
            [phasesamp(1+(inc>0)), bf(1+(inc>0)), relphase(1+(inc>0))] = ...
                lfp_nhtPhaseChange(freqfilenum, ...
                burstsamp, inc, bfpts, ...
                'searchbound', boundsamp, ...
                phase_opts{:});
        end
        % Note that <phasesamp> can be NaN, and thus must not be used
        % directly as an index.
        if isnan(relphase(1+(inc>0)))
            finalboundsamp = boundsamp;
        else
            finalboundsamp = phasesamp(1+(inc>0));
        end
        if verboseflag
            fprintf(  ' inc:%2d; errtoobig @ %d\n', ...
                inc, errtoobigsamp(1+(inc>0)) - burstsamp);
            if ~isnan(nansample)
                fprintf('  Found NaN at %d\n', nansample - burstsamp);
            end
            if badzeroxing
                fprintf('  badzeroxing @ %d\n', ...
                    badzeroxingsamp - burstsamp);
            end
            if ~isnan(relphase(1+(inc>0)))
                fprintf('   phase shift @ %d, bf=%.1f relphase=%.2f\n', ...
                    phasesamp(1+(inc>0)) - burstsamp, ...
                    bf(1+(inc>0))/lfp_SamplePeriod, relphase(1+(inc>0)));
            end
        end
        % Found boundary of burst.  Trim burst back from the sample that
        % disqualified to the last zeroxing.  We use 'sign', so a sample
        % that is exactly zero is itself a zeroxing if its neighbor is not
        % zero.
        testsamp = finalboundsamp;
        if testsamp > 2 && testsamp < numel(lfp_Samples{IMFfilenum}) - 1 ...
                && sign(lfp_Samples{IMFfilenum}(testsamp)) ~= ...
                sign(lfp_Samples{IMFfilenum}(testsamp + inc))
            % testsamp is the zeroxing
            thisburst(1+(inc>0)) = testsamp; %#ok<*AGROW>
        else
            while testsamp > 1 && testsamp < numel(lfp_Samples{IMFfilenum})
                if sign(lfp_Samples{IMFfilenum}(testsamp)) ~= ...
                        sign(lfp_Samples{IMFfilenum}(testsamp - inc))
                    thisburst(1+(inc>0)) = testsamp - inc;
                    break
                end
                testsamp = testsamp - inc;
            end
        end
    end
    if isnan(thisburst(1)) || isnan(thisburst(2))
        % thisburst(:) is no good because a bound is missing.
        if isnan(thisburst(1))
            boundmissingstr = 'start';
        else
            boundmissingstr = 'end';
        end
        if verboseflag
            fprintf('     bad burst: %d to %d\n', ...
                thisburst(1) - burstsamp, thisburst(2) - burstsamp);
        end
        warning('lfp_findBurstBounds:boundmissing', ...
            'Input burst #%d has no %s', burstidx, boundmissingstr);
        continue
    end
    % Found both boundaries of <thisburst>.  If the burst no longer exists,
    % move on without adding it to burstbounds.  "No longer exists" means
    % either the end is not after the beginning, or the nominal burst peak
    % is not between the beginning and end, or the burst does not contain
    % at least one full cycle of oscillation.
    if thisburst(1) >= thisburst(2) || ...
            (burstsamp < thisburst(1)) || ...
            (burstsamp > thisburst(2))
        if verboseflag
            fprintf('     bad burst: %d to %d\n', ...
                thisburst(1) - burstsamp, thisburst(2) - burstsamp);
        end
        warning('lfp_findBurstBounds:zombie', ...
            'Burst #%d has no duration (%d to %d)', burstidx, ...
            thisburst(1) - burstsamp, thisburst(2) - burstsamp);
        zerodurcount = zerodurcount + 1;
    else
        % The burst might still disqualify by containing less than a full
        % cycle.  Check that.  Because measuring the frequency is somewhat
        % CPU intensive, we assume here that any burst that is at least as
        % long as <bfpts> contains at least one full cycle.
        if thisburst(2) - thisburst(1) + 1 < bfpts
            numcycles = sum(lfp_Samples{freqfilenum}( ...
                thisburst(1):thisburst(2) )) * lfp_SamplePeriod;
            numsamples = thisburst(2) - thisburst(1) + 1;
            minfrac = (numsamples - 1) / numsamples;
            if numcycles < minfrac
                if verboseflag
                    fprintf('     bad burst: %d to %d\n', ...
                        thisburst(1) - burstsamp, thisburst(2) - burstsamp);
                end
                warning('lfp_findBurstBounds:tooshort', ...
                    'Burst #%d has < %.2f cycle duration (%d to %d = %.2f cycle)', ...
                    burstidx, minfrac, ...
                    thisburst(1) - burstsamp, thisburst(2) - burstsamp, ...
                    numcycles);
                shortdurcount = shortdurcount + 1;
                continue % to next putative burst
            end
        end
        % The rest of this 'else' clause is devoted to checking for overlap
        % with previous bursts ("merging") and reporting the results.
        % <burstsamp> points to the latest burst examined so far, but it
        % may begin within the bounds of a previous good burst, in which
        % case they must be merged.  (This could actually be done earlier?)
        thisburstTS = lfp_index2time(thisburst);
        mergedburst = false;
        if size(burstbounds,1) > 0 && ...
                thisburstTS(1) < burstbounds(end, 2)
            mergedburst = true;
            burstbounds(end, 1) = min( ...
                burstbounds(end, 1), thisburstTS(1) );
            burstbounds(end, 2) = max( ...
                burstbounds(end, 2), thisburstTS(2) );
        end
        % Now check preceding good bursts to see if any retrospective
        % merging is required.
        while size(burstbounds,1) > 1
            if burstbounds(end-1, 2) > burstbounds(end, 1)
                % merge to previous good burst, updating burstbounds
                mergedburst = true;
                burstbounds(end - 1, 2) = max( ...
                    burstbounds(end - 1, 2), thisburstTS(2) );
                burstbounds(end, :) = [];
            else
                % All the bursts in <burstbounds> are supposed to be
                % disjoint, so we are done checking for merge.
                break
            end
        end
        if ~mergedburst
            % Found new isolated (thus far) burst.
            burstbounds(end + 1, :) = [thisburstTS(1) ...
                thisburstTS(2)]; %#ok<AGROW>
        end
        if verboseflag
            fprintf('     good burst: %d to %d\n', ...
                thisburst(1) - burstsamp, thisburst(2) - burstsamp);
            fprintf('     merged:%d\n', ...
                mergedburst);
        end
    end
end
fprintf('%d bursts (%.2g percent) had zero duration\n', zerodurcount, ...
    100*zerodurcount/length(burstTS));
fprintf('%d bursts (%.2g percent) had < 1 cycle duration\n', shortdurcount, ...
    100*shortdurcount/length(burstTS));
fprintf('%d bursts (%.2g percent) had zeroxings on the burst sample\n', badburstcount, ...
    100*badburstcount/length(burstTS));
fprintf('%d bursts (%.2g percent) had burst bounds within 1 sample of the burst sample\n', closeboundcount, ...
    100*closeboundcount/length(burstTS));
end

