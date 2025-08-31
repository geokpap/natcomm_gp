function [phaseshift, lastpeakpair] = lfp_p2pPhaseChange(...
    filenum, startsample, endsample)
% Returns a boolean marking an "abrupt" change in phase (> 90 degrees)
% between any trio of successive peaks in the HHT waveform between
% startsample and endsample.  The period between the first two peaks is
% used to extrapolate the predicted time of the third peak, and if the
% extrapolated peak falls outside of the bounds of the zero crossings of
% the actual third peak, then a <true> value is returned.  This comparison
% is always carried out in forwards time, and so <startsample> must be less
% than <endsample>.
%INPUTS
% filenum: contains the HHT waveform.
% startsample, endsample: indices into lfp_Samples{filenum}.  They are both
%   assumed to point at peaks, and are thus included in the peaks list
%   <peakrelsamples> (see code).
%OUTPUTS
% phaseshift: true if there is a phaseshift, false otherwise.
% lastpeakpair: The sample indices of the two peaks whose extrapolation
%   revealed a phase shift in the third peak.  If there are fewer than
%   three peaks, meaning that phaseshift undefined, [] is returned for
%   <lastpeakpair> and <phaseshift> is <false>.
%NOTES
% BEFORE checking phaseshift, need to make sure that the maxs don't overlap
% and that there is a zero xing between successive maxs

%$Rev: 186 $
%$Date: 2010-12-22 21:13:51 -0500 (Wed, 22 Dec 2010) $
%$Author: dgibson $

global lfp_Samples;

phaseshift = false;
lastpeakpair = [];
if startsample >= endsample
    error('lfp_p2pPhaseChange:indices', ...
        '<startsample> must be less than <endsample>: %d %d', ...
        startsample, endsample);
end

% find the peaks in the interval:
if lfp_Samples{filenum}(startsample) > 0 %#ok<*USENS>
    peakrelsamples = [1 ...
        dg_findpks(lfp_Samples{filenum}(startsample:endsample)) ...
        endsample - startsample + 1];
else
    % the value at startsample is negative:
    peakrelsamples = [1 ...
        dg_findpks(-lfp_Samples{filenum}(startsample:endsample)) ...
        endsample - startsample + 1];
end
if length(peakrelsamples) < 3
    warning('lfp_p2pPhaseChange:nopeak', ...
        'There are fewer than three peaks in the specified sample range %d:%d.', ...
        startsample, endsample);
    return
end
oppmaxs = sign(lfp_Samples{filenum}(startsample)) * ...
    sign(lfp_Samples{filenum}(endsample)) ~= 1;

for j = 1 : length(peakrelsamples) - 1
    if j == 1
        prevPeak = 1;
    else
        prevPeak = peakrelsamples(j - 1);
    end
    thispeak = peakrelsamples(j);
    % period in samples from previous peak to the current peak:
    period = thispeak - prevPeak;
    % add the period to the current peak:
    extrapPeak = period + thispeak;
    % Check whether the phase changes by more than 90 degrees from one peak in
    % the HHT waveform to the next:
    myrelsample = peakrelsamples(j + 1);
    if extrapPeak > peakrelsamples(j + 1)
        % compare the zeroxing after next peak
        while myrelsample < (numel(lfp_Samples{filenum}) - startsample + 1)
            if (sign(lfp_Samples{filenum}(startsample + myrelsample - 1)) * ...
                    sign(lfp_Samples{filenum}(startsample + myrelsample - 2))) ...
                    ~= 1
                break
            else myrelsample = myrelsample + 1;
            end
        end
        phaseshift = extrapPeak >= myrelsample;
    else % compare to the zeroxing before next peak
        while myrelsample > 1
            if (sign(lfp_Samples{filenum}(startsample + myrelsample - 1)) * ...
                    sign(lfp_Samples{filenum}(startsample + myrelsample - 2))) ...
                    ~= 1
                break
            else myrelsample = myrelsample - 1;
            end
        end
        phaseshift = extrapPeak < myrelsample;
    end
    if phaseshift
        break
    end
end

% now deal with the max at endsample, which is not returned by dg_findpks above:
if ~phaseshift
    if length(peakrelsamples) == 1
        prevPeak = 1;
    else
        prevPeak = peakrelsamples(end - 1);
    end
    thispeak = peakrelsamples(end);
    period = thispeak - prevPeak;
    
    if oppmaxs
        extrapPeak = .5 * period + thispeak;
    else
        extrapPeak = period + thispeak;
    end
    % Check whether the phase changes by more than 90 degrees from one peak in
    % the HHT waveform to the next:
    myrelsample = endsample - startsample + 1;
    if extrapPeak > endsample - startsample + 1;
        % compare the zeroxing after next peak
        while myrelsample < (numel(lfp_Samples{filenum}) - startsample + 1)
            if (sign(lfp_Samples{filenum}(startsample + myrelsample - 1)) * ...
                    sign(lfp_Samples{filenum}(startsample + myrelsample - 2))) ...
                    ~= 1
                break
            else myrelsample = myrelsample + 1;
            end
        end
        phaseshift = extrapPeak >= myrelsample;
    else % compare to the zeroxing before next peak
        while myrelsample > 1
            if (sign(lfp_Samples{filenum}(startsample + myrelsample - 1)) * ...
                    sign(lfp_Samples{filenum}(startsample + myrelsample - 2))) ...
                    ~= 1
                break
            else myrelsample = myrelsample - 1;
            end
        end
        phaseshift = extrapPeak < myrelsample;
    end
end

lastpeakpair = [prevPeak thispeak] + startsample - 1;