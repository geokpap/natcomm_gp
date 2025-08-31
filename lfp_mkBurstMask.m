function mask = lfp_mkBurstMask(burstbounds)

% Makes a mask containing <true> for each sample that is part of a burst
% and <false> for all other samples.
%INPUTS
% burstbounds: cell contents of one element of <burstbounds> as returned by
%   lfp_findBurstBounds, containing timestamps of the start and end of each
%   burst.
% filenum: a CSC channel containing the data from which <burstbounds> was
%   calculated.
%OUTPUTS
% mask: a logical array the same size as lfp_Samples{lfp_ActiveFilenums(1)}
%   which is true for samples that belong to bursts and false elsewhere.
%   Bursts are defined as extending from the last zero crossing or
%   "zeroxing failure" (i.e. a local maximum whose value is negative or a
%   local minimum whose value is positive) before the first peak in the
%   burst (as listed in <burstbounds>) to the first zeroxing or zeroxing
%   failure after the last peak in the burst.  If there is no such zero
%   crossing, then the burst extends to the beginning of the data or to the
%   end of the data.  <mask> is suitable for use with lfp_wavemask.
%NOTES
%   Example usage:
% sessiondir = 'HH092807';
% fname = 'e19_hht_13-20_hz.mat';
% implant_str = 'HH0407';
% brain_structure_str = 'CN'; 
% [bursts, goodburstsidx, burstbounds] = hhtBurstDG( ...
%     sessiondir, fname, implant_str, brain_structure_str);
% fragnum = 1;
% mask = mkBurstMask(bursts{fragnum}(burstbounds{fragnum}, 1);
% lfp_createWave(@lfp_wavemask, 1, mask);
% lfp_disp(100, [], [], 'avg', 'ovr')


%$Rev: 303 $
%$Date: 2013-08-26 18:51:03 -0400 (Mon, 26 Aug 2013) $
%$Author: dgibson $

global lfp_Samples lfp_ActiveFilenums

firstTS = lfp_index2time(1);
lastTS = lfp_index2time(numel(lfp_Samples{lfp_ActiveFilenums(1)}));
mask = false(size(lfp_Samples{lfp_ActiveFilenums(1)}));
for burstnum = 1:size(burstbounds, 1)
    if burstbounds(burstnum, 1) > lastTS || ...
            burstbounds(burstnum, 2) < firstTS
        continue
    end
    mask( lfp_time2index(burstbounds(burstnum, 1)) : ...
        lfp_time2index(burstbounds(burstnum, 2)) ) = true;
end
end

