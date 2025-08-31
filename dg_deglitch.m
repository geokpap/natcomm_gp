function [samples, numglitches] = dg_deglitch(samples, thresh)
%DG_DEGLITCH removes single-sample spikes from sampled data
%samples, numglitches = dg_deglitch(samples, thresh)
%  Returns the same sample data that was passed in, modified by
%  replacing single-sample glitches with the average of the two surrounding
%  samples.  A single-sample glitch is defined as a sample that differs
%  from its predecessor by more than <thresh>, but which is followed by a
%  sample that differs by less than <thresh>.  <numglitches> is the number
%  of glitches that were removed.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

dif1 = samples(2:end-1) - samples(1:end-2);
dif2 = samples(3:end) - samples(1:end-2);

glitches = find(abs(dif1) > thresh & abs(dif2) < thresh);
numglitches = numel(glitches);
samples(glitches+1) = (samples(glitches) + samples(glitches+2))/2;