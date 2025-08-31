function [IMFavg, IMFsem, commonstart, composite, freqavg, freqsem, ...
    compositefreq] = dg_makeCompositeIMF(...
    imfs, hht_f, hht_plot, startpt, freqlim, pow)
%[IMFavg, IMFsem, commonstart, composite, ...
%   freqavg, freqsem, compositefreq] = dg_makeCompositeIMF( ...
%   imfs, hht_f, startpt, freqlim, pow)
% Creates a composite IMF by computing a weighted sum of <imfs>, where the
% weighting is determined by a function analogous to the Gaussian but using
% a higher power to give it a flatter top and steeper sides.  Note that the
% first point of each composite is always NaN, because the instantaneous
% frequency at the first sample is be definition NaN.
%
%INPUTS
% mfs: as returned by lfp_hht
% hht_f: as returned by lfp_hht
% startpt: as returned by lfp_hht
% freqlim: a two-element vector specifying the frequency range of interest
%   in Hz.  The weighting function reaches 1/sqrt(e) at both frequencies in
%   <freqlim>.
% pow: the exponent in the "supergauss" function.  This arg is optional,
%   and defaults to 4.
%OUTPUTS
% commonstart: time of first sample that is in common to all jobs in
%   <hht_f>, in samples relative to the alignment ref.
% composite: the sum of the IMFs that have hht_plot{jobidx} set true,
%   weighted according to the "supergauss" function.
% compositefreq: the average of hht_f{jobidx} for the IMFs that have
%   hht_plot{jobidx} set true, weighted according to the "supergauss"
%   function.

%$Rev: 76 $
%$Date: 2010-08-23 20:25:11 -0400 (Mon, 23 Aug 2010) $
%$Author: dgibson $

if nargin < 6
    pow = 4;
end

commonstart = max(startpt);
endpt = [];
for jobidx = 1:length(imfs)
    endpt(jobidx) = length(imfs{jobidx}) + startpt(jobidx) - 1;
end
commonend = min(endpt);
mu = mean(freqlim);
sigma = mu - freqlim(1);
numsamples = commonend - commonstart + 1;
composite = NaN(length(imfs), numsamples);
compositefreq = NaN(length(hht_f), numsamples);
for jobidx = 1:length(imfs)
    thisstartpt = commonstart - startpt(jobidx) + 1;
    samprange = thisstartpt : (thisstartpt + numsamples - 1); 
    IMFweight = exp( ...
        -(hht_f{jobidx}(hht_plot{jobidx},samprange) - mu) .^ pow ...
        / (2 * sigma^pow) );
    composite(jobidx, :) = sum( ...
        imfs{jobidx}(hht_plot{jobidx},samprange) .* IMFweight, 1);
    compositefreq(jobidx, :) = sum( ...
        hht_f{jobidx}(hht_plot{jobidx},samprange) .* IMFweight, 1) ...
        ./ sum(IMFweight, 1);
end
IMFavg = mean(composite, 1);
IMFsem = std(composite, [], 1) / sqrt(length(imfs));
freqavg = mean(compositefreq, 1);
freqsem = std(compositefreq, [], 1) / sqrt(length(hht_f));
