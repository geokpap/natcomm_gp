function samples = dg_rmFSCVart(samples, artidx, nbefore, nafter)
%INPUTS
% samples: array of samples.
% artidx: index into <samples> of all samples that have been identified
%   with an instance of an FSCV artifact, e.g. validated threshold crossing
%   samples.
% nbefore: number of samples to interpolate before the index sample.
% nafter: number of samples to interpolate after the index sample.
%OUTPUTS
% samples: a copy of the input <samples> with FSCV artifacts interpolate
%   away.
%NOTES
% This is a trivial addendum to dg_fitFSCVart.

%$Rev: 276 $
%$Date: 2021-08-06 22:39:12 -0400 (Fri, 06 Aug 2021) $
%$Author: dgibson $

for artnum = 1:length(artidx)
    idx = (-nbefore:nafter) + artidx(artnum);
    samples(idx) = linspace( ...
        samples(idx(1)), samples(idx(end)), length(idx) );
end

    