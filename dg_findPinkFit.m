function [pinkBL, cfun, f1] = dg_findPinkFit(BL, f1, f2, fstep, errthresh)
%INPUTS
% BL: as returned by lfp_BLspectrum; see dg_fitPinkSpec.
% f1: initial lower limit of fit.
% f2: static upper limit of fit.
% fstep: amount by which to decrement <f1> on each iteration.
% errthresh: quitting criterion; when the mean squared error between the
%   logarithm of the fitted spectrum and the logarithm of BL.sum increases
%   by more than <errthresh> times the mean squared error of the first fit,
%   the iteration has gone too far and the previous values are returned.
%OUTPUTS
% pinkBL: fitted pink spectrum; see dg_fitPinkSpec.
% cfun: parameters of fitted pink spectrum; see dg_fitPinkSpec.
% f1: the last (lowest) value of <f1> before termination of the iteration.
%NOTES
% Starts with the result of dg_fitPinkSpec(...'logfit', 'freqlim', [f1 f2])
% and progressively lowers <f1> by <fstep> until the mean squared error
% between the log of pinkBL and the log of BL increases by more than
% <errthresh> times the mean squared error on the first fit.  Since the
% size of the increase in error will usually be larger when the frequency
% steps are larger, <fstep>, <errthresh>, and the width of the initial band
% <f2>-<f1> all interact with each other. However, due to the fact that LFP
% spectra always have bands that are much more powerful than the best-fit
% pink spectrum, surprisingly large values of <errthresh>  can be required
% to make the fit cover what looks "reasonably pink" to the eye.  In a case
% like this:
%   dg_findPinkFit(BL, 5, 12, 0.1, 6);
% where the step size is only 70 times smaller than the initial bandwidth,
% a very high <errthresh> may be needed.  On the other hand, this:
%   dg_findPinkFit(BL, 100, 400, 1, 0.1);
% does a pretty good job of finding the upper edge of the gamma band and
% fitting the frequencies above it.  Note that the step size is 300 times
% smaller than the initial bandwidth, so much smaller <errthresh> values
% perform well.
%   The Grand Caveat: "It all depends on what you're trying to do." -DG

%$Rev: 210 $
%$Date: 2015-02-09 22:45:59 -0500 (Mon, 09 Feb 2015) $
%$Author: dgibson $

prev_mse = [];
while f1 > 0
    [newpinkBL, newcfun] = dg_fitPinkSpec(BL, 'logfit', 'freqlim', [f1 f2]);
    inrange = BL.f >= f1 & BL.f < f2;
    logmse = mean((log10(BL.sum(inrange)/BL.N) - ...
        log10(newpinkBL.sum(inrange))).^2);
    if isempty(prev_mse) || logmse - prev_mse < errthresh * mse1;
        if isempty(prev_mse)
            mse1 = logmse;
        end
        prev_mse = logmse;
        pinkBL = newpinkBL;
        cfun = newcfun;
        prev_f1 = f1;
        f1 = max(0, f1 - fstep);
    else
        break
    end
end
if logmse - prev_mse < errthresh * mse1
    f1 = prev_f1;
end
