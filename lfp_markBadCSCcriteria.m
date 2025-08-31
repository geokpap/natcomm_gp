function [isbadfile, isdeadfile, allcriteria] = lfp_markBadCSCcriteria(...
    medians, clipping, fmax, numbadpeaks, choppiness, LFpwr, ...
    binwidths, signalSD, bincounts, winwidth)
%NOTES
% <binthresh> has arbitrarily been made symmetrical by starting with the
% actual minimum and maximum normcounts from the good files in the
% hand-classified sessions in
% 2010-11-18_18-28-38_+_2011-07-12_19-30-52down16.mat and then narrowing
% the "non-fishy" region for the bin where it was wider.


%$Rev: 343 $
%$Date: 2015-04-05 16:41:55 -0400 (Sun, 05 Apr 2015) $
%$Author: dgibson $

%% Thresholds of badness
badmedianmax = 1e-4;
badclippingmax = 2;
badfmaxmax = 15;
badchoppinessmax = 11;
badLFpwrmax = 37;
badbinwidthmax = 2e-4;
deadfilemin = 1e-6;
badsignalSDmax = 3e-4;
badbincountmax = [NaN
    0.15
    0.30
    0.69
    1
    1.17
    0.69
    0.30
    0.15
    NaN];
badbincountmin = [NaN
    0.05
    0.16
    0.40
    1
    0.90
    0.40
    0.16
    0.05
    NaN];

%%

binnums = [2 3 4 6 7 8 9];
fishybins = false(length(medians), length(binnums));
for binidx = 1:length(binnums)
    bn = binnums(binidx);
    binthresh(binidx,:) = [badbincountmin(bn) badbincountmax(bn)]; %#ok<AGROW>
end
normcounts = bincounts ./ repmat(bincounts(:,5), 1, size(bincounts,2));
for binidx = 1:length(binnums)
    fishybins(:,binidx) = normcounts(:, binnums(binidx)) < binthresh(binidx,1) ...
        | normcounts(:, binnums(binidx)) > binthresh(binidx,2);
end

% Thresholds of badness
isbadmedian = abs(medians) > badmedianmax;
isbadclipping = clipping > badclippingmax;
isbadfmax = fmax > max(badfmaxmax, 3/winwidth);
isbadchoppiness = choppiness > badchoppinessmax;
isbadLFpwr = LFpwr > badLFpwrmax;
isbadbinwidth = binwidths > badbinwidthmax;
isbadsignalSD = signalSD > badsignalSDmax;
isbadbincount = sum(fishybins, 2) > 5;

% Thresholds of suspiciousness
isdeadfile = signalSD < deadfilemin;
suspiciousness = false(length(medians), 10);
suspiciousness(:,2) = clipping > 0.5;
suspiciousness(:,4) = choppiness > 7;
suspiciousness(:,5) = LFpwr > 34;
suspiciousness(:,9) = sum(fishybins, 2) > 3;
suspiciousness(:,10) = isdeadfile;

% The final answers
isbadfile = isbadmedian | isbadclipping | isbadfmax | isbadchoppiness ...
    | isbadLFpwr | isbadbinwidth | isbadsignalSD ...
    | isbadbincount | numbadpeaks > 0 | sum(suspiciousness, 2) > 1;
allcriteria = [isbadmedian  isbadclipping  isbadfmax  isbadchoppiness ...
    isbadLFpwr isbadbinwidth  isbadsignalSD  isbadbincount ...
    numbadpeaks sum(suspiciousness, 2) ];

