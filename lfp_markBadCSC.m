function [filenames, isbadfile, isdeadfile] = lfp_markBadCSC( ...
    sessiondir, varargin)
%INPUT
% sessiondir: a string specifying the relative or absolute pathname to a
%   directory containing CSC files and an events file to scan for bad CSC
%   files.
%OUTPUT
% filenames: the list of files in <sessiondir> that were analyzed.
% isbadfile: logical array of same size as <filenames>; true for files that
%   are probably bad and thus should not be included in local average
%   references.
% isdeadfile: logical array of same size as <filenames>; true for files
%   that might be from dead channels.  These are not marked <true> in
%   <isbadfile> because they have much higher chances of actually being
%   good channels than the ones marked in <isbadfile>.
% allcriteria: logical array of same number of rows as <filenames>; columns
%   are: isbadmedian isbadclipping isbadLFpwr isbadbinwidth isbadsignalSD
%   isbadbincount.  A simple graded score of how bad each file is can be
%   computed by counting the number of reasons for which it was marked bad,
%   e.g.:
%       numreasons = sum(allcriteria, 2)
% Output is also saved to the file 'lfp_markBadCSC.mat' in <sessiondir>.
% 'lfp_markBadCSC.mat' contains the following variables:
%   filenames, isbadfile, isdeadfile: as returned.
%   medians: col vector of doubles containing median for each file.
%   binwidths: the difference between the median and 25 prctile for each
%       file.  This is used as the binwidth for a low-resolution histogram
%       of sample values.
%   bincounts: from -5 * binwidth to 5 * binwidth, in files x bins format,
%       meaning it has 10 columns (zero is a bin edge, not a bin center).
%   clipping: an estimate of how severe the clipping is.  Note there may be
%       a lot of clipping during other task times (e.g. ITIs) that do not
%       matter to a particular analysis, in which case these files can
%       still be analyzed even if they are marked as containing excessive
%       levels of clipping.
%   signalSD: if a middle peak is found in the Gauss fitting, then this is
%     the sigma of that peak; otherwise, it is the standard deviation of
%     the samples that have values in <ismid>, which is "middle range"
%     determined by various heursitcis (see code). 
%   fmax: the frequency of absolute maximum power in the spectrum (see code
%       beginning with "BL = lfp_BLspectrum" for parameters of spectrum
%       computation).
%   numbadpeaks: spectral peaks are "bad" (i.e. most likely caused by
%       electrical noise) if the width of the peak is less than 1/10 of the
%       peak frequency, the power integrated over the width of the peak is
%       more than 80% of the power integrated from 5 Hz to 30 Hz, and the
%       frequency of the peak is between 5 and 100 Hz.
%   choppiness: an index of the degree to which the spectrum swings between
%       high and low values, rather than being smooth (as is typical of
%       clean LFP spectra).
%   LFpwr: the maximum power below 5 Hz relative to a pink (power-law)
%       spectrum fitted from 5 to 55 Hz that is extrapolated linearly to
%       lower frequencies.
%   BLs: the raw power spectrum for each file, computed by lfp_BLspectrum
%       (see code beginning with "BL = lfp_BLspectrum" for parameters of
%       spectrum computation).
%   freqlim: the frequency band over which most spectral analysis is done,
%       5 - 100 Hz by default.
%   winwidth: window width for lfp_BLspectrum.
%OPTIONS
% 'destdir', destdir - saves output file to <destdir> with the session ID
%	coded into the filename like lfp_markbadCSC_<sessionID>.mat where
%	<sessionID> is the final directory name in <sessiondir>.  (This is
%	mainly meant for testing and debugging purposes.)
% 'eventfile', funchandle - to handle the case (e.g. Simon) where each CSC
%   file has its own separate events file, or for some other reason
%   dg_findEventsFile will not be able to find the right events file for a
%   given CSC file, <funchandle> provides a way to specify a function that
%   will find the right events file.  <funchandle> must return a string
%   containing the event file name (without directory), and must accept an
%   argument whose value will be the CSC filename (but of course
%   <funchandle> is free to do nothing with that value).
% 'freqlim', freqlim - <freqlim> is a two-element numeric vector
%   specifying the lower limit and upper limit of the frequency range in
%   Hertz to fit with a pink spectrum to compute LFpwr.  Default [5 100].
%   Note that use of this option likely requires making changes to
%   lfp_markBadCSCcriteria, which does not have any corresponding option,
%   and therefore use of this option is NOT RECOMMENDED.
% 'his' - only compute sample-value-histogram-based criteria.
% 'matchExpr', matchExpr - specify a different regular expression to use
%   for finding files to analyze.  Default is 'csc\d+_down\d+.mat'.  The
%   matching is done case-insensitively.  As for dg_localAvgRef,
%   <matchExpr> may alternatively be a cell string array containing the
%   literal filenames including extensions.
% 'spec' - only compute spectrum-based criteria.
% 'winwidth', winwidth - <winwidth> is a numeric scalar that specifies the
%   width in seconds of the time window used for computing BLs.  Default 3.
%   Note that use of this option likely requires making changes to
%   lfp_markBadCSCcriteria, which does not have any corresponding option,
%   and is NOT RECOMMENDED. See NOTES for caveats.
%NOTES
%   isdeadfile:
% The (signalSD < 3e-5) criterion produces some false alarms, i.e. files
% that have low signal amplitudes but nonetheless look physiological in the
% spectrograms.  So someday those should be analyzed spectrogrammatically
% rather than just calling them all bad.  At this point, they are just
% returned separately in <isdeadfile>.
%   Criteria for isbadfile:
% See lfp_markBadCSCcriteria.  Parameters that are tested include median
% value, width of second quartile, clipping, spectral "choppiness",
% narrowness of spectral peaks, frequency of maximum power (note that this
% can raise false alarms in very strongly oscillatory data at frequencies
% above 15 Hz or 3/winwidth, e.g. parkinsonian beta), low frequency power,
% variance, overall distribution of sample values.
%   Choppiness:
% Note that this measurement does not even attempt to take into account the
% effects of changing winwidth, so that makes the winwidth option pretty
% useless.  In fact all the spectral measure thresholds do.
%   Spectral Peaks
% Physiological peaks tend to be rather broad, whereas artifactual peaks
% from electrical equipment tend to be narrow.  The ratio of a peak's
% half-power bandwidth to its center frequency is used the measure of peak
% width, and if it is less than 0.1 it counts as a "bad peak".  Peaks with
% maxima at less than one hundredth the maximum of the whole spectrum (i.e.
% less than -20 dB) are not counted.

%$Rev: 391 $ $Date: 2017-08-02 14:38:47 -0400 (Wed, 02 Aug 2017) $ $Author:
%dgibson $

global lfp_NominalTrialStart lfp_NominalTrialEnd lfp_Samples ...
    lfp_ActiveFilenums lfp_FreqLim

winwidth = 3;
freqlim = [5 100];
han5 = hanning(5);
maskthresh = 0.02;

% Flags are: his, spec
destdir = '';
evtfunc = [];
flags = true(2,1); % his
matchExpr = 'csc\d+_down\d+.mat';
argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'destdir'
                argnum = argnum + 1;
                destdir = varargin{argnum};
            case 'eventfile'
                argnum = argnum + 1;
                evtfunc = varargin{argnum};
            case 'freqlim'
                argnum = argnum + 1;
                freqlim = varargin{argnum};
            case 'his'
                flags(:) = false;
                flags(1) = true;
            case 'matchExpr'
                argnum = argnum + 1;
                matchExpr = varargin{argnum};
            case 'spec'
                flags(:) = false;
                flags(2) = true;
            case 'winwidth'
                argnum = argnum + 1;
                winwidth = varargin{argnum};
            otherwise
                error('lfp_markBadCSC:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('lfp_markBadCSC:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if isempty(evtfunc)
    evtfname = dg_findEventsFile(sessiondir);
end
files = dir(sessiondir);
filenames = reshape({files(~cell2mat({files.isdir})).name}, [], 1);
if iscell(matchExpr)
    fnidx = find(ismember(filenames, matchExpr));
else
    fnidx = find(~cellfun(@isempty, regexpi(filenames, matchExpr)));
end
if isempty(fnidx)
    error('lfp_markBadCSC:nomatch', ...
        'There is no match in %s for ''%s''', ...
        sessiondir, dg_thing2str(matchExpr));
end
filenames = filenames(fnidx);

% Initialize all values used for thresholding so the formulas won't crash
% if a computation gets skipped:
medians = NaN(length(filenames),1);
binwidths = NaN(length(filenames),1);
bincounts = NaN(length(filenames),10);
clipping = NaN(length(filenames),1);
signalSD = NaN(length(filenames),1);
LFpwr = NaN(length(filenames),1);
fmax = NaN(length(filenames),1);
choppiness = NaN(length(filenames),1);
numbadpeaks = NaN(length(filenames),1);
BLs = [];

for fidx = 1:length(filenames)
    if ~isempty(evtfunc)
        evtfname = feval(evtfunc, filenames{fidx});
    end
    lfp_read2('preset', sessiondir, {evtfname, filenames{fidx}});
    fn = lfp_ActiveFilenums(1);
    if flags(1)
        % Coarse median-centered histogram for distribution shape; start by
        % finding first 2 quartiles.
        prcts = prctile(lfp_Samples{fn}(:), [25 50]);
        medians(fidx) = prcts(2);
        binwidths(fidx) = medians(fidx) - prcts(1);
        counts = histc(lfp_Samples{fn}(:), ...
            (-5:5)*binwidths(fidx) + medians(fidx));
        bincounts(fidx,:) = counts(1:end-1);
        % Fine full-range histogram for clipping detection. <scale> is the
        % actual range of the data, exapanded by 5%:
        rawscale = [min(lfp_Samples{fn}(:)), max(lfp_Samples{fn}(:))];
        rawscalerange = diff(rawscale);
        scale(1) = rawscale(1) - 0.02 * rawscalerange;
        scale(2) = rawscale(2) + 0.02 * rawscalerange;
        numbins = 500;
        finebinwidth = diff(scale)/numbins;
        binctrs = ((0:numbins)' * diff(scale)) / numbins + scale(1);
        n = reshape(hist(lfp_Samples{fn}(:), binctrs), [], 1);
        % smooth the histogram slightly:
        n = conv(n, han5, 'same');
        % mask out "long tails" if there are any:
        nmasked = n;
        [maxn, maxbin] = max(n);
        nmasked(n < maxn/100) = 0;
        
        % Begin the ordeal of trying to fit the histogram.  Rather than
        % attempt to control 'fit' with 'fitoptions' (which is difficult),
        % we just run the default fitter piecewise. First we find the range
        % of bins that encompass the 60% points to the left and to the
        % right of the tallest peak.
        leftidx = find(nmasked(1:maxbin) < 0.6 * maxn, 1, 'last');
        rightidx = find(nmasked(maxbin:end) < 0.6 * maxn, 1) + maxbin - 1;
        thisrange = max(leftidx, 1) : min(rightidx, length(nmasked));
        residual = nmasked;
        newpeaknum = 1;
        emptyranges = {[binctrs(1), binctrs(end) + finebinwidth]};
        rangeidx = 1; % index into <emptyranges>
        % Fit first gaussian and compute residual
        try
            fitobject = fit(binctrs(thisrange), ...
                residual(thisrange), 'gauss1');
            c = coeffvalues(fitobject); 
        catch
            c = [];
        end
        if isempty(c)
            % forget the whole curve fitting concept and compute simple
            % statistics
            signalSD(fidx) = std(lfp_Samples{fn}(:));
            clipping(fidx) = 0;
        else
            % compute the fitted function and the residual:
            amp = c(1);
            mu = c(2);
            sigma = c(3) / sqrt(2);
            Y = amp*exp(-((binctrs-mu)/c(3)).^2);
            residual = nmasked - Y;
            % re-mask the residual:
            residual(residual < maxn * maskthresh) = 0;
            % update <emptyranges> list:
            mumin = max( binctrs(1), ...
                mu(newpeaknum) - max(finebinwidth, sigma(newpeaknum)) );
            mumax = min( ...
                mu(newpeaknum) + max(finebinwidth, sigma(newpeaknum)), ...
                binctrs(end) + finebinwidth );
            emptyranges{end+1} = [emptyranges{rangeidx}(1) mumin]; %#ok<AGROW>
            emptyranges{rangeidx} = [mumax emptyranges{rangeidx}(2)];
            % Put the remaining peak(s) in the longest empty range(s),
            % provided there is enough residual in the range to be worth
            % fitting (i.e. at least 5% of <maxn>).
            newpeaknum = newpeaknum + 1;
            while newpeaknum <= 9
                if isempty(emptyranges)
                    % There are no unfitted ranges left to fit, so stop
                    % fitting.
                    break
                end
                rangelengths = cellfun(@diff, emptyranges) + 1;
                [~, rangeidx] = max(rangelengths);
                thisrange = binctrs >= emptyranges{rangeidx}(1) & ...
                    binctrs < emptyranges{rangeidx}(2);
                if max(residual(thisrange)) / maxn <= maskthresh
                    % Forget that emptyrange forever, try next
                    emptyranges(rangeidx) = [];
                    continue
                else
                    try
                        fitobject = fit(binctrs(thisrange), ...
                            residual(thisrange), 'gauss1');
                        c = coeffvalues(fitobject);
                    catch
                        % Forget that emptyrange forever, try next
                        emptyranges(rangeidx) = [];
                        continue
                    end
                    amp(newpeaknum) = c(1);
                    mu(newpeaknum) = c(2);
                    sigma(newpeaknum) = c(3) / sqrt(2);
                    Y = amp(newpeaknum) * exp( -( ...
                        (binctrs-mu(newpeaknum)) / c(3) ).^2 );
                    residual = residual - Y;
                    residual(residual < maxn * maskthresh) = 0;
                    % update <emptyranges> list:
                    mumin = max( binctrs(1), ...
                        mu(newpeaknum) - max( finebinwidth, ...
                        sigma(newpeaknum) ) );
                    mumax = min( ...
                        mu(newpeaknum) + max(finebinwidth, sigma(newpeaknum)), ...
                        binctrs(end) + finebinwidth );
                    emptyranges{end+1} = [emptyranges{rangeidx}(1) mumin]; %#ok<AGROW>
                    emptyranges{rangeidx} = [mumax emptyranges{rangeidx}(2)];
                    newpeaknum = newpeaknum + 1;
                end
            end
            % We don't want any peaks that are outside the range of
            % <scale>, which 'fit' is perfectly capable of generating.
            isbogus = mu < scale(1) | mu > scale(2);
            amp(isbogus) = [];
            mu(isbogus) = [];
            sigma(isbogus) = [];
            % Find <islow>, <ismid>, and <ishigh>, which are logical
            % indices into <binctrs>.  <ismid> should cover the normal
            % signal range for the channel.  The other two should cover the
            % clipping peaks if there are any, or the tails of the expected
            % normalish distribution if there aren't any clipping peaks.
            % (This definition should produce slightly larger than normal
            % values of <clipping> for signals that have "up and down
            % states".)
            %   Also, there can "real peaks", which produce minima in
            % between, and "bogus peaks", which don't.  If there are no
            % minima at all, then we just want the "biggest" peak, which I
            % hereby define as the one that has the largest product of
            % amplitude and sigma. If there are minima, life gets
            % complicated.
            peaksize = amp .* sigma;
            if length(amp) == 1
                hasminimum = false;
            else
                fittedcurve = amp(1) * exp( -((binctrs-mu(1)) / ...
                    (sqrt(2)*sigma(1)) ).^2 );
                for k = 2:length(amp)
                    fittedcurve = fittedcurve + ...
                        amp(k) * exp( -((binctrs-mu(k)) / ...
                        (sqrt(2)*sigma(k)) ).^2 );
                end
                isrepeat = fittedcurve(2:end) == fittedcurve(1:end-1);
                % there could potentially be long strings of zeros that
                % would mess up the minimum-finding algorithm:
                nnorepeat = fittedcurve(~isrepeat);
                binctrsnorepeat = binctrs(~isrepeat);
                isincreasing = nnorepeat(2:end) > nnorepeat(1:end-1);
                hasminimum = any( ~isincreasing(1:end-1) ...
                    & isincreasing(2:end) );
            end
            if hasminimum
                % <minidx> points to the last bin before the minimum:
                minidx = find( ~isincreasing(1:end-1) & ...
                    isincreasing(2:end) );
                % Each minimum has a maximum on each side.  Each maximum
                % can contain only one "real peak".  If a maximum contains
                % more than one <mu>, then we must throw away the "smaller"
                % sized peak.  We process the left maximum at each minimum
                % until we get to the last, where we also process the right
                % maximum.
                isbogus = false(size(amp));
                for minnum = 1:length(minidx)
                    if minnum == 1
                        leftbound = 1;
                    else
                        leftbound = minidx(minnum - 1) + 1;
                    end
                    rightbound = minidx(minnum);
                    peakidxinmax = find( mu >= binctrsnorepeat(leftbound) ...
                        & mu < binctrsnorepeat(rightbound) );
                    if length(peakidxinmax) > 1
                        [~, biggest] = max(peaksize(peakidxinmax));
                        isbogus(setdiff(peakidxinmax, biggest)) = true;
                    end
                    if minnum == length(minidx)
                        leftbound = minidx(minnum);
                        rightbound = length(binctrsnorepeat);
                        peakidxinmax = find( ...
                            mu >= binctrsnorepeat(leftbound) ...
                            & mu < binctrsnorepeat(rightbound) );
                        if length(peakidxinmax) > 1
                            [~, biggest] = max(peaksize(peakidxinmax));
                            isbogus(setdiff(peakidxinmax, biggest)) = true;
                        end
                    end
                end
                amp(isbogus) = [];
                mu(isbogus) = [];
                sigma(isbogus) = [];
            else
                [~, biggest] = max(peaksize);
                amp = amp(biggest);
                mu = mu(biggest);
                sigma = sigma(biggest);
            end
            % We now have now have a fit which sums one to three gaussians
            % that all correspond to "real peaks".
            numpeaks = length(amp);
            MIDidx = [];
            switch(numpeaks)
                case 0
                    % hopeless, do nothing
                case 1
                    % We define the single peak as the "middle" one:
                    MIDidx = 1;
                    % Use the single fitted peak for signalSD:
                    signalSD(fidx) = sigma;
                    % Set <islow>, <ismid>, and <ishigh> using fitted
                    % parameters:
                    islow = binctrs >= scale(1) & binctrs < ...
                        mu - 2*sigma;
                    ismid = binctrs >= mu - 2*sigma & ...
                        binctrs < mu + 2*sigma;
                    ishigh = binctrs >= mu + 2*sigma;
                case {2 3 4 5 6 7 8 9}
                    % <Lidx> and <Ridx> point to the "clipping peaks":
                    [~, Lidx] = min(mu);
                    [~, Ridx] = max(mu);
                    % Clipping peaks do not remotely include zero, whereas
                    % the "middle" peak should include zero (even if there
                    % is a lot of DC offset).  At five sigma from the mean,
                    % 'normcdf' is p<1e-6, so we put the dividing line at
                    % 10 sigma just to be liberal, and re-assign peaks
                    % that qualify as "middle".
                    numsigmas = abs(mu) ./ sigma;
                    if numsigmas(Lidx) < 10
                        MIDidx = Lidx;
                        Lidx = [];
                    end
                    if numsigmas(Ridx) < 10
                        MIDidx = Ridx;
                        Ridx = [];
                    end
                    % It is now logically impossible to have empty values
                    % in both MIDidx *and* Lidx or Ridx.  We set <islow>,
                    % <ismid>, and <ishigh> according to which peaks are
                    % defined (i.e. have non-empty pointers):
                    if isempty(Lidx)
                        islow = binctrs >= scale(1) & ...
                            binctrs < mu(MIDidx) - 2*sigma(MIDidx);
                        ismid = binctrs >= mu(MIDidx) - 2*sigma(MIDidx);
                    else
                        islow = binctrs >= scale(1) & ...
                            binctrs < mu(Lidx) + 2*sigma(Lidx);
                        ismid = binctrs >= mu(Lidx) + 2*sigma(Lidx);
                    end
                    if isempty(Ridx)
                        ishigh = binctrs <= scale(2) & ...
                            binctrs >= mu(MIDidx) + 2*sigma(MIDidx);
                        ismid = ismid & ...
                            binctrs < mu(MIDidx) + 2*sigma(MIDidx);
                    else
                        ishigh = binctrs <= scale(2) & ...
                            binctrs >= mu(Ridx) - 2*sigma(Ridx);
                        ismid = ismid & ...
                            binctrs < mu(Ridx) - 2*sigma(Ridx);
                    end
                otherwise
                    error('lfp_markBadCSC:badfix1', ...
                        'This is logically impossible.');
            end
            % islow should now encompass the left peak, if there is one,
            % and ishigh the right peak, if there is one, and ismid the
            % middle peak.  If any of these ranges is empty, that means the
            % fitting didn't work.
            if sum(islow)==0 || sum(ismid)==0 || sum(ishigh)==0
                % retrench to using the quartiles to set the clipping
                % scale. In the normal distribution, the second quartile
                % extends from -0.67 sigma to zero, which means 2 sigmas
                % would be 3 times the width of the second quartile; so we
                % replace "2 * sigma" with the quartile formula.  Note that
                % prcts(2) is the median sample value.
                qrt2width = prcts(2) - prcts(1);
                islow = binctrs < prcts(2) - 3 * qrt2width;
                ismid = binctrs >= prcts(2) - 3 * qrt2width & ...
                    binctrs < prcts(1) + 3 * qrt2width;
                ishigh = binctrs >= prcts(1) + 3 * qrt2width;
                % <islow> and <ishigh> are not guaranteed to include any
                % bins, so if they include nothing, then they should be
                % fudged to at least include the extreme bins.
                if ~any(islow)
                    islow(1) = true;
                end
                if ~any(ishigh)
                    ishigh(end) = true;
                end
            end
            if sum(islow)==0 || sum(ismid)==0 || sum(ishigh)==0
                error('lfp_markBadCSC:badfix2', ...
                    'This is logically impossible.');
            end
            if max(nmasked(ismid)) > 0
                % clipping is defined as the ratio of maximum count in
                % the clipping peaks to maximum count in the middle.
                clipping(fidx) = max(nmasked(islow | ishigh)) / ...
                    max(nmasked(ismid));
            else
                % For some reason, Inf is bad here so we assign
                % 'realmax':
                clipping(fidx) = realmax;
            end
            if isempty(MIDidx)
                signalSD(fidx) = std(lfp_Samples{fn}( ...
                    lfp_Samples{fn}(:) >= binctrs(find(ismid, 1)) ...
                    & lfp_Samples{fn}(:) < binctrs(find(ismid, 1, 'last')) ));
            else
                % Consider the
                signalSD(fidx) = sigma(MIDidx);
            end
        end
    end
    if flags(2)
        % compute full-band smoothed spectrum:
        lfp_FreqLim = [];
        BL = lfp_BLspectrum([], fn, ...
            {lfp_NominalTrialStart lfp_NominalTrialEnd}, [0 0], ...
            winwidth, 'rmdc');
        smoothBL = BL;
        smoothBL.sum = conv(BL.sum, han5, 'same') / sum(han5);
        % compute <choppiness>:
        isinband = BL.f >= freqlim(1) & BL.f < freqlim(2);
        choppiness(fidx) = sum(abs(diff(log10(smoothBL.sum(isinband)))));
        % find and count bad peaks:
        [freq, pwr, width, indices] = lfp_findSpectralPeaks(BL.f, ...
            smoothBL.sum, 10);
        normpwr = NaN(size(pwr));
        sub30pwr = sum(smoothBL.sum(smoothBL.f>freqlim(1) & smoothBL.f<30));
        validwidth = ~isnan(width);
        normwidth = width ./ freq;
        isbadpeak = normwidth < 0.1;
        peakinband = freq >= freqlim(1) & freq < freqlim(2);
        peakidx = find(validwidth & peakinband & isbadpeak);
        for peakidx2 = 1:length(peakidx)
            normpwr(peakidx(peakidx2)) = sum(smoothBL.sum( ...
                smoothBL.f>freq(peakidx(peakidx2))-width(peakidx(peakidx2)) ...
                & smoothBL.f<freq(peakidx(peakidx2))+width(peakidx(peakidx2)) ...
                ))/sub30pwr;
        end
        numbadpeaks(fidx) = ...
            sum(isbadpeak & validwidth & normpwr>0.8 & peakinband);
        % interpolate away peaks and fit pink spectrum:
        peakidx = find(validwidth & freq >= freqlim(1));
        for peakidx2 = 1:length(peakidx)
            idx1 = indices(peakidx(peakidx2), 1);
            idx2 = indices(peakidx(peakidx2), 2);
            smoothBL.sum(idx1 : idx2) = linspace( smoothBL.sum(idx1), ...
                smoothBL.sum(idx2), idx2 - idx1 + 1 );
        end
        [~, fmaxidx] = max(BL.sum);
        fmax(fidx) = BL.f(fmaxidx);
        pinkBL = dg_fitPinkSpec(BL, 'logfit', ...
            'freqlim', [40 freqlim(2)]);
        diffspec = 10 * (log10(BL.sum/BL.N) - log10(pinkBL.sum));
        LFpwr(fidx) = max(diffspec(1:find(BL.f < freqlim(1), 1, 'last')));
        if isempty(BLs)
            BLs = BL;
        else
            BLs(fidx) = BL; %#ok<AGROW>
        end
    end
end
% Make <clipping> tractable for 'classify' function:
clipping(isnan(clipping)) = max(clipping(:));

[isbadfile, isdeadfile] = lfp_markBadCSCcriteria(medians, clipping, ...
    fmax, numbadpeaks, choppiness, LFpwr, binwidths, signalSD, ...
    bincounts, winwidth);

if isempty(destdir)
    matfilename = fullfile(sessiondir, 'lfp_markBadCSC.mat');
else
    [~, sessionID] = fileparts(sessiondir);
    matfilename = fullfile(destdir, ...
        sprintf('lfp_markBadCSC_%s.mat', sessionID) );
end
save(matfilename, 'filenames', ...
    'isbadfile', 'isdeadfile', 'medians', 'binwidths', 'bincounts', ...
    'clipping', 'signalSD', 'fmax', 'numbadpeaks', 'choppiness', ...
    'LFpwr', 'BLs', 'freqlim', 'winwidth');
fprintf('Saved %s\n', matfilename);

