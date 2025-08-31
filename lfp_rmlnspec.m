function [result, units] = lfp_rmlnspec(filenum, F, varargin)
%lfp_rmlnspec(filenum, F)  For use with lfp_createWave.
% Removes frequency components that are within 0.2 Hz of any harmonic of
% <F>.  This is done recording segment by recording segment so as not to
% introduce any new discontinuities.
%OPTIONS
% 'auto' - finds the actual frequency of the nominal <F>-cycle noise
%   fundamental, estimates its width, and uses the results to set values of
%   <F> and <ftol>.  The specified value of <F> must be accurate to within
%   1 part in 50 (which is the half-width of the band searched, unless
%   'searchwidth' is invoked).  If no peak at all can be found in the
%   search band, then error 'lfp_rmlnspec:nopeak' is raised. If the peak is
%   not big enough to satisfy two statistical tests, then warning
%   'lfp_rmlnspec:notsignif' or 'lfp_rmlnspec:notsignif2' is raised.
% 'ftol', ftol - removes components that are within <ftol> of <F>.  Default
%   = 0.2 Hz.  <ftol> is in Hz.  It is assumed that <ftol> is less than
%   F/2.
% 'searchwidth', searchwidth - this only has an effect when 'auto' is also
%   invoked.  <searchwidth> specifes the half-width of the band to be
%   searched as a fraction of <F>. The search band is centered on <F>.
%   Default value is 0.02.
% 'trials' - does the removal trial by trial instead of rec segment by rec
%   segment; this is less demanding of memory, but may introduce
%   artifactual discontinuities at trial boundaries.  Also note that when
%   trials are on the order of 5 seconds or less in duration, the default
%   0.2 Hz value of <ftol> may not be sufficient to eliminate an entire
%   line in the spectrum even if it works fine at rec segment level. Also,
%   the 0.02 default value of <searchwidth> may not be enough to find the
%   intended peak.
%NOTES
% This works better than the other implementations of line noise filtering
% (i.e. lfp_rmlinesc and dg_ch_rmlinesc).  It's hard to imagine what one
% would want those versions for; if a recording segment is so short that
% frequency resolution becomes a problem, then it is surely also so short
% that discontinuities at the splice points would create even worse
% problems if several recording segments were concatenated.
% -DG 11-Mar-2016.

%$Rev: 385 $
%$Date: 2016-09-13 18:44:57 -0400 (Tue, 13 Sep 2016) $
%$Author: dgibson $

global lfp_SamplesUnits lfp_Samples lfp_TrialIndex lfp_RecSegments ...
    lfp_SamplePeriod lfp_ActiveFilenums

autoflag = false;
ftol = 0.2;
searchwidth = 0.02;
trialsflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'auto'
            autoflag = true;
        case 'ftol'
            argnum = argnum + 1;
            ftol = varargin{argnum};
        case 'searchwidth'
            argnum = argnum + 1;
            searchwidth = varargin{argnum};
        case 'trials'
            trialsflag = true;
        otherwise
            error('lfp_rmlnspec:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

if ~ismember(filenum, lfp_ActiveFilenums)
    error('lfp_rmlnspec:filenum', ...
        'No such filenum: %d', filenum )
end;

units = lfp_SamplesUnits{filenum};
result = NaN(size(lfp_Samples{filenum}));

Fs = 1/lfp_SamplePeriod;
N = floor(Fs/F/2); % number of harmonics before Nyquist

if trialsflag
    numsegments = size(lfp_TrialIndex, 1);
    segstart = lfp_TrialIndex(:, 3);
    segend = lfp_TrialIndex(:, 4);
else
    numsegments = size(lfp_RecSegments, 1);
    segstart = lfp_RecSegments(:, 1);
    segend = lfp_RecSegments(:, 2);
end

for segnum = 1:numsegments
    L = segend(segnum) - segstart(segnum) + 1;
    NFFT = 2^nextpow2(L);
    Y = fft(lfp_Samples{filenum}(segstart(segnum):segend(segnum)), NFFT);
    freqperpoint = Fs / NFFT;
    if autoflag
        % Find noise fundamental and width.
        freqs = (0:NFFT-1)*freqperpoint;
        numfreqpts = round(searchwidth * F / freqperpoint);
        [~, idx] = min(abs(freqs - F));
        P = Y((-numfreqpts:numfreqpts) + idx) ...
            .* conj(Y((-numfreqpts:numfreqpts) + idx));
        [~, idx2] = max(P);
        if idx2 == 1 || idx2 == length(P)
            error( 'lfp_rmlnspec:nopeak', ...
                'There is no peak in the %d - %d Hz range.', ...
                freqs(idx - numfreqpts), freqs(idx + numfreqpts) );
        end
        % We look for p < .05 Bonferroni-corrected for the actual number of
        % points we are looking at.  This is actually a z-test, but I don't
        % like the Matlab ztest help page's insistence on <x> being a
        % vector, so although the answer comes out the same, I prefer to
        % spell it out in terms of normcdf.  I use the symmetry about <mu>
        % of normcdf to get around the roundoff problems associated with
        % dealing with values greater than the mean.  Note that this test
        % will produce a lot of false positives (i.e. it will admit a lot
        % of peaks that are truly nothing special) because the raw power is 
        % really not normally distributed at all.  I try to ameliorate that
        % by running the stats on amplitude instead of power, but amplitude
        % isn't really normally distributed either.
        plev = 0.05 / (2 * numfreqpts + 1);
        A = sqrt(P);
        mu = mean(A);
        sigma = std(A);
        peakdiff = A(idx2) - mu;
        p = normcdf(mu-peakdiff, mu, sigma);
        if p > plev
            warning( 'lfp_rmlnspec:notsignif', ...
                'segnum=%d, Raw mean amplitude=%d, std=%d, peak=%d', ...
                mu, sigma, A(idx2) );
        end
        F = freqs(idx - numfreqpts + idx2 - 1);
        % Smooth with kernel 100 times narrower than the search band to
        % find the peak's width, which is defined as the bandwidth where
        % smoothed amplitude is halfway between the median and the peak
        % smoothed amplitude. Assuming the kernel width is considerably
        % more than one point, the central limit theorem predicts that the
        % smoothed values will be much closer to a normal distribution, so
        % the statistical test should work much better here.
        smoothpts = max(1, round(F/(100*50*freqperpoint)));
        Asmooth = conv(A, hanning(smoothpts), 'same');
        [peakamp, pkAidx] = max(Asmooth);
        mu = mean(Asmooth);
        sigma = std(Asmooth);
        peakdiff = peakamp - mu;
        p = normcdf(mu-peakdiff, mu, sigma);
        if p > plev
            warning( 'lfp_rmlnspec:notsignif2', ...
                'segnum=%d, Smoothed mean amplitude=%d, std=%d, peak=%d', ...
                mu, sigma, peakamp );
        end
        halfway = (peakamp + mu)/2;
        Rside = find(Asmooth(pkAidx:end) < halfway, 1) + pkAidx - 1;
        Lside = find(Asmooth(1:pkAidx) < halfway, 1, 'last');
        ftol = ( Rside - Lside + 1 ) * freqperpoint;
    end
    % All the index calculations are "+1" because the index of 0 Hz is 1.
    % First we calculate all the bad frequency points that are below
    % Nyquist. We pre-allocate enough memory to hold all the points in the
    % sidebands, which is roughly the Nth triangular number times the
    % number of points that we want to remove around the fundamental, and
    % then we add increments later on as needed.  The very last harmonic is
    % a special case because the region of bad points has to end by the
    % time we hit Nyquist.
    endharmonicidx = round((N * F) / freqperpoint) + 1;
    badsidebandpts = floor(ftol/freqperpoint);
    numpts = badsidebandpts * N * (N + 1) / 2;
    badfreqidx = zeros(numpts, 1);
    badfreqidx2 = 1; % index into <badfreqidx>
    for k = 1 : (N - 1)
        % The size of the band to remove must scale up in proportion to the
        % center frequency:
        badsidebandpts = floor(k*ftol/freqperpoint);
        % If we need to expand badfreqidx, we want to do it in big chunks,
        % in this case chunks of size numpts/4.  Unless <N> is very small,
        % it will only happen once:
        if badfreqidx2 + 2*badsidebandpts > length(badfreqidx)
            badfreqidx(end + fix(numpts/4)) = 0;
        end
        badfreqidx(badfreqidx2 : badfreqidx2 + 2*badsidebandpts) = ... 
            (-badsidebandpts:badsidebandpts)' + ...
            round(k*F/freqperpoint) + 1;
        badfreqidx2 = badfreqidx2 + 2*badsidebandpts + 1;
    end
    % Now remove the extra zeros at the end:
    badfreqidx(badfreqidx == 0) = [];
    % At this point, all(Y(2:NFFT/2) == conj(Y(NFFT:-1:NFFT/2+2))) because
    % the input was real.  The middle point of the spectrum, Y(NFFT/2+1),
    % is the Nyquist frequency, where the value can only be interpreted as
    % the known minimum quantity of aliased crap (which might be greater
    % but we can't know how much).
    badsidebandpts = floor(N*ftol/freqperpoint);
    badfreqidx = [ badfreqidx
        ( endharmonicidx - badsidebandpts : min( ...
        endharmonicidx + badsidebandpts, ...
        NFFT/2 ) )' ]; %#ok<AGROW>
    % Now we eliminate the bad frequencies in complex conjugate pairs
    % where 2 matches NFFT, 3 matches NFFT - 1, etc., and so (badfreqidx)
    % matches (NFFT - badfreqidx + 2).
    Y([badfreqidx; NFFT - badfreqidx + 2]) = 0;
    y = ifft(Y);
    result(segstart(segnum):segend(segnum)) = ...
        y(1:(segend(segnum) - segstart(segnum) + 1));
end


