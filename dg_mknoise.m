function wave = dg_mknoise(f, P, N, Fs)
% Creates random noise to fit a specified spectrum.  This is done by
% spline interpolating the specified spectrum to yield <N> points,
% assigning uniformly distributed random phase to each frequency, and then
% inverse Fourier transforming to produce the time series.
%INPUTS
% f: vector of frequency points.  The first frequency must be 0, and the
%   last frequency must be <Fs>/2, but the spacing of points in between is
%   arbitrary and can be irregular.
% P: vector of spectral power values (linear, not log scaled).
% N: number of samples of data to produce.  Must be even.
% Fs: sample rate of output.
%OUTPUTS
% wave: vector of <N> samples in time domain at <Fs>.
%NOTES
%   The first element of the Matlab FFT corresponds to frequency 0.  There
% are as many elements in the FFT as there are in the equivalent time
% series, and frequency points are spaced at an interval of Fs/N, so the
% last element in the FFT corresponds to the last frequency point before
% Fs. However, we need to constrain the time series to be real, and this
% places constraints on the FFT: the values at 0 frequency and at the
% Nyquist frequency (Fs/2) must be real, and the remaining values must be
% conjugate symmetric around the Nyquist frequency.  Setting the amplitude
% at the Nyquist freq to zero has the effect of truncating the harmonic
% series, so the remaining values in the spectrum are still the best
% possible approximation to the given time series, if there were a given
% time series; but in this case we are just trying to generate a synthetic
% time series, and we really don't care if the bandwidth is one frequency
% point less than the maximum possible.  Indeed, a case could be made for
% ensuring that the specified spectrum asymptotes smoothly to zero as it
% approaches Nyquist.
%   It's up to the user to ensure that the spectrum specified by <f> and
% <P> does not become ridiculous when interpolated.  The Matlab function
% 'spline' is used to do the interpolation, and if the points in <f> are
% too coarsely spaced with respect to the size of spectral features
% specified, then 'spline' produces large overshoots.  Some sort of
% automation is clearly needed here, possibly involving transforming a raw
% spectrum to a logarithmic frequency scale, smoothing that and setting
% frequency points with spacing inversely proportional to the second
% derivative, then transforming back. 

%$Rev: 267 $
%$Date: 2019-09-19 16:17:47 -0400 (Thu, 19 Sep 2019) $
%$Author: dgibson $

if f(1) ~= 0 || abs(f(end)/(Fs/2) - 1) > 2 * eps
    error('dg_mknoise:f', ...
        'Please see header comments on values of <f>.');
end
if mod(N, 2)
    error('dg_mknoise:N', ...
        '<N> must be even.');
end

fgrid = Fs * (0 : (N-1))' / N;
Nyqidx = N/2 + 1;
A = spline(f, sqrt(P), fgrid); % amplitude (magnitude) of spectrum
A(Nyqidx) = 0;
A(2:Nyqidx-1) = A(2:Nyqidx-1) .* exp(1i*2*pi*rand(Nyqidx-2, 1));
A(Nyqidx+1:end) = conj(A((Nyqidx-1):-1:2));
% A desirable effect of using dpss tapers in chronux is that they produce
% spectral density values that are independent of the number of samples.
% This is because the vector magnitude of each taper is normalized to 1, so
% as N increases, the peak magnitude of the taper has to decrease by
% sqrt(N).  So applying the taper to the signal not only tapers the signal,
% it also decreases its amplitude by sqrt(N).  Chronux also converts
% spectral densities from (signalunits)^2 / freqpoint to (signalunits)^2 /
% Hz by dividing the signal by sqrt(Fs). Here we do the inverse of both of
% those scalings, for compatibility with chronux-style spectra:
wave = ifft(A) * sqrt(N * Fs);
if any(imag(wave))
    error('dg_mknoise:argh', ...
        'Matlab is broken!  Hardware failure!  Everyone to get from street!');
end

