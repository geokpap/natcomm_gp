function [bwHz, dur, sn, c, P] = dg_getTaperStats(nw, k, T, plotlims)
%DG_GETTAPERSTATS returns/displays useful statistics for Slepian tapers.
%[bwHz, dur, sn, c, P] = dg_getTaperStats(nw, k)
% Window width T defaults to 1.
%[bwHz, dur, sn, c, P] = dg_getTaperStats(nw, k, T)
% bandwidth: half-power full-bandwidth of the spectral power blur function
% duration: between half-amplitude points of last taper (% of T)
%[bwHz, dur, sn, c, P] = dg_getTaperStats(nw, k, T, plotlims)

% <nw>, <k> are inputs to Matlab function dpss(n, nw, k).
% <T> is the duration in seconds of the window width; default = 1.0.
% <plotlims> is an optional argument which, if present, controls how the
% blur function is plotted.  Plotting is skipped entirely if <plotlims> =
% 0; otherwise the first row of plotlims sets the x limits (remember that
% frequencies can also be negative) and the second row sets y limits (the
% plot is in dB re: max, so the values are all <= 0).
% <bwHz> is the half-power bandwidth of the spectral power blur function.
% <dur> is the duration between the most extreme points where the last
% taper attains an absolute value of at least half its maximum absolute
% value.
% <sn> is the ratio in dB between the top of the spectral power blur
% function and the top of the side lobes, which are here defined to include
% everything to the right of the first minimum beyond the half-power
% frequency.
% <c> is the average concentration of the taper set.
% <P> is the spectral power blur function.
%
% The number of sample points used is always 512, so the sampling rate is
% always T/512 - but this only matters at high frequencies that are way out
% on the tails of the blur function anyway.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin <2
    error('dg_getTaperStats:insufficientArgs', ...
        'Not enough arguments; type "help dg_getTaperStats" for details.' );
end
if nargin < 3
    T = 1;
end
if nargin < 4
    plotlims(1,:) = [-10 10]/T;
    plotlims(2,:) = [-50 0];
end

npts = 512;
[tapers, c] = dpss(npts,nw,k);
padfactor = 8;
S = fft(tapers, npts * padfactor);
P = mean(S.*conj(S),2);
[pmax, pmaxidx] = max(P);
P = P / pmax;

if ~isequal(plotlims, 0)
    % Plot the blur function
    freqs = (0:npts*padfactor/2)/(padfactor*T);
    figure;
    % Nyquist freq is plotted at beginning and again at end
    plot([-freqs(end:-1:2)'; freqs'], ...
        10*log10([P((1+end/2):-1:2); P(1:(1+end/2))]));
    title(sprintf('Blur Function for\nnw = %.2f, k = %.0f, T = %d s', nw, k, T));
    xlabel('Hz');
    ylabel('dB');
    grid on;
    xlim(plotlims(1,:));
    ylim(plotlims(2,:));
end

bwPerSample = fzero(@(x) (interp1(P, x) - 1/2), [1 ceil(length(P)/2)]);
bwHz = 2 * (bwPerSample-1)/(padfactor*T);

abstaper =abs(tapers(:,end));
abstaper = abstaper/max(abstaper);
includedpts = find(abstaper > sqrt(.5));
endptsfrac = includedpts([1 end]) / npts;
dur = T * (endptsfrac(2) - endptsfrac(1));

bwidx = round(bwPerSample);
increasing = find(P(bwidx+1:end) > P(bwidx:end-1)) + bwidx - 1;
sn = 10 * log10(1/max(P(increasing(1):round(length(P)/2))));

c = mean(c);

if nargout == 0
    disp(sprintf(...
        'bandwidth: %1.2g Hz\nduration: %2.2d%% (%1.2g s)\nS/N: %1.1d dB', ...
        bwHz, round(100*(endptsfrac(2) - endptsfrac(1))), dur, round(sn) ));
    disp(sprintf(...
        'concentration: %1.2d%%', round(100*c) ));
end
