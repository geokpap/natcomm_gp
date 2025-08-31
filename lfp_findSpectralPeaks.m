function [freq, pwr, width, indices] = lfp_findSpectralPeaks(f, P, N)
% Peaks are returned in order of decreasing power.
%INPUTS
% f: frequencies
% P: power, same length as f
% N: number of peaks to find
%OUTPUTS
% freq: a list of <N> center frequencies of peaks.  May contain NaN if
%   fewer than <N> peaks were present.
% pwr: peak maximum power, matching <freq>.  May contain NaN if
%   fewer than <N> peaks were present.
% width: peak width at half maximum power, linearly interpolated between
%   points in <f>.  May contain NaN if fewer than <N> peaks were present,
%   and also in the cases where a neighboring peak made it impossible to
%   follow the measured peak down to half power.
% indices: the indices into f at which each peak starts and ends.
%NOTES
% No attempt is made to interpolate any particular peak shape, so <pwr> and
% consequently <width> could be substantially inaccurate for sufficiently
% sharp peaks.  Neither is any attempt made to ignore peaks that fail to
% descend to at least half power on both sides; these must be eliminated a
% priori by smoothing, or else they will show up in the output list with
% totally misleading width measurements.  Peaks are found in order of
% decreasing power by simply setting P=0 at all points between the
% <indices> belonging to the peak, and repeating the process.

%$Rev: 337 $
%$Date: 2015-02-09 22:51:32 -0500 (Mon, 09 Feb 2015) $
%$Author: dgibson $

freq = NaN(N,1);
pwr = NaN(N,1);
width = NaN(N,1);
indices = NaN(N,2);

for k = 1:N
    [fidx, pwr(k), width(k)] = findOnePeak(f, P);
    freq(k) = f(fidx);
    indices(k,:) = findPeakRange(fidx, P);
    P(indices(k,1):indices(k,2)) = 0;
    if all(P==0)
        return
    end
end
end

function [fidx, pwr, width] = findOnePeak(f, P)
% fidx: index into <P> and <f>.
[pwr, fidx] = max(P);
idx1 = find(P(1:fidx-1) < pwr/2, 1, 'last');
if isempty(idx1)
    idx1 = 1;
end
idx2 = find(P(fidx+1:end) < pwr/2, 1, 'first') + fidx;
if isempty(idx2)
    idx2 = length(P);
end
m = (P(idx1+1) - P(idx1)) / (f(idx1+1) - f(idx1));
f1 = f(idx1) + (pwr/2 - P(idx1)) / m;
m = (P(idx2) - P(idx2-1)) / (f(idx2) - f(idx2-1));
f2 = f(idx2-1) + (pwr/2 - P(idx2-1)) / m;
if P(idx1) == 0 || P(idx2) == 0
    width = NaN;
else
    width = f2 - f1;
end
end

function indices = findPeakRange(fidx, P)
% Find the range of indices into <P> where <P> is monotonic decreasing as
% one moves away from <fidx>.
idx1 = find(P(2:fidx) <= P(1:fidx-1), 1, 'last') + 1;
if isempty(idx1)
    idx1 = 1;
end
idx2 = find(P(fidx+1:end) >= P(fidx:end-1), 1) + fidx;
if isempty(idx2)
    idx2 = length(P);
end
indices = [idx1 idx2];
end



