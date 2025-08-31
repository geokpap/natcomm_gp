function [pks, smoothings, matches] = dg_findPeaksRobust(curve, thresh, tol)
% Performs majority vote peak detection using the raw <curve>, a
% 3-point smoothed version, and a 5-point smoothed version.
%INPUTS
% curve: numeric vector containing the curve to find peaks in.
% thresh: min acceptable difference in value between peak and nearest
%   valley in smoothed curve.
% tol: number of points (indices into <curve>) of mismatch to tolerate
%   between peak locations in different smoothings.  A warning is raised if
%   the mismatch exceeds <tol> by only one point.
%OUTPUT
% pks: table of peaks that pass the majority vote logic, one row per
%   peak, col. 1 = index into the more smoothed version of <curve> where
%   the peak was found, col. 2 = peak value in that smoothing.  The order
%   of the peaks is determined by the vagaries of which smoothing revealed
%   which peaks, so it is effectively random.  If the same index is listed
%   more than once in col. 1, then this means the peak was considered to
%   match multiple peaks in one or both of the other smoothings.  One way
%   to resolve the "true" peak position is to use the matches that were
%   found in the least smoothed version of <curve> (see <matches>), and if
%   there is more than one such match, use the one that is closest to the
%   value in pks(:,1).  However, in some applications it might be
%   preferable simply to use the peak index in the more-smoothed curve,
%   i.e. use unique(pks(:,1)) as the list of peaks' positions.
% smoothings: same number of rows as <pks>, 2 columns containing the
%   smoothing indices (1, 2, or 3) indicating which pair of smoothings the
%   corresponding row of <pks> was found between.
% matches: the index into the less smoothed version of <curve> where the
%   matching peak was found.  Column vector with same number of rows as
%   <pks>.
%NOTES
% It's hard to predict when this will get to be a performance problem, but
% for the sake of clarity, I have used a simple brute force exhaustive
% pairwise comparison between the peak lists for each smoothing.  That
% should be fine as long as there aren't too many peaks in each list (see
% code, <pklist1idx>).

%$Rev: 244 $
%$Date: 2016-07-28 18:16:26 -0400 (Thu, 28 Jul 2016) $
%$Author: dgibson $

peakidx = cell(3, 1); % each cell contains 2-col peak table like <pks>
% 2*<smpts> + 1 points of smoothing:
for smpts = 0 : 2
    if smpts > 0
        smcurve = dg_convpropa(curve, ones(2 * smpts + 1, 1));
    else
        smcurve = curve;
    end
    isincreasing = diff(smcurve) > 0;
    pkidx = find(isincreasing(1:end-1) & ~isincreasing(2:end)) + 1;
    % Remove any peaks that fail the <thresh> dB nearest valley test:
    peaks2rm = [];
    for peakidx2 = 1:size(pkidx, 1)
        idx = pkidx(peakidx2) - 1;
        while idx > 1 && smcurve(idx - 1) < smcurve(idx)
            idx = idx - 1;
        end
        Lvalley = smcurve(idx);
        idx = pkidx(peakidx2) + 1;
        while idx < length(smcurve) && ...
                smcurve(idx + 1) < smcurve(idx)
            idx = idx + 1;
        end
        Rvalley = smcurve(idx);
        if smcurve(pkidx(peakidx2)) < Lvalley + thresh || ...
                smcurve(pkidx(peakidx2)) < Rvalley + thresh
            peaks2rm(end+1) = peakidx2;  %#ok<AGROW>
        end
    end
    pkidx(peaks2rm) = [];
    if ~isempty(pkidx)
        peakidx{smpts+1} = reshape(pkidx, [], 1);
        peakidx{smpts+1}(:, 2) = smcurve(pkidx);
    end
end
% There are three comparisons being made in order to find matching peaks:
% peakidx{1} vs. peakidx{2}, peakidx{1} vs. peakidx{3}, and peakidx{2} vs.
% peakidx{3}. The majority vote rule thus translates into:
%   A peak will be included in the output if it has a match in at least one
%   other list (i.e. it was found in at least two out of three smoothings).
% The index into the less smoothed curve is returned in <pks>.
pks = zeros(0,2);
smoothings = zeros(0,2);
matches = zeros(0,1);
for pklist1idx = 1:2
    % For each element in peak list # <pklist1idx>:
    for pklist1idx2 = 1:size(peakidx{pklist1idx}, 1)
        for pklist2idx = pklist1idx + 1 : 3
            if isempty(peakidx{pklist2idx})
                continue
            end
            absdiffs = abs( peakidx{pklist2idx}(:, 1) ...
                - peakidx{pklist1idx}(pklist1idx2, 1) );
            % logical indices into peakidx{pklist2idx}:
            ismatch = absdiffs <= tol;
            isclose = absdiffs == tol + 1;
            % The logic for raising warnings becomes hideously complicated
            % if one attempts to avoid redundant warnings (i.e. where a
            % peak has a match in one list and is close but not matching in
            % the other), so I'm not attempting that.
            if any(isclose) && ~any(ismatch)
                warning( 'dg_findPeaksRobust:close', ...
                    'Almost-match between %d and %d point smoothings: %d, %d', ...
                    (pklist1idx - 1) * 2 + 1, (pklist2idx - 1) * 2 + 1, ...
                    peakidx{pklist1idx}(pklist1idx2), ...
                    peakidx{pklist2idx}(isclose) );
            end
            if any(ismatch)
                if sum(ismatch) > 1
                    % Choose the closest one:
                    [~, idx] = min(absdiffs(ismatch));
                    ismatch(:) = false;
                    ismatch(idx) = true;
                end
                % Add new row to outputs:
                pks(end+1, :) = peakidx{pklist2idx}(ismatch, :); %#ok<AGROW>
                smoothings(end+1, :) = [pklist1idx pklist2idx]; %#ok<AGROW>
                matches(end+1, 1) = peakidx{pklist1idx}(pklist1idx2, 1); %#ok<AGROW>
            end
        end
    end
end

end
