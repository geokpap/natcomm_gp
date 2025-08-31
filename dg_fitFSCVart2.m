function startidx = dg_fitFSCVart2(ts, samples, thresh)
%INPUTS
% ts: timestamps in seconds of all samples, same size as <samples>.
% samples: array of sampled values, one per timestamp.
% thresh: positive-going threshold to trigger interpolation.
%OUTPUTS
% startidx: the index into <samples> of the start of each identified FSCV
%   artifact.
%NOTES
%   Inspired by the gross deficiencies of 'dg_fitFSCVart' in real life,
% this function uses a different strategy.  Positive-going threshold
% crossings are used as a rough first approximation to the FSCV scan times,
% so the best place for the threshold to be is right below the clipping
% level of the FSCV artifacts.  If the slope of the two samples at and
% immediately before the threshold crossing is less than 2% of <thresh> per
% sample, then the threshold crossing is assumed not to actually be an FSCV
% artifact and is removed from the list.  (Note that this criterion also
% removes crossings where the baseline is so close to threshold that there
% is less than 300 microseconds between the start of the scan and the
% threshold crossing.) For the remaining threshold crossings, the scanperiod 5
% ms before the threshold crossing is fitted with a two-part piecewise
% linear function, and the intersection between the two lines is taken to
% be the start of the scan.
%   A simple linear interpolator is provided in dg_rmFSCVart.
%   The sample scanperiod is estimated based only on the first 10,000 samples.

%$Rev: 269 $
%$Date: 2019-10-07 17:40:04 -0400 (Mon, 07 Oct 2019) $
%$Author: dgibson $

sampleperiod = median(diff(ts(1:1e4)));
% These constants might have to become arguments:
fitwindow = floor(3e-3 / sampleperiod); % for piecewise linear fit
scanperiod = 100e-3 / sampleperiod; % FSCV scanning period
minscansamp = round(.001/sampleperiod); % min number suprathreshold samples
periodicitytol = 0.1; % max fractional deviation from the excpected period

% Find threshold crossings ("xings"), scanperiod, and start of first scanperiod.
% <xingidx> points to the first sample that is above <thresh> in each
% putative FSCV artifact.
isxing = [ false
    reshape( samples(2:end) >= thresh & samples(1:end-1) < thresh, ...
    [], 1 ) ];
xingidx = reshape(find(isxing), [], 1 );
startidx = zeros(size(xingidx));
lastgoodxingidx = [];
for idx2 = 1:length(xingidx)
    slope = samples(xingidx(idx2)) - samples(xingidx(idx2) - 1);
    if slope < 0.02 * thresh
        % Bogus xing, skip:
        continue
    end
    % Check to make sure there really is something that could possibly be
    % an FSCV artifact:
    if sum(samples(xingidx(idx2) + (0:(2*fitwindow))) > thresh) ...
            < minscansamp
        % Ignore this one, it's too brief.
        continue
    end
    % Check periodicity:
    if ~isempty(lastgoodxingidx)
        interval = xingidx(idx2) - lastgoodxingidx;
        remainder = rem(interval, scanperiod);
        if (interval < (1 - periodicitytol) * scanperiod) || ...
                ( remainder > periodicitytol * scanperiod && ...
                remainder < (1 - periodicitytol) * scanperiod )
            % More than 10% off periodicity
            continue
        end
    end
    % Extract the segment of waveform to analyze.
    if xingidx(idx2) - fitwindow < 1
        startsamp = 1;
    else
        startsamp = xingidx(idx2) - fitwindow;
    end
    % <samp2fit> is a column vector of sample values to fit:
    samp2fit = reshape(samples(startsamp:xingidx(idx2)), [], 1);
    % <isartifact> is true for any point preceded by a high slope:
    isartifact = [false; diff(samp2fit) >= 0.02 * thresh];
    % Check to make sure we're not looking at oscillatory junk:
    numtrailingedge = sum(isartifact(1:end-1) & ~isartifact(2:end));
    if numtrailingedge > 10
        % This requires manual inspection!
        warning('dg_fitFSCVart2:ambiguous', ...
            'Numtrailingedge=%d for xingidx at sample %d, skipping', ...
            numtrailingedge, xingidx(idx2));
        continue
    end
    % Note that <isartifact(end)> is already known to be <true>:
    artstart = find(~isartifact, 1, 'last') + 1;
    artend = length(samp2fit);
    % Iteratively move <artstart> backwards so that any little holes of 3
    % or fewer <false> in a row get included in the "artifact":
    while artstart > 5
        % if there are three or fewer <false> values and at least one
        % <true> value in the 4 points preceding artstart, we can move
        % <artstart> back to the first <true> value:
        if sum(isartifact(artstart-4:artstart-1)) > 0 ...
                && sum(~isartifact(artstart-4:artstart-1)) < 3
            artstart = find(isartifact(artstart-4:artstart-1), 1) ...
                + artstart - 4;
        else
            break
        end
    end
    % At this point, <artstart> either points at the first point preceded
    % by a high slope, or its value is less than 5.  4 is close enough to
    % the beginning of <samp2fit> not to keep worrying about it.
    if artstart > 1
        % include the sample that is the beginning of the first
        % <isartifact> interval:
        artstart = artstart - 1;
    end
    % At this point, <artstart> points at the first sample of "artifact".
    % We deliberately leave one sample unassigned in between the end of the
    % clean LFP <lfpend> and the start of the artifact, since that's
    % exactly the murky region we are trying to define better.
    lfpend = artstart - 2;
    if lfpend < 1
        % The artifact extends all the way back to at least the third
        % sample of <samp2fit>, so there is no hope of fitting anything,
        % and we just treat <startsamp> as the start of the artifact:
        startidx(idx2) = startsamp;
        continue
    end
    % compute slopes and offsets for both linear fits.
    % <B(1,:)> is the y-intercepts, <B(2,:)> is the slopes:
    B = [ones(lfpend, 1) (1:lfpend)'] \ samp2fit(1:lfpend);
    B(:, 2) = [ones(artend - artstart + 1, 1) (artstart:artend)'] ...
        \ samp2fit(artstart:artend);
    % However, the curve always bows upward due to the capacitive
    % component, so we'll get a better extrapolation to the start of the
    % artifact if we line up the end of the fitted segment:
    vshift = samp2fit(artend) - (B(1,2) + B(2,2) * artend);
    % Find the intersection point:
    xisect = round((B(1,1) - (B(1,2) + vshift)) / (B(2,2) - B(2,1)));
    startidx(idx2) = xisect + startsamp - 1;
    lastgoodxingidx = xingidx(idx2);
end
% Delete the entries for the bogus xings:
startidx(startidx == 0) = [];



