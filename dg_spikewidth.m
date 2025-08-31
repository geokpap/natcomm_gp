function [poswidth, delay, negwidth] = dg_spikewidth(waveform, varargin)
%width = dg_spikewidth(waveform)
% Returns the width in points of the spike waveform in <waveform>.  It is
% assumed that the absolute maximum of the spike waveform represents the
% peak of the main positive-going component of the spike, and the
% absolute minimum of the samples after the peak represents the depth of
% the valley following the peak.
%INPUTS
% <waveform> is a vector of data points representing the spike waveform.
%   For width calculations, it is linearly interpolated.
%OUTPUTS
% <poswidth> is the half-height width (with respect to zero) of the
%	positive phase.  If the positive phase runs off the end of <waveform>,
%	then <poswidth> is NaN.
% <delay> is the width of the interval between the half-height end of the
%   peak and the half-depth start of the valley.
% <negwidth> is the half-height width (with respect to zero) of the
%	positive phase.  If the negative phase runs off the end of <waveform>,
%	then <negwidth> is NaN.
%OPTIONS
% 'firstsampleref' - uses the value of the first sample instead of zero as
%   the reference for computing the half-height and half-depth points.
% 'allpoints' - each value returned becomes a three-element column vector
%   whose first element is the width as normally returned, and the second
%   and third elements are the start and end points of the relevant time
%   period.

%$Rev: 84 $
%$Date: 2010-09-29 18:59:20 -0400 (Wed, 29 Sep 2010) $
%$Author: dgibson $

refval = 0;
argnum = 1;
allpointsflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'allpoints'
            allpointsflag = true;
        case 'firstsampleref'
            refval = waveform(1);
        otherwise
            error('dg_spikewidth:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

poswidth = NaN;
delay = NaN;
negwidth = NaN;
peakstart = NaN;
peakend = NaN;
valstart = NaN;
valend = NaN;

[peak, peakidx] = max(waveform);
peakheight = peak - refval;
peaksamples = find(waveform > (peakheight/2 + refval));
if any(diff(peaksamples) > 1)
    % There are gaps in the "peak".  Remove the disconnected bits.
    myidx = peakidx;
    while myidx > 1 && ismember(myidx-1, peaksamples)
        myidx = myidx - 1;
    end
    if myidx == 1 && ismember(myidx, peaksamples)
        % the peak extends back to the first sample; do nothing
    else
        peaksamples = peaksamples(peaksamples > myidx - 1);
    end
    myidx = peakidx;
    while myidx < length(waveform) && ismember(myidx+1, peaksamples)
        myidx = myidx + 1;
    end
    if myidx == length(waveform) && ismember(myidx, peaksamples)
        % the peak extends to the last sample; do nothing
    else
        peaksamples = peaksamples(peaksamples < myidx + 1);
    end
end
if peakidx < length(waveform) ...
        && ~isempty(peaksamples) ...
        && peaksamples(end) < length(waveform)
    if peaksamples(1) > 1
        peakstart = interp1(waveform([peaksamples(1) - 1, peaksamples(1)]), ...
            [peaksamples(1) - 1, peaksamples(1)], peakheight/2 + refval );
    else
        peakstart = 1;
        warning('dg_spikewidth:earlypeak', ...
            'The first sample is part of the peak; <poswidth> is thus not accurate.');
    end
    peakend = interp1(waveform([peaksamples(end), peaksamples(end) + 1]), ...
        [peaksamples(end), peaksamples(end) + 1], peakheight/2 + refval );
    poswidth = peakend - peakstart;
    [valley, validx] = min(waveform(peakidx:end));
    validx = validx + peakidx - 1;
    if validx == length(waveform)
        warning('dg_spikewidth:novalley2', ...
                'The valley is at the last point.' );
    else
        valdepth = valley - refval;
        valsamples = find(waveform(peakidx:end) < (valdepth/2 + refval)) ...
            + peakidx - 1;
        if any(diff(valsamples) > 1)
            % There are gaps in the "valley".  Remove the disconnected bits.
            myidx = validx;
            while myidx > 1 && ismember(myidx-1, valsamples)
                myidx = myidx - 1;
            end
            if myidx == 1 && ismember(myidx, valsamples)
                % the valley extends back to the first sample; that can't
                % happen.
                error('dg_spikewidth:bug', ...
                    'Something impossible has happened.  Please panic.' );
            else
                valsamples = valsamples(valsamples > myidx - 1);
            end
            myidx = validx;
            while myidx < length(waveform) && ismember(myidx+1, valsamples)
                myidx = myidx + 1;
            end
            if myidx == length(waveform) && ismember(myidx, valsamples)
                % the valley extends to the last sample; do nothing
            else
                valsamples = valsamples(valsamples < myidx + 1);
            end
        end
        if isempty(valsamples)
            warning('dg_spikewidth:novalley', ...
                'There is no valley.' );
        else
            valstart = interp1(waveform([valsamples(1) - 1, valsamples(1)]), ...
                [valsamples(1) - 1, valsamples(1)], valdepth/2 + refval );
            delay = valstart - peakend;
            if valsamples(end) < length(waveform)
                valend = interp1(waveform([valsamples(end), valsamples(end) + 1]), ...
                    [valsamples(end), valsamples(end) + 1], valdepth/2 + refval );
                negwidth = valend - valstart;
            end
        end
    end
end

if allpointsflag
    poswidth(2:3, 1) = [peakstart; peakend];
    delay(2:3, 1) = [peakend; valstart];
    negwidth(2:3, 1) = [valstart; valend];
end
