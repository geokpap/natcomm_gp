function lfp_interpolateWave_result = ...
    lfp_interpolateWave(waveTS, wavedata)
% Fill each CSC frame with interpolated data; ts is timestamp
% of the frame we are filling, and trkidx (tracker index)
% points at the first of the two tracker time points between
% which we are interpolating.
% waveTS:  timestamps in sec of wavedata
% wavedata:  values at stated timestamps as column vector

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
global lfp_interpolateWave_result;
lfp_interpolateWave_result = zeros(1, ...
    length(lfp_TimeStamps) * lfp_SamplesPerFrame );
hWaitBar = waitbar(0, '', 'Name', 'interpolating data to match CSC sampled format');
trkidx = dg_binsearch(waveTS, lfp_TimeStamps(1));
framenum = 1;
if trkidx > 1
    trkidx = trkidx - 1;
else
    % Initial frame(s) must be filled with padding up to the
    % time point where tracker data start
    ts = lfp_TimeStamps(framenum);
    while waveTS(trkidx) > ts
        % Fill beginning part with padding
        timepoints = ts + (0:511)*lfp_SamplePeriod;
        offset = (framenum-1)*lfp_SamplesPerFrame;
        points2fill = find(timepoints <= waveTS(trkidx));
        lfp_interpolateWave_result(offset+points2fill) = wavedata(1);
        % Interpolate ending part
        while (waveTS(trkidx) <= timepoints(end)) ...
                && (trkidx < length(waveTS))
            points2intrp = find(timepoints > waveTS(trkidx) ...
                & timepoints <= waveTS(trkidx+1) );
            intrpTrkData(offset, timepoints, ...
                points2intrp, 1, wavedata, trkidx, waveTS);
            trkidx = trkidx + 1;
        end
        if trkidx >= length(waveTS)
            % Ran out of wavedata, must pad rest of frames
            padEndTrkData(waveTS, wavedata, 1);
            framenum = length(lfp_TimeStamps);
        else
            framenum = framenum + 1;
            trkidx = max(1, trkidx - 1);    % Backtrack to start next frame
            ts = lfp_TimeStamps(framenum);
        end
    end
end
while framenum < length(lfp_TimeStamps)
    waitbar(framenum/length(lfp_TimeStamps), hWaitBar );
    ts = lfp_TimeStamps(framenum);
    timepoints = ts + (0:511)*lfp_SamplePeriod;
    offset = (framenum-1)*lfp_SamplesPerFrame;
    % interpolate frame # framenum for X and Y channels
    while (waveTS(trkidx) <= timepoints(end)) ...
            && (trkidx < length(waveTS))
        points2intrp = find(timepoints > waveTS(trkidx) ...
            & timepoints <= waveTS(trkidx+1) );
        intrpTrkData(offset, timepoints, ...
            points2intrp, 1, wavedata, trkidx, waveTS);
        trkidx = trkidx + 1;
    end
    if trkidx >= length(waveTS)
        % Ran out of wavedata, must pad rest of frame(s)
        padEndTrkData(waveTS, wavedata, 1);
        framenum = length(lfp_TimeStamps);
    else
        framenum = framenum + 1;
        trkidx = trkidx - 1;    % Backtrack to start next frame
    end
end
close(hWaitBar);


function intrpTrkData(offset, tpoints, points2intrp, ...
    trkcolnum, trkdata, trkidx, trkTS)
% points2intrp can be empty if there is a gap between frames, in which case we
% do nothing.
if isempty(points2intrp) return
end
global lfp_Samples;
global lfp_interpolateWave_result;
m = (trkdata(trkidx+1,trkcolnum)-trkdata(trkidx,trkcolnum)) ...
    / (trkTS(trkidx+1)-trkTS(trkidx)) ;
delta = ( tpoints(points2intrp) ...
    - trkTS(trkidx) ) * m;
lfp_interpolateWave_result(offset + points2intrp) ...
    = trkdata(trkidx,trkcolnum) + delta;


function padEndTrkData(waveTS, wavedata, tkcolnum)
lfp_declareGlobals;
global lfp_interpolateWave_result;
lastTrackerPoint = lfp_time2index(waveTS(end));
lastsample = length(lfp_TimeStamps) * lfp_SamplesPerFrame;
lfp_interpolateWave_result(lastTrackerPoint:lastsample) ...
    = wavedata(end, tkcolnum);
