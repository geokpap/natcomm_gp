function [calibrated, units] = lfp_calibEyeTrace(filenum, axis)
%calibrated = lfp_calibEyeTrace(filenum, axis)
%   Returns a column vector containing a calibrated eye trace using
%   filenum as raw data.  If the output value is not
%   used, then the calibrated waveform is constructed directly in the
%   filenum, replacing the raw values, and nothing is returned.
%   <axis> is a string that may have any of the values 'x', 'X', 'y', 'Y'.
%   Eye trace data recorded before the first calib use the first calib;
%   other data use the last calib before the recording time.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals

Nsamples = numel(lfp_Samples{filenum});
units = 'pix';

inplaceflag = false;
if nargout == 0
    inplaceflag = true;
else
    calibrated = zeros(Nsamples, 1);
end

% sort the calibrations in temporal order:
calibtable = lfp_getEyeCalibTable;
calibtimes = cell2mat(calibtable(:,1));
sortidx = sortrows([calibtimes (1:length(calibtimes))']);
sortedtimes = calibtimes(sortidx(:,2));
sortedcalibs = calibtable(sortidx(:,2), 2);

% delete any calibs with bad timestamps, except for the very first one
maxTime = lfp_index2time(Nsamples);
bad = find(sortedtimes > maxTime);
bad(find(bad==1)) = [];
sortedtimes(bad) = [];
sortedcalibs(bad) = [];

endsample = 0;
for idx = 1 : length(sortedtimes)
    startsample = endsample + 1;
    if idx < length(sortedtimes)
        endsample = lfp_time2index(sortedtimes(idx+1)) - 1;
    else
        endsample = Nsamples;
    end
    if inplaceflag
        lfp_Samples{filenum}(startsample : endsample) = ...
            lfp_calibEye(lfp_Samples{filenum}(startsample:endsample), ...
            axis, sortedcalibs{idx} );
        lfp_SamplesUnits{filenum} = units;
    else
        calibrated(startsample:endsample) = ...
            lfp_calibEye(lfp_Samples{filenum}(startsample:endsample), ...
            axis, sortedcalibs{idx} );
    end
end