function result = lfp_findPosThresh2(filenum, thresh, option)
%lfp_findPosThresh for lfp_createEvents finds positive-going threshold
% crossings
%result = lfp_findPosThresh(filenum, thresh)
%result = lfp_findPosThresh(filenum, thresh, 'sample')

% Finds the positive-going crossings of <thresh> in the waveform in
% <filenum>.  For each crossing, returns a timestamp midway between the two
% samples straddling the threshold.  If <option> is given and is 'sample',
% then returns sample index numbers of last sample before crossing.
% Uses >, <=, so can trigger on waveform that just touches threshold.

%$Rev: 225 $
%$Date: 2011-05-19 20:01:52 -0400 (Thu, 19 May 2011) $
%$Author: dgibson $

lfp_declareGlobals;

sampleflag = false;
if nargin >= 3
    switch option
        case 'sample'
            sampleflag = true;
    end
end

result = find( ...
    (lfp_Samples{filenum}(2 : end) > thresh) ...
    & (lfp_Samples{filenum}(1 : end - 1) <= thresh) );
if ~sampleflag
    result = lfp_index2time(result) + lfp_SamplePeriod/2;
end