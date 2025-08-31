function result = lfp_findNegThresh2(filenum, thresh, option)
%lfp_findNegThresh for lfp_createEvents finds negative-going threshold
% crossings
%result = lfp_findNegThresh(filenum, thresh)
%result = lfp_findNegThresh(filenum, thresh, 'sample')

% Finds the negative-going crossings of <thresh> in the waveform in
% <filenum>.  For each crossing, returns a timestamp midway between the two
% samples straddling the threshold.  If <option> is given and is 'sample',
% then returns sample index numbers of last sample before crossing.
% Uses <, >=, so can trigger on waveform that just touches threshold.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
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
    (lfp_Samples{filenum}(2 : end) < thresh) ...
    & (lfp_Samples{filenum}(1 : end - 1) >= thresh) );
if ~sampleflag
    result = lfp_index2time(result) + lfp_SamplePeriod/2;
end