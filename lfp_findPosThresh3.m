function result = lfp_findPosThresh3(filenum, thresh, option)
%lfp_findPosThresh for lfp_createEvents finds positive-going threshold
% crossings
%result = lfp_findPosThresh(filenum, thresh)
%result = lfp_findPosThresh(filenum, thresh, 'sample')

% Finds the positive-going crossings of <thresh> in the waveform in
% <filenum>.  For each crossing, returns a timestamp midway between the two
% samples straddling the threshold.  If <option> is given and is 'sample',
% then returns sample index numbers instead.  Uses >, <=, but explicitly
% removes any false triggers caused by local minima that just touch thresh.
% <result> is a row vector.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

sampleflag = false;
if nargin >= 3
    switch option
        case 'sample'
            sampleflag = true;
        otherwise
            error('lfp_findPosThresh3:badopt', ...
                'No such option: %s', dg_thing2str(option));
    end
end

result = find( ...
    (lfp_Samples{filenum}(2 : end) > thresh) ...
    & (lfp_Samples{filenum}(1 : end - 1) <= thresh) );

% Remove false triggers
badtriggers = [];
borderlineidx = find(lfp_Samples{filenum}(result) == thresh);
for k = 1:length(borderlineidx)
    if result(borderlineidx(k)) == 1
        badtriggers(end+1) = borderlineidx(k);
    end
    idx = result(borderlineidx(k)) - 1;
    while (idx > 0)
        if lfp_Samples{filenum}(idx) > thresh
            badtriggers(end+1) = borderlineidx(k);
            break
        elseif lfp_Samples{filenum}(idx) < thresh
            % good trigger
            break
        end
        idx = idx - 1;
    end
    if idx == 0
        badtriggers(end+1) = borderlineidx(k);
    end
end
result(badtriggers) = [];

if ~sampleflag
    result = lfp_index2time(result) + lfp_SamplePeriod/2;
end

