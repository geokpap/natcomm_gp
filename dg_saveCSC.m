function dg_saveCSC(matfilename, dg_Nlx2Mat_Timestamps, ...
    dg_Nlx2Mat_Samples, dg_Nlx2Mat_SamplesUnits)
%dg_saveCSC(matfilename, TS, Samples)
% Saves CSC data in lfp_lib-readable format. Does not check to see if file
% exists first, so it will overwrite any such file.
%INPUTS
% matfilename:  absolute or relative pathname to output file; if you do not
%   include an extension, Matlab will add its default.
% TS:  timestamp in seconds for each frame (every 512th sample).
% Samples:  this will be padded with trailing NaNs if there are not an
%   integral multiple of 512 samples.

%$Rev: 201 $
%$Date: 2014-06-26 13:35:18 -0400 (Thu, 26 Jun 2014) $
%$Author: dgibson $

if nargin < 4 || isempty(dg_Nlx2Mat_SamplesUnits)
    dg_Nlx2Mat_SamplesUnits = 'arbs'; %#ok<NASGU>
end
numframes = ceil(numel(Samples) / 512);
if length(dg_Nlx2Mat_Timestamps) ~= numframes
    error('dg_saveCSC:numframes', ...
        'There must be a timestamp for each 512-sample frame.');
end
if mod(numel(Samples), 512) ~= 0
    dg_Nlx2Mat_Samples(end+1 : (numframes*512)) = NaN; %#ok<NASGU>
end
save(matfilename, 'dg_Nlx2Mat_Timestamps', ...
    'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
