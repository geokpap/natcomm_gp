function [result, units] = lfp_deglitch(filenum, threshold)
%LFP_DEGLITCH is a function intended for use with lfp_createWave.
%result = lfp_deglitch(filenum, threshold)

%result = lfp_deglitch(filenum, threshold)
%  Function intended for use with lfp_createWave.
%  Calls dg_deglitch with threshold on filenum. 

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~isequal(size(filenum), [1 1])
    error('lfp_deglitch:badfilenums', ...
        '<filenum> must be a single number' );
end
if ~isequal(size(threshold), [1 1])
    error('lfp_deglitch:badthreshold', ...
        '<threshold> must be a single number' );
elseif threshold <= 0
    error('lfp_deglitch:badthreshold3', ...
        '<threshold> must be greater than zero' );
end

units = lfp_SamplesUnits{filenum};
[result, numglitch] = dg_deglitch(lfp_Samples{filenum}, threshold);
disp(sprintf(...
    'Removed %d glitches from %d samples (%d%% of samples)', ...
    numglitch, numel(lfp_Samples{filenum}), ...
    round(100*numglitch/numel(lfp_Samples{filenum})) ));
