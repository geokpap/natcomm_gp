function timestamps = lfp_markClipping(filenum, clipvals)
%timestamps = lfp_markClipping(filenum, clipvals)
% Function intended for use with lfp_createEvents.  Returns a timestamp at
% each sample in <filenum> that is exactly equal to one of <clipvals>.
%DEFAULTS
% clipvals = [-2048 2047]
%
%WARNING:
% The default clipvals will not work if voltage conversion is enabled.
% See documentation for lfp_SamplesUnits.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 2
    clipvals = [-2048 2047];
end

if ~isequal(size(filenum), [1 1])
    error('lfp_markClipping:badfilenums', ...
        '<filenum> must have exactly 1 element' );
end

timestamps = lfp_index2time( ...
    find(ismember( ...
    reshape(lfp_Samples{filenum}, [], 1), clipvals )));

