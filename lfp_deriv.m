function result = lfp_deriv(filenum)
%result = lfp_deriv(filenum)

%result = lfp_deriv(filenum)
%  Function intended for use with lfp_createWave.
%  Computes the first-order difference of the wave in filenum.  Since this
%  shortens the wave by one point, the first point is duplicated at the
%  beginning of the result wave.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~isequal(size(filenum), [1 1])
    error('lfp_deriv:badfilenums', ...
        '<filenum> must have exactly 1 element' );
end
result(2:numel(lfp_Samples{filenum})) = ...
    lfp_Samples{filenum}(2:end) - lfp_Samples{filenum}(1:end-1);
result(1) = result(2);
