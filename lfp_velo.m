function [result, units] = lfp_velo(filenums, smoothing, varargin)
%lfp_velo:  Function intended for use with lfp_createWave.
%result = lfp_velo(filenums, smoothing)
%result = lfp_velo(..., 'gain')
%result = lfp_velo(..., 'raw')
%result = lfp_velo(..., 'seconds')
%result = lfp_velo(..., 'order', N)

%result = lfp_velo(filenums)
%  Smooths the waves in lfp_Samples{filenums(1:2)} using a Hanning window
%  2*<smoothing>+1 bins wide, and calculates the magnitude of the first order
%  difference vector, i.e. the magnitude of the velocity when filenums(1:2)
%  represent position.
%result = lfp_velo(..., 'gain')
%  Compensates for the gain of the Hanning window so that gain is unity for
%  constant (DC) input.  This is default.
%result = lfp_velo(..., 'raw')
%  No gain compensation.
%result = lfp_velo(..., 'seconds')
%  Result is in <units> per second instead of per sample.
%result = lfp_velo(..., 'order', N)
%  Instead of the first-order diff, does the nth-order diff, i.e.:
%       diff(x, N)

%$Rev: 422 $
%$Date: 2023-09-08 14:25:38 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

global lfp_SamplesUnits lfp_Samples lfp_SamplePeriod

if ~isequal(lfp_SamplesUnits{filenums(1)}, lfp_SamplesUnits{filenums(2)})
    warning('lfp_velo:units', ...
        'Files are in different units' );
end

if ~isequal(size(filenums), [1 2])
    error('lfp_velo:badfilenums', ...
        '<filenums> must have exactly 2 elements on one row' );
end

argnum = 1;
gainflag = true;
N = 1;
secondsflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'gain'
            % Default; do nothing
        case 'order'
            argnum = argnum + 1;
            N = varargin{argnum};
        case 'raw'
            gainflag = false;
        case 'seconds'
            secondsflag = true;
    end
    argnum = argnum + 1;
end

if secondsflag
    units = sprintf('%s/sec', lfp_SamplesUnits{filenums(1)});
else
    units = sprintf('%s/samp', lfp_SamplesUnits{filenums(1)});
end

if smoothing == 0
    result = sqrt((lfp_Samples{filenums(1)}(2:end)-lfp_Samples{filenums(1)}(1:end-1)).^2 + ...
        (lfp_Samples{filenums(2)}(2:end)-lfp_Samples{filenums(2)}(1:end-1)).^2);
else
    hw = hanning(2*smoothing+1);
    % s1 and s2 have an "extra" <smoothing> points on each end, so the first
    % "valid" point of each is (smoothing+1).
    s1 = conv(reshape(lfp_Samples{filenums(1)}, 1 ,[]), hw, 'same');
    s2 = conv(reshape(lfp_Samples{filenums(2)}, 1 ,[]), hw, 'same');
    result = sqrt(diff(s1, N).^2 + diff(s2, N).^2);
end
clear s1 s2;
% But result has <N> fewer points than s1 and s2, so pad it with
% duplicates of the beginning and end points:
ptsbefore = ceil(N/2);
ptsafter = floor(N/2);
result = [ repmat(result(1), 1, ptsbefore) ...
    reshape(result, 1, []) ...
    repmat(result(end), 1, ptsafter) ];
if gainflag && (smoothing > 0)
    result = result / sum(hw);
end
if secondsflag
    result = result / lfp_SamplePeriod;
end
