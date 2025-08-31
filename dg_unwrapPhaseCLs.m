function unwrapped = dg_unwrapPhaseCLs(phi, varargin)
%unwrapped = dg_unwrapPhaseCLs(phi)
%INPUT
% phi: phase values in freqs X {mid|lowCL|highCL} format, in radians.
%OUTPUT
% unwrapped: data from <phi> unwrapped such that the CLs are still the same
%   distance from the unwrapped "mid" column as in the original data.
%OPTIONS
% 'degrees' - converts data from degrees to radians at entry, and back
%   again to degrees at exit.

%$Rev: 154 $
%$Date: 2012-07-24 19:04:58 -0400 (Tue, 24 Jul 2012) $
%$Author: dgibson $

argnum = 1;
degreesflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'degrees'
            degreesflag = true;
        otherwise
            error('dg_unwrapPhaseCLs:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if degreesflag
    phi = 2 * pi * phi / 360;
end

CLs = [phi(:,1)-phi(:,2) phi(:,3)-phi(:,1)];

unwrapped = unwrap(phi(:,1));
unwrapped(:,2) = unwrapped(:,1) - CLs(:,1);
unwrapped(:,3) = unwrapped(:,1) + CLs(:,2);

if degreesflag
    unwrapped = 360 * unwrapped / (2*pi);
end
