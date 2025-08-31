function tapers=dpsschk(tapers,N)
% Helper function to calculate tapers and, if precalculated tapers are supplied, 
% to check that they (the precalculated tapers) the same length in time as
% the time series being studied. The length of the time series is specified
% as the second input argument N. Thus if precalculated tapers have
% dimensions [N1 K], we require that N1=N.
% Usage: tapers=dpsschk(tapers,N)
% Inputs:
% tapers        (tapers in the form of: 
%                                   (i) precalculated tapers or,
%                                   (ii) [NW K] - time-bandwidth product, number of tapers) 
% N             (number of samples)
% Outputs: 
% tapers        (calculated or precalculated tapers)
if nargin < 2; error('Need all arguments'); end
sz=size(tapers);
if sz(1)==1 & sz(2)==2;
    tapers=dpss(N,tapers(1),tapers(2));
elseif N~=sz(1);
    error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
end;