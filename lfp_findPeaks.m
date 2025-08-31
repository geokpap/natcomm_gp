function vmaxidx = lfp_findPeaks(fn, thresh1, thresh2)
%vmaxidx = lfp_findPeaks(fn, thresh1, thresh2)
%   <thresh2> is optional. Returns a column vector of indices into
%   lfp_Samples{fn} of 3-point maxima that are greater than <thresh1> and
%   (if given) less than <thresh2>.  The middle point must be strictly
%   greater than either of its neighbors to qualify as a 3-point maximum.

%$Rev: 123 $
%$Date: 2010-05-06 19:56:19 -0400 (Thu, 06 May 2010) $
%$Author: dgibson $

lfp_declareGlobals;
allocsize = 2^14;

% Start by qualifying small regions of the trace to search for maxima:
qualifiedpts(:,1) = reshape(lfp_Samples{fn}(:) > thresh1, [], 1);
% Since the first and last point of each 3-point maximum do not have to
% meet the minimum threshold criterion, we must expand the qualified
% regions by one point on each side:
qualifiedpts = qualifiedpts | [0; qualifiedpts(1:end-1)] ...
    | [ qualifiedpts(2:end); 0] ;
if nargin > 2
    qualifiedpts = qualifiedpts & ...
        reshape(lfp_Samples{fn}(:) < thresh2, [], 1);
end
qualstarts = [ qualifiedpts(1)
    qualifiedpts(2:end) & ~qualifiedpts(1:end-1) ];
qualends = [ qualifiedpts(1:end-1) & ~qualifiedpts(2:end)
    qualifiedpts(end) ];
startidx = find(qualstarts);  % index into lfp_Samples{fn}
endidx = find(qualends);  % index into lfp_Samples{fn}
if length(startidx) ~= length(endidx)
    error('lfp_findPeaks:badqual', ...
        'Internal error - call the geek squad!');
end

% Now search for 3-point maximum values:
vmaxidx = zeros(allocsize,1);
vmaxidxptr = 0;
for k = 1:length(startidx)
    ix1 = startidx(k);
    ix2 = endidx(k);
    localmaxidx = find([ false
        reshape(lfp_Samples{fn}(ix1+1:ix2-1), [], 1) > ...
        reshape(lfp_Samples{fn}(ix1:ix2-2), [], 1) ...
        & reshape(lfp_Samples{fn}(ix1+1:ix2-1), [], 1) > ...
        reshape(lfp_Samples{fn}(ix1+2:ix2), [], 1)
        false ]);
    if vmaxidxptr + length(localmaxidx) > size(vmaxidx,1)
        vmaxidx = [ vmaxidx; zeros(allocsize,1) ];
    end
    vmaxidx(vmaxidxptr + (1:length(localmaxidx))) = ...
        localmaxidx + ix1 - 1;   % sample index
    vmaxidxptr = vmaxidxptr + length(localmaxidx);
end
vmaxidx(vmaxidxptr + 1 : end, :) = [];
