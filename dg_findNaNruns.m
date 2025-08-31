function [nanidx nonnanidx] = dg_findNaNruns(A)
%[nanidx nonnanidx] = dg_findNaNruns(A)
%   Returns two vectors: <nanidx> contains the indices of the first NaN in
%   each run of consecutive NaNs in <A>, and <nonnanidx> contains the
%   indices of the first non-NaN in each run of consecutive non-NaNs.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

nanidx = [];
if ~isnan(A(1))
    nonnanidx = 1;
else
    nonnanidx = [];
end

startpt = 1;
while startpt <= numel(A)
    nans = find(isnan(A(startpt:end)));
    if ~isempty(nans)
        nanidx(end+1) = startpt + nans(1) - 1;
        % find next non-NaN
        nonnan = find(~isnan(A( nanidx(end)+1 : end )));
        if ~isempty(nonnan)
            nonnanidx(end+1) = nanidx(end) + nonnan(1);
            startpt = nonnanidx(end);
        else
            break;
        end
    else
        break;
    end
end
