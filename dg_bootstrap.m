function [bootsam, varargout] = dg_bootstrap(nretn, nboot, bootfun, ...
    varargin)
%   A simpler but more complicated version of 'bootstrp'.
% WARNING: uses 'parfor', may be bizarrely tweaky.
%INPUTS
% nretn: number of return values to accept from bootfun
% nboot: number of bootstrap sets to run
% bootfun: function handle to run on each bootstrap set; its first argument
%   must be the one that is getting bootstrapped, and it gets bootstrapped
%   over its rows (so if it is a vector, it must be a column vector).
% varargin: all of the arguments to pass to <bootfun>.  The first one gets
%   resampled to form the bootstrap sets, and the remainder are passed in
%   verbatim.  As in 'bootstrp', the sampling is done on a row-by-row
%   basis.
%OUTPUTS
% bootsam:  Each column in <bootsam> contains indices of the values that
%   were drawn from varargin{1} (i.e. the first argument to <bootfun> to
%   constitute the corresponding bootstrap sample. 
% v1, ..., v<nretn>:  <nretn> values are returned, each of which is a cell
%   column vector containing <nboot> instances of one of the values
%   returned by <bootfun>.  <v1> contains the first return value, v<nretn>
%   contains the last (i.e., the <nretn>th) value.  The <nboot> instances
%   are in the same order for each return value, i.e. v<k>{bootidx} is the
%   kth return value for bootstrap set number <bootidx>.

%$Rev: 305 $
%$Date: 2023-09-08 14:19:01 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

numrows = size(varargin{1},1);
bootsam = zeros(numrows, nboot);
bootdata = varargin{1};
otherargs = varargin(2:end);
parfor bootnum = 1:nboot
    bootsam(:,bootnum) = ceil(numrows * rand(numrows,1));
    retnvals = cell(1, nretn);
    [retnvals{:}] = feval(bootfun, ...
        bootdata(bootsam(:,bootnum)), otherargs{:}); %#ok<PFBNS>
    retnvalarray{bootnum} = retnvals;
end
% retnvalarray is now a cell vector of length <nboot>, each element of
% which is a cell vector of length <nretn>.  That is exactly inside out of
% how we want varargout.
varargout = repmat(cell(nboot,1), 1, nretn);
for valnum = 1:nretn
    for bootnum = 1:nboot
        varargout{valnum}{bootnum,1} = retnvalarray{bootnum}{valnum};
    end
end

