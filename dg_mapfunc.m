function result = dg_mapfunc(fhandle, varargin)
% result = dg_mapfunc(fhandle, arg1, arg2, ...)
%The function denoted by <fhandle> must take as many arguments as are
%provided after <fhandle>; at least one argument must be provided. All args
%must be of the same size, and all must be cell arrays. The result is a
%cell array of the same size as the args, such that element j is the result
%of applying <fhandle> to element j of each of the arguments.
%NOTE: Matlab's "cellfun" is probably better.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

result = cell(size(varargin{1}));
args = cell(size(varargin));
for k = 1:numel(varargin{1})
    for a = 1:length(varargin)
        args{a} = varargin{a}{k};
    end
    result{k} = feval(fhandle, args{:});
end