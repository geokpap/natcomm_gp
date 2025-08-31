function dg_mkdir_p(dirpath)
% Equivalent to unix "mkdir -p".  Creates whatever intermediate directories
% are required in order to create <dirpath>.

%$Rev: 288 $
%$Date: 2022-02-04 17:39:43 -0500 (Fri, 04 Feb 2022) $
%$Author: dgibson $

parent = fileparts(dirpath);
if ~exist(parent, 'dir')
    dg_mkdir_p(parent);
end
mkdir(dirpath);

