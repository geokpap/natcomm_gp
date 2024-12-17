function lfp_changeSetup(name)
%lfp_changeSetup(name)
% Changes the name of the file that is used for setup.  If <name> is [],
% then lfp_setup.m is used.  Otherwise, lfp_setup_<name>.m is used.  For
% example,
%   lfp_changeSetup('bill')
% causes the file lfp_setup_bill.m to be run.

%$Rev: 42 $
%$Date: 2008-12-29 18:09:15 -0500 (Mon, 29 Dec 2008) $
%$Author: dgibson $

global lfp_SetupName

lfp_SetupName = name;
filename = 'lfp_SetupName.mat';
path = which(filename);
if isempty(path)
    path = filename;
end
save(path, 'lfp_SetupName');
