function lfp_showSetup
%LFP_SHOWSETUP displays the current setup name in the command window.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

global lfp_SetupName

% Get setup name:
if isempty(lfp_SetupName)
    if exist('lfp_SetupName.mat') == 2
        load('lfp_SetupName.mat');
        disp(['Using ' lfp_SetupName ' setup']);
    else
        disp('No setup name saved');
    end
else
    disp(['Using ' lfp_SetupName ' setup']);
end
