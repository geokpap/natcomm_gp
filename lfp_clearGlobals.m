%lfp_clearGlobals
% A script that clears all the variables listed in the file
% 'lfp_globals.mat', except for 'lfp_GlobalNames'.  See
% lfp_declareGlobals.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

load('lfp_globals.mat');
lfp_clearGlobals_globals_string = '';
for lfp_clearGlobals_names = lfp_declareGlobals_global_names
    lfp_clearGlobals_globals_string = ...
        [lfp_clearGlobals_globals_string ' ' char(lfp_clearGlobals_names)];
end
lfp_clearGlobals_cmd = ['clear global ' lfp_clearGlobals_globals_string];
eval(lfp_clearGlobals_cmd);
lfp_GlobalNames = lfp_declareGlobals_global_names;
clear lfp_clearGlobals_globals_string lfp_clearGlobals_names ...
lfp_clearGlobals_cmd lfp_declareGlobals_global_names