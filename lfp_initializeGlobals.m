%lfp_initializeGlobals
% A script that sets all the variables listed in the file
% 'lfp_globals.mat' to null ([]).  See lfp_declareGlobals.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

load('lfp_globals.mat');
lfp_initializeGlobals_globals_string = '';
for lfp_initializeGlobals_name = lfp_declareGlobals_global_names
    lfp_initializeGlobals_cmd = [ char(lfp_initializeGlobals_name) ...
            '= [];' ];
    eval(lfp_initializeGlobals_cmd);
end
lfp_GlobalNames = lfp_declareGlobals_global_names;
clear lfp_initializeGlobals_globals_string lfp_initializeGlobals_name ...
lfp_initializeGlobals_cmd lfp_declareGlobals_global_names