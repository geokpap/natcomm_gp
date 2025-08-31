function lfp_savePersistentGlobals

%$Rev: 56 $
%$Date: 2009-04-10 13:51:07 -0400 (Fri, 10 Apr 2009) $
%$Author: dgibson $

lfp_declareGlobals;
filename = 'lfp_PersistentGlobalValues.mat';
path = which(filename);
if isempty(path)
    path = filename;
end
save(path, 'lfp_DataDir', 'lfp_XLimAll', 'lfp_ClipBoard', ...
    'lfp_FreqLim', 'lfp_CLimAll', 'lfp_YLimAll', 'lfp_AlignmentRef', ...
    'lfp_XTicksOnAll', 'lfp_NoWaitbar');
