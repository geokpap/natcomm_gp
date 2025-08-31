function [match, total] = dg_scoreTrials(filepath, condition)
%   Reads a Yasuo-format events file, returns number of trials that satisfy
%   condition and total number of trials.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

[header, trialdata] = dg_readYasuoEvents(filepath);
match = 0;
total = length(trialdata);
condstr = dg_ParseSelectionSpec2(condition);
for trialnum = 1:length(trialdata)
    try
        trialmatches = eval(condstr);
    catch
        disp(lasterror);
        disp(sprintf('Aborting dg_scoreTrials for "%s" %s', ...
            filepath, condition));
    end
    if trialmatches
        match = match + 1;
    end
end