function table2 = lfp_restoreEmptyTrials(table, trialIDs)
%lfp_restoreEmptyTrials
%   In dg_readEyeTab results, restores empty cells to that were discarded
%   when writing a fixations or saccades <table> to a text file.  This is
%   done on the basis of the assumption that the trial selection state is
%   currently the same as it was when the table was written to the file,
%   and that the trials are stored in the same order in memory as when the
%   table was written. There is no way to check the validity of these
%   assumptions (except that if <table> contains more entries than there
%   are selected trials, then the first assumption is clearly false).
%   <table> and <trialIDs> are as returned by dg_readEyeTab.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if length(trials) < length(table)
    error('lfp_restoreEmptyTrials:badSelect', ...
        'Trial selection state has changed since generating the table' );
end

table2 = {};
trialIDidx = 1;
for trial = trials
    if trialIDidx >length(trialIDs)
        length2 = length(table2);
        table2(1, (length2+1):length(trials)) = ...
            cell(1, length(trials) - length2);
        break
    end
    if isequal(lfp_getTrialID(trial), trialIDs{trialIDidx})
        table2(1, trial) = table(trialIDidx);
        trialIDidx = trialIDidx + 1;
    else
        table2(1, trial) = {[]};
    end
end

    