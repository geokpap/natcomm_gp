function lfp_addSeqs2DB(db, seqs)
% Adds the target IDs in <seqs> to the 'fixations' table in database <db>.
% <seqs> is as returned by lfp_TargSeq, <db> is as required by dg_runSQL.
% For large numbers of fixations, a significant time savings will be
% realized if <db> is a connection object rather than a string.  The
% 'fixations' table in <db> must already contain a numeric 'targID' column.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if length(trials) ~= numel(seqs)
    error('lfp_addSeqs2DB:badSelect', ...
        'Trial selection state has changed since generating <seqs>' );
end

for trialidx = 1:numel(seqs)
    trialnum = trials(trialidx);
    trialID = lfp_getTrialID(trialnum);
    for seqrow = 1:size(seqs{trialidx}, 1)
        dg_runSQL(db, sprintf([ 'UPDATE fixations SET targID = %d ' ...
            'WHERE TrialID = ''%s'' AND Start = %.6f' ], ...
            seqs{trialidx}(seqrow, 2), trialID, seqs{trialidx}(seqrow, 1) ));
    end
end
