function s2 = lfp_RecoverTargSeq(s1)
%LFP_RECOVERTARGSEQ renumbers ordered-target sequences to reflect the 
% spatial location of targets that the monkey is intended to capture.
% (Joey-specific function)

%s2 = lfp_RenumberTargSeq(s1)
%   <s1> and <s2> are both cell vectors of target fixations in the same
%   format as is returned by lfp_TargSeq.  Any fixation on a target which
%   is part of the sequence that the monkey is intended to capture is given
%   a new target number equal to <offset> plus the target's position in the
%   intended sequence (so if <offset> is 10, then the first target is 11,
%   the second is 12, etc.).  <offset> should be chosen so that the
%   resulting values do not conflict with any of the target ID numbers in
%   s1.  The renumbered sequence table is returned in s2.  It is assumed
%   that lfp_SelectedTrials has the same value it did when
%   lfp_EyeTabulation was run; if its length is different, then an error is
%   signalled, but it could have the same length and still be different.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if length(trials) ~= length(s1)
    error('lfp_RenumberTargSeq:badSelect', ...
        'Trial selection state has changed since generating fixations table' );
end

s2 = cell(size(s1));

for trialidx = 1:length(s1)
    if ~isempty(s1{trialidx})
        trial = trials(trialidx);
        lfp_getTaskParams;
        fixations = s1{trialidx}(:,2);
        fixations(find(fixations == 11)) = Target1;
        fixations(find(fixations == 12)) = Target2;
        fixations(find(fixations == 13)) = Target3;
        s2{trialidx} = [s1{trialidx}(:,1) fixations];
    end
end