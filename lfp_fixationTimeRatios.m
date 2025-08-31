function ratios = lfp_fixationTimeRatios(varargin)
%ratios = lfp_fixationTimeRatios(trials, fixations, seqs, targIDs1, targIDs2)
%ratios = lfp_fixationTimeRatios(trials, seqs, targIDs1, targIDs2)
%   For each trial in <trials>, calculates the ratio of the total amount of
%   time spent fixating the targets in <targIDs1> to that spent fixating
%   targets in <targIDs2>.  <fixations> is as returned by
%   lfp_EyeTabulation, and <seqs> is as returned by lfp_TargSeq.  <trials>
%   here does NOT necessarily contain trial numbers, but rather contains
%   indices into <fixations> and <seqs>. Those indices will be the same the
%   same as trial numbers ONLY IF <fixations> and <seqs> were calculated
%   for ALL trials.  If <trials> is [], then all trials <fixations> and
%   <seqs> in will be used.  The target IDs in <targIDs1> and <targIDs2>
%   may overlap to any degree desired.  Note that it is possible for there
%   to zero time spent fixating <targIDs2>, in which case a divide-by-zero
%   warning will be raised and Inf or NaN returned.  NaN is returned if
%   <seqs> or <fixations> is empty for that trial.  <ratios> is a column
%   vector.
%   The 5-argument form is an archaic one for the pre-7/25/05 2-column
%   format of <seqs>.  The 4-argument form is preferred.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

if length(varargin) == 5
    [trials, fixations, seqs, targIDs1, targIDs2] = deal(varargin{:});
    if numel(fixations) ~= numel(seqs)
        error('lfp_fixationTimeRatios:badLengths', ...
            '<fixations> and <seqs> must contain the same trials' );
end
elseif length(varargin) == 4
    [trials, seqs, targIDs1, targIDs2] = deal(varargin{:});
else
    error('Type "help lfp_fixationTimeRatios" for calling syntax.');
end

if isempty(trials)
    trials = 1:numel(seqs);
end

ratios = repmat(NaN, numel(trials), 1);
for trialidx=1:numel(trials)
    seqtable = seqs{trials(trialidx)};
    if length(varargin) == 5
        fixtable = fixations{trials(trialidx)};
        seqtable = seqtable(:,1:2);
        if ~isempty(fixtable) && ~isempty(seqtable)
            % Inner-join the two tables on start time; any pairs found by dg_zip where
            % the start times differ by more than .5 microsecond must be different
            % fixations.  The end result is to add the target IDs as a sixth column
            % to the data in fixtable.
            pairs = dg_zip(fixtable(:,3), seqtable(:,1));
            remove = find(abs(fixtable(pairs(:,1), 3) - seqtable(pairs(:,2), 2)) < 0.5e-6);
            pairs(remove,:) = [];
            joinedtable = [fixtable(pairs(:,1), :) seqtable(pairs(:,2), 2)];
            fixtime1 = sum(joinedtable(find(ismember(joinedtable(:,6), targIDs1)), 4));
            fixtime2 = sum(joinedtable(find(ismember(joinedtable(:,6), targIDs2)), 4));
        end
    else
        seqidx1 = ismember(seqtable(:,2), targIDs1);
        seqidx2 = ismember(seqtable(:,2), targIDs2);
        fixtime1 = sum(seqtable(seqidx1, 3) - seqtable(seqidx1, 1));
        fixtime2 = sum(seqtable(seqidx2, 3) - seqtable(seqidx2, 1));
    end
    ratios(trialidx) = fixtime1/fixtime2;
end

