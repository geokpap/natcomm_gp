function [p, ambiguous, unmatchable] = dg_resolveAmbig(p, a)
%[pairs, ambiguous] = lfp_resolveCtxAmbig(trialpairs, a, u)
% Uses the no-crossover constraint to eliminate ambiguous trial matches
% where possible.  <p> is a list of pairs like <trialpairs> from
% lfp_getCtxTrialMap. <ambiguous>, <a> like <ambiguous> from
% lfp_getCtxTrialMap.  Repeats the basic ambiguity removal procedure until
% no more ambiguities can be removed.  <unmatchable> contains any rows
% from a where all the possibilities would involve crossovers.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

ambiguous = a;
old_a = {};
unmatchable = cell(0,2);
while ~isequal(ambiguous, old_a)
    removedidx = [];
    for ambidx = 1:size(ambiguous,1)
        predecessoridx = find(p(:,1) < ambiguous{ambidx, 1});
        if ~isempty(predecessoridx)
            predecessoridx = predecessoridx(end);
        end
        successoridx = find(p(:,1) > ambiguous{ambidx, 1});
        if ~isempty(successoridx)
            successoridx = successoridx(1);
        end
        possibilities = ambiguous{ambidx, 2};
        if ~isempty(predecessoridx)
            possibilities(possibilities <= p(predecessoridx, 2)) = [];
        end
        if ~isempty(successoridx)
            possibilities(possibilities >= p(successoridx, 2)) = [];
        end
        if length(possibilities) == 1
            p = sortrows([p; ambiguous{ambidx, 1} possibilities]);
            removedidx(end+1) = ambidx;
        elseif length(possibilities) == 0
            unmatchable(end+1,:) = ambiguous(ambidx, :);
            removedidx(end+1) = ambidx;
        end
    end
    old_a = ambiguous;
    ambiguous(removedidx, :) = [];
end
