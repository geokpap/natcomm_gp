function [shufidx, numfailed] = dg_multiShuffle(numtr, numch, ...
    heuristic, getdesperate)
% Given <numch> channels of data containing <numtr> trials, generate a
% shuffled trial index set s.t. no two channels in any given shuffled trial
% are from the same original trial.  Depending on the relative numbers of
% channels and trials, this can be a difficult problem to solve by random
% search, which is effectively what we are doing here.  In the extreme case
% where numch = numtr, I believe that setting the values for the first
% channel determines that the only solutions are the permutations of the
% <numtr> rotations of column 1, which is a very small number compared to
% the space of all possible unconstrained shuffles, and when numch > numtr
% there are clearly no solutions at all.  Therefore, after <getdesperate>
% unsuccessful attempts at finding a solution randomly, we start using some
% compute-intensive heuristics to prune bad choices early in the search.
% After <maxfailures> unsuccessful attempts, we give up and throw an error.
%INPUTS
% numtr: number of trials
% numch: number of channels
% heuristic: how aggressively to apply heuristic pruning, which is
%   initially powerfully constrained by the ordering of the first column
%   and thus causes the second column to be a noisy version of the reverse
%   of the first.  0 (default): no heuristics. N > 0: start heuristics
%   after assigning trials to column <N> (so N=1 implies that the very
%   first set of randomly chosen trials will be pruned to avoid conflicts
%   with the perfectly sequential column 1).
% getdesperate: number of failures after which to apply even more
%   compute-intensive pruning heuristics which "should be" more successful.
%   A value of 0 indicates never to do that.  Has no effect when heuristic
%   = 0.
%OUTPUTS
% shufidx: a <numtr> X <numch> array of trial indices.  shufidx(:,1) is
%   always unshuffled, i.e. (1:<numtr>)'.  (If you don't like that, you can
%   shuffle the rows of <shufidx> after it's returned.)
% numrepeats: an integer that counts the total number of times across all
%   channels that the simple single-channel shuffling had to be repeated.
%NOTES
% Whenever it is possible to use dg_multiBoot instead of dg_multiShuffle,
% it runs much faster.

%$Rev: 191 $
%$Date: 2014-04-10 19:20:37 -0400 (Thu, 10 Apr 2014) $
%$Author: dgibson $

maxfailures = 50; % number of failed attempts after which we give up
% Random number generator initialization.  This value is hard coded for now
% because it should probably always be 1 for testing and 'shuffle' for
% actual use.
%   0: none (which means that the random number generator is always
%       initialized to the same value when starting Matlab, so the sequence
%       produced dg_multiShuffle depends on the history of the Matlab
%       process).
%   <nonnegative integer>: fixed constant initialization on entering
%       dg_multiShuffle.
%   'shuffle': seeds the random number generator based on the current time.
initalization = 'shuffle';

if nargin < 3
    heuristic = 0;
end
if ~isequal(initalization, 0)
    rng(initalization); 
end
shufidx = zeros(numtr, numch);
shufidx(:, 1) = (1:numtr)';
numfailed = 0;
broken = true;
while broken
    broken = false;
    shufidx = zeros(numtr, numch);
    shufidx(:, 1) = (1:numtr)';
    for ch = 2:numch
        for tr = 1:numtr
            availtr = setdiff( 1:numtr, ...
                [shufidx(tr, 1:ch-1) shufidx(1:tr-1, ch)'] );
            if getdesperate > 0 && numfailed >= getdesperate && ...
                    heuristic > 0 && ch > heuristic
                % Use desperate measures, i.e. an even stronger versions of
                % "help avoid painting ourselves into a corner".  Note that
                % these can be made yet stronger by constructing
                % second-preferred, third=preferred etc. lists when
                % <preferredtr> comes up empty.
                if tr < numtr
                    % Strategy #1: the more frequently a trial has been
                    % used already in all the later trials on all the
                    % preceding channels, the more we want to use it for
                    % the current trial.
                    freqused = histc( ...
                        reshape(shufidx(tr+1:end, 1:ch-1), 1, []), 1:numtr+1 );
                    [~, ix] = sort(freqused, 'descend');
                    topidx = ix(freqused(ix) == freqused(ix(1)));
                    % This kludge is necessary because at some point we have
                    % already assigned all of the top-scoring original trial
                    % numbers, and we must assign the rest
                    % arbitrarily:
                    preferredtr = intersect(availtr, topidx);
                    if ~isempty(preferredtr)
                        availtr = preferredtr;
                    end
                end
                % Strategy #2: the more we reduce the number of possible
                % solutions on later trials, the more likely we are to
                % eventually have no possible solutions.  Therefore, we
                % avoid any choices that will make the solution set even
                % smaller on trials that already have the smallest solution
                % sets.
                solutionsets = cell(numtr,1);
                setlengths = NaN(numtr,1);
                for latertr = tr+1:numtr
                    solutionsets{latertr} = setdiff(1:numtr,  ...
                        [shufidx(latertr, 1:ch-1) shufidx(1:tr-1, ch)'] );
                    setlengths(latertr) = length(solutionsets{latertr});
                end
                minlength = min(setlengths);
                minlntr = setlengths == minlength;
                % eliminate any of <availtr> that make the solution sets of
                % any of the <minlntr> smaller - i.e. any members of those
                % solution sets.
                reservedtr = unique(cell2mat(solutionsets(minlntr)'));
                preferredtr = setdiff(availtr, reservedtr);
                if ~isempty(preferredtr)
                    availtr = preferredtr;
                end
            elseif heuristic > 0 && ch > heuristic
                % To help avoid painting ourselves into a corner, we want to
                % choose trials early on that have already been used at later
                % positions on the preceding channels:
                preferredtr = intersect( availtr, ...
                    reshape(shufidx(tr+1:end, 1:ch-1), 1, []) );
                if ~isempty(preferredtr)
                    availtr = preferredtr;
                end
            end
            
            if isempty(availtr)
                % We have already chosen all the trials so far either on
                % this channel, or on preceding channels for this trial;
                % i.e. we have painted ourselves into a corner. Start all
                % over again.
                warning('dg_multiShuffle:gridlock', ...
                    'Shuffle became gridlocked; trying again');
                numfailed = numfailed + 1;
                if numfailed == getdesperate
                    warning('dg_multiShuffle:desperate', ...
                        'invoking desperate measures');
                end
                broken = true;
                break
            end
            % Matlab doc states that the result of 'rand' is "drawn from the
            % standard uniform distribution on the open interval (0,1)":
            idx = ceil(rand * length(availtr));
            shufidx(tr, ch) = availtr(idx);
        end
        if broken
            if numfailed >= maxfailures
                error('dg_multiShuffle:darn', ...
                    'ch=%d, tr=%d, %d failures; I give up!', ...
                    ch, tr, numfailed);
            end
            break
        end
    end
end
