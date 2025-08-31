function [trialpairs, ambiguous, unmatched, imperfect, ...
    ctxtimescale_est] = lfp_getCtxTrialMap(...
    ctxdata, ctxtrials, lfptrials, showfigs)
%[trialpairs, ambiguous, unmatched, imperfect, ...
%     ctxtimescale_est] = lfp_getCtxTrialMap(...
%     ctxdata, ctxtrials, lfptrials, showfigs)
% Finds correct mapping between Cortex trial numbers and lfp_lib trial
% numbers by matching longest common subsequences of event IDs and then
% looking for matching inter-event intervals.  <ctxdata> is the entire
% cortex file as returned by Spiketools' cortex_read function.  <ctxtrials>
% is the range of trial numbers in the Cortex file that will be matched (if
% possible) to the range <lfptrials> in lfp_lib.  If <showfigs> is true,
% then histograms of the summary stats are displayed.
%
% Preprocessing of Cortex data such as removal of spike events and events
% #0 should be done first. <trialpairs> is a two-column array of cortex
% trial numbers in col. 1 and their matched counterparts in lfp_lib in col.
% 2.  It is in sorted order in both columns.  <ambiguous> is a 2-column
% cell array with cortex trial numbers in col. 1 and lists of potential
% matching counterparts in lfp_lib in col. 2.  <unmatched> contains ctx
% trial numbers and intermediate results for trials that didn't match at
% all, as a cell array with col 1 = ctxtrial, col 2 = structure with
% fields: 
%   lfptrials
%   fracs
%   diff
%   ctxints
%   lfpints
%   evtpairs
% <imperfect>, <unmatched> and <ambiguous> are sorted by cortex trial
% number (col 1).  <ctxtimescale_est> is a vector of estimates of
% ctxtimescale calculated from each of <trialpairs>.
%NOTE:
% This function assumes that the matching session is already loaded in
% lfp_lib.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

if nargin < 4
    showfigs = false;
end

trialpairs = zeros(0,2);
ambiguous = cell(0,2);
unmatched = cell(0,2);
imperfect = cell(0,2);
ctxtimescale_est = zeros(0,1);

if isempty(ctxtrials) || isempty(lfptrials)
    return
end

lfp_declareGlobals;

ctxtrials = reshape(ctxtrials, 1, []);

maxcode = 115;  % maximum Cortex event code value
match = diag(ones(maxcode, 1));
match(match==0) = -Inf;
indel = zeros(maxcode,1);

% Find best sets of event matches between ctx trials and nlx trials based
% on event ID alone:
sim = zeros(length(ctxdata.codes), size(lfp_TrialIndex,1));
evtpairs = cell(length(ctxdata.codes), size(lfp_TrialIndex,1));
for ctxtrl = ctxtrials
    for lfp_trial = lfptrials
        [sim(ctxtrl, lfp_trial), evtpairs{ctxtrl, lfp_trial}] = ...
            dg_seqAlign3(ctxdata.codes{ctxtrl}, ...
            lfp_Events(lfp_TrialIndex(lfp_trial,1):lfp_TrialIndex(lfp_trial,2), 2));
    end
end

% Calculate inter-event intervals for matched events, and calculate
% differences in inter-event intervals between ctx trials and lfp trials
% (differences of more than one Cortex clock tick ~ 1 ms indicate that at
% least one event pairing is wrong):
ctxtimescale = 9.7652e-4;  % measured on session 2006-5-4_14-45-37
diff = cell(size(evtpairs));
ctxints = cell(size(evtpairs));
lfpints = cell(size(evtpairs));
for ctxtrl = ctxtrials
    for lfp_trial = 1:size(lfp_TrialIndex,1)
        if ~isempty(evtpairs{ctxtrl, lfp_trial})
            ctxidx = evtpairs{ctxtrl, lfp_trial}(:,1);
            lfpidx = evtpairs{ctxtrl, lfp_trial}(:,2);
            ctxints{ctxtrl, lfp_trial} = ( ctxdata.codetimes{ctxtrl}(ctxidx(2:end)) ...
                - ctxdata.codetimes{ctxtrl}(ctxidx(1:end-1)) ) * ctxtimescale;
            lfpints{ctxtrl, lfp_trial} = ...
                lfp_Events(lfpidx(2:end) + lfp_TrialIndex(lfp_trial,1) - 1, 1) ...
                - lfp_Events(lfpidx(1:end-1) + lfp_TrialIndex(lfp_trial,1) - 1, 1);
            diff{ctxtrl, lfp_trial} = lfpints{ctxtrl, lfp_trial} - ctxints{ctxtrl, lfp_trial};
        end
    end
end

% Find the ctx-lfp trial evtpairs with acceptably low inter-event interval
% differences:
numdiffs = NaN(length(ctxtrials), size(lfp_TrialIndex,1));
smalldiffs = NaN(length(ctxtrials), size(lfp_TrialIndex,1));
for ctxtrl = ctxtrials
    for lfp_trial = 1:size(lfp_TrialIndex,1)
        smalldiffs(ctxtrl, lfp_trial) = sum(abs(diff{ctxtrl, lfp_trial}) < .001);
        numdiffs(ctxtrl, lfp_trial) = length(diff{ctxtrl, lfp_trial});
    end
end

% For each ctx trial, find the lfp trial with the highest fraction of
% well-matched inter-event intervals.  Use total number of paired events to
% break ties.  Put any remaining ties on the "ambiguous" list.  Give up
% entirely on any ctx trials where the best fraction is below 80%.
fracs = NaN(length(ctxtrials), size(lfp_TrialIndex,1));
for ctxtrl = ctxtrials
    s = warning('off', 'MATLAB:divideByZero');
    fracs(ctxtrl,:) = smalldiffs(ctxtrl, :)./numdiffs(ctxtrl, :);
    s = warning(s.state, 'MATLAB:divideByZero');
    bestfrac = max(fracs(ctxtrl,:));
    if bestfrac < 1
        matchinfo.lfptrials = lfptrials;
        matchinfo.fracs = fracs(ctxtrl,:);
        matchinfo.diff = diff(ctxtrl, :);
        matchinfo.ctxints = ctxints(ctxtrl, :);
        matchinfo.lfpints = lfpints(ctxtrl, :);
        matchinfo.evtpairs = evtpairs(ctxtrl, :);
        if bestfrac < 0.8
            unmatched{end+1,1} = ctxtrl;
            unmatched{end, 2} = matchinfo;
            continue
        else
            imperfect{end+1,1} = ctxtrl;
            imperfect{end,2} = matchinfo;
        end
    end
    bfix = find(fracs(ctxtrl,:) == bestfrac);
    lfp_trial = [];
    if length(bfix) > 1
        bestsim = max(sim(ctxtrl, bfix));
        bsix = find(sim(ctxtrl, bfix) == bestsim);
        if length(bsix) > 1
            ambiguous(end+1, :) = {ctxtrl bfix(bsix)};
        else
            lfp_trial = bfix(bsix);
        end
    else
        lfp_trial = bfix;
        bestsim = sim(ctxtrl, lfp_trial);
    end
    if lfp_trial
        trialpairs(end+1, :) = [ctxtrl, lfp_trial];
    end
end

% Remove mismatches according to the heuristic that we expect long runs of
% trials to be matched in sequence, so that the offset from ct trial num to
% lfp trial num remains constant.  (Note that there may be some trials
% missing, but it will be the same number - usually just one - in both
% sequences of trial nums.)  Therefore, any trial pair whose offset is not
% equal to the preceding offset nor the succeeding offset is considered a
% mismatch.  (There may be step changes in offset, so equality is not
% required on both sides.)  The ends are special cases because they have
% only one neighbor, and are more prone to mismatches than interior
% trialpairs.  Also, since the ctx trials are guaranteed to be in monotonic
% increasing order, the matched nlx trials must also be in increasing
% order; therefore it would in principle be a mismatch if the following
% pair has a nlx trial number that is less than or equal to the nlx trial
% number of the current pair.  However, life gets more complicated, because
% a single erroneous match to an earlier nlx trial at the end of a long
% run of constant offset pairs would cause the entire long run to be
% invalidated.  What we want instead is majority vote logic, where the
% pair whose offset is in the majority wins.
done = false;
while ~done && size(trialpairs, 1) > 2
    offset = trialpairs(:,2) - trialpairs(:,1);
    % First the offset criterion:
    mismatch = [ offset(1) ~= offset(2)
        offset(2:end-1) ~= offset(1:end-2) & offset(2:end-1) ~= offset(3:end)
        offset(end) ~= offset(end-1) ];
    % OR in the majority-vote monotonicity criterion:
    nonmono = [ trialpairs(2:end,2) <= trialpairs(1:end-1,2) ];
    losers = logical(zeros(size(mismatch)));
    for k = 1:length(nonmono)
        thisvote = sum(offset == offset(k));
        thatvote = sum(offset == offset(k+1));
        if thisvote < thatvote
            losers(k) = true;
        else
            losers(k+1) = true;
        end
    end
    mismatch = mismatch | losers ;

    % For each mismatched trial, check to see if there are any acceptable
    % matches with offsets that match the neighbors.  "Acceptable" means
    % (fracs>=0.8) because there could be several corrupted event IDs that got
    % mismatched.
    remove = [];
    for badpairidx = find(mismatch)'
        ctxtrl = trialpairs(badpairidx,1);
        lfp_candidates = find(fracs(ctxtrl,:) >= 0.8);
        targetoffsets = [];
        if badpairidx > 1
            targetoffsets = offset(badpairidx-1);
        end
        if badpairidx < length(offset)
            targetoffsets = [targetoffsets offset(badpairidx+1)];
        end
        candidateoffsets = lfp_candidates - ctxtrl;
        lfp_candidates(~ismember(candidateoffsets, targetoffsets)) = [];
        matchinfo.lfptrials = lfptrials;
        matchinfo.fracs = fracs(ctxtrl,:);
        matchinfo.diff = diff(ctxtrl, :);
        matchinfo.ctxints = ctxints(ctxtrl, :);
        matchinfo.lfpints = lfpints(ctxtrl, :);
        matchinfo.evtpairs = evtpairs(ctxtrl, :);
        if length(lfp_candidates) == 1
            % Good, replace the mismatch
            trialpairs(badpairidx, 2) = lfp_candidates;
            if fracs(ctxtrl,lfp_candidates) < 1
                imperfect{end+1,1} = ctxtrl;
                imperfect{end,2} = matchinfo;
            end
        else
            if length(lfp_candidates) == 0
                % No good; it's unmatchable
                remove(end+1) = badpairidx;
                unmatched{end+1,1} = ctxtrl;
                unmatched{end, 2} = matchinfo;
            else
                % replacement for mismatch is ambiguous
                remove(end+1) = badpairidx;
                ambiguous(end+1, :) = {ctxtrl lfp_candidates};
            end
        end
    end
    oldtrialpairs = trialpairs;
    trialpairs(remove,:) = [];
    if isequal(oldtrialpairs, trialpairs)
        done = true;
    end
end
if ~isequal(sort(trialpairs(:,1)), trialpairs(:,1)) ...
        || ~isequal(sort(trialpairs(:,2)), trialpairs(:,2))
    error('lfp_getCtxTrialMap:mismatchremoval', 'Internal error');
end

[a, ix] = sort(cell2mat(imperfect(:,1)));
imperfect = imperfect(ix, :);
[a, ix] = sort(cell2mat(unmatched(:,1)));
unmatched = unmatched(ix, :);
[a, ix] = sort(cell2mat(ambiguous(:,1)));
ambiguous = ambiguous(ix, :);

% Compute summary stats re: trial pairs
numevtpairs = NaN(size(trialpairs,1),1);
smalldifffracs = NaN(size(trialpairs,1),1);
ctxtimescale_est = zeros(size(trialpairs,1),1);
for pairidx = 1:size(trialpairs,1)
    ctxtrl = trialpairs(pairidx,1);
    lfp_trial = trialpairs(pairidx,2);
    numevtpairs(pairidx) = size(evtpairs{ctxtrl, lfp_trial}, 1);
    smalldifffracs(pairidx) = fracs(ctxtrl,lfp_trial);
    % Estimate ctxtimescale from the first and last well-matched event
    % pairs in each trial:
    smalldiffidx = find(abs(diff{ctxtrl, lfp_trial}) < .001);
    evtpairidx1 = smalldiffidx(1);
    evtpairidx2 = smalldiffidx(end) + 1;
    ctxTS1 = ctxdata.codetimes{ctxtrl}( ...
        evtpairs{ctxtrl, lfp_trial}(evtpairidx1, 1) );
    ctxTS2 = ctxdata.codetimes{ctxtrl}( ...
        evtpairs{ctxtrl, lfp_trial}(evtpairidx2, 1) );
    lfpTS1 = lfp_Events( ...
        evtpairs{ctxtrl, lfp_trial}(evtpairidx1, 2) ...
        + lfp_TrialIndex(lfp_trial,1) - 1, 1 );
    lfpTS2 = lfp_Events( ...
        evtpairs{ctxtrl, lfp_trial}(evtpairidx2, 2) ...
        + lfp_TrialIndex(lfp_trial,1) - 1, 1 );
    ctxtimescale_est(pairidx) = (lfpTS2- lfpTS1)/(ctxTS2 - ctxTS1);
end
s = warning('off', 'MATLAB:divideByZero');
disp(sprintf('Mean event pairs = %d', mean(numevtpairs)));
disp(sprintf('Mean fraction of small diffs = %d', mean(smalldifffracs)));
s = warning(s.state, 'MATLAB:divideByZero');
if showfigs
    figure;
    subplot(2,1,1);
    hist(numevtpairs, min(numevtpairs):max(numevtpairs));
    title('Number of Event Pairs');
    subplot(2,1,2);
    hist(smalldifffracs, ((75:100)-.5)/100);
    title('Fraction of Small Interevent Interval Differences');
end

