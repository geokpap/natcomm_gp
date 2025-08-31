function [trialpairs, ambiguous, unmatched, imperfect, ...
    ctxtimescale] = lfp_assembleCtxTrialMap(ctxdata)
%[trialpairs, ambiguous, unmatched] = lfp_assembleCtxTrialMap(ctxdata)
% Invoke lfp_getCtxTrialMap iteratively starting in the middle of the ctx
% session and working outwards until the entire session is mapped.  Any
% runs of consecutive unmapped ctx trials at beginning or end of session
% are ignored.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

t = cputime; 
trialpairs = zeros(0, 2);
ambiguous = cell(0, 2);
unmatched = cell(0, 2);
imperfect = [];
ctxtimescale_est = [];
batchsize = 50;
ctxstart = fix(length(ctxdata.codes)/2);
ctxend = min(ctxstart+batchsize-1, length(ctxdata.codes));
lfp_log(sprintf('Finding seed matches %d:%d -> %d:%d.', ...
    ctxstart, ctxend, 1, size(lfp_TrialIndex,1) ));
[p, a, u, i, c] = ...
    lfp_getCtxTrialMap(ctxdata, ctxstart:ctxend, 1:size(lfp_TrialIndex,1));
% Unless there is a substantial fraction of trials matched, the matches are
% most likely false positives:
if size(p, 1) > batchsize/4
    trialpairs = p;
    ambiguous = a;
    unmatched = u;
    imperfect = i;
    ctxtimescale_est = c;
end

lfp_log('done.');

% Match end of session in batches
lfpend = size(lfp_TrialIndex,1);
while ctxend < length(ctxdata.codes)
    ctxstart= ctxend + 1;
    ctxend = min(ctxstart+batchsize-1, length(ctxdata.codes));
    if isempty(trialpairs)
        lfpstart = 1;
    else
        lfpstart = trialpairs(end, 2) + 1;
    end
    lfp_log(sprintf('Finding tail matches %d:%d -> %d:%d.', ...
        ctxstart, ctxend, lfpstart, lfpend ));
    [p, a, u, i, c] = lfp_getCtxTrialMap(ctxdata, ...
        ctxstart:ctxend, lfpstart:lfpend);
    lfp_log('done.');
    if size(p, 1) > (lfpend-lfpstart+1)/4
        trialpairs = [trialpairs; p];
        ambiguous = [ambiguous; a];
        unmatched = [unmatched; u];
        imperfect = [imperfect; i];
        ctxtimescale_est = [ctxtimescale_est; c];
    end
end

% Match beginning of session in batches
lfpstart = 1;
ctxstart = fix(length(ctxdata.codes)/2);
while ctxstart > 1
    ctxend= ctxstart - 1;
    ctxstart = max(ctxend-batchsize+1, 1);
    if isempty(trialpairs)
        lfpend = size(lfp_TrialIndex,1);
    else
        lfpend = trialpairs(1, 2) - 1;
    end
    lfp_log(sprintf('Finding head matches %d:%d -> %d:%d.', ...
        ctxstart, ctxend, lfpstart, lfpend ));
    [p, a, u, i, c] = lfp_getCtxTrialMap(ctxdata, ...
        ctxstart:ctxend, lfpstart:lfpend);
    lfp_log('done.');
    if size(p, 1) > (lfpend-lfpstart+1)/4
        trialpairs = [p; trialpairs];
        ambiguous = [a; ambiguous];
        unmatched = [u; unmatched];
        imperfect = [i; imperfect];
        ctxtimescale_est = [c; ctxtimescale_est];
    end
end

ctxtimescale = mean(ctxtimescale_est);
lfp_log(sprintf('ctxtimescale = %d s, 2*SEM = %.0d %%', ...
    ctxtimescale, ...
    100 * 2 * std(ctxtimescale_est) / ...
    (ctxtimescale * sqrt(length(ctxtimescale_est))) ));

% Unmatchable trials returned by dg_resolveAmbig do not have any match info
% associated with them, so use empty elements in col. 2:
[trialpairs, ambiguous, u] = dg_resolveAmbig(trialpairs, ambiguous);
unmatched = [unmatched; u(:,1) cell(size(u,1), 1)];
[a, ix] = sort(cell2mat(unmatched(:,1)));
unmatched = unmatched(ix, :);

% We can safely assume that any runs of missing trial matches that are
% preceded and followed by the same offset are also matched at the same
% offset, because ctx trials are the reference standard and lfp trials
% cannot be inserted (so there cannot be balanced insertions and
% deletions). 
offset = trialpairs(:,2) - trialpairs(:,1);
rmun = [];
rmamb = [];
for gapidx = find((trialpairs(2:end, 1) ~= trialpairs(1:end-1,1) + 1)')
    % gapidx points to the pair before the gap.
    if (offset(gapidx+1) == offset(gapidx))
        ctxtrials = ...
            ((trialpairs(gapidx,1) + 1) : (trialpairs(gapidx+1,1) - 1))';
        lfptrials = ...
            ((trialpairs(gapidx,2) + 1) : (trialpairs(gapidx+1,2) - 1))';
        trialpairs = [ trialpairs
            ctxtrials lfptrials
            ];
        rmun = [ rmun
            find(ismember(cell2mat(unmatched(:,1)), ctxtrials)) ];
        rmamb = [ rmamb
            find(ismember(cell2mat(ambiguous(:,1)), ctxtrials)) ];
    end
end
trialpairs = sortrows(trialpairs);
unmatched(rmun, :) = [];
ambiguous(rmamb, :) = [];

elapsed = cputime-t;
lfp_log(sprintf('Finished in %.2f min.', elapsed/60));
