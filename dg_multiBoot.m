function shufidx = dg_multiBoot(numtr, numch)
% As the number of channels starts to become comparable to the number of
% trials (e.g. 20 channels and 2000 trials, or 1%), the problem of
% constructing a true multichannel shuffle (where no original trial is
% repeated on any individual channel and no two channels come from the same
% original trial on any shuffled trial) becomes difficult, slow, and
% cumbersome.  However, for purposes of constructing network analysis null
% hypothesis data sets the constraint that no trial be repeated in any
% individual is not really necessary.  As long as no two channels come from
% the same original trial, the proxy data set is a legitimate instantiation
% of the idea that the within-channel statistics remain (essentially) the
% same while the across-channel correlations are destroyed.  However, we
% can preserve the within-channel statistics a little better than a
% standard sample-and-replace bootstrap by shuffling once and then
% replacing any repeated original trial numbers with a random selection
% from original trials that have not yet been used on the current shuffled
% trial.

% First, estimate the running time for chunkyII (see "lfp_lib Maintenance
% Vol 8.doc", preceding "3/21/2014 9:13 PM".  Note that this model
% substantially overestimates for smallish numbers of channels (<200),
% getting worse for largish numbers of trials (>1000).  The ratio of model
% estimation to actual running time on chunkyII at several points:
%   numch   numtr   ratio
%   950     1400    1.02
%   100     1000    6.4
%   100     2000    8.0
%   950     1000    0.88
%   950     2000    1.2
a =   0.0001696;
b =      0.7738;
c =      -182.5;
d =     -0.5771;
modeltime = (a * numch ^ b) * (numtr + c) + d;
if modeltime > 30
    warning('dg_multiBoot:time', ...
        'Estimated running time for generating this multiBoot is %.0f s', ...
        modeltime);
end

% Do the shufloid:
shufidx = zeros(numtr, numch);
shufidx(:,1) = (1:numtr)';
for ch = 2:numch
    [~, shufidx(:,ch)] = sort(rand(numtr,1));
    for tr = 1:numtr
        if any(shufidx(tr, 1:ch-1) == shufidx(tr,ch))
            isavailtr = true(numtr,1);
            isavailtr(shufidx(tr,1:ch-1)) = false;
            availtr = find(isavailtr);
            shufidx(tr, ch) = availtr(ceil(rand * length(availtr)));
        end
    end
end


