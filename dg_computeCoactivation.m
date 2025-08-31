function [highscore, lowscore, highratio, rawhist, shufhist, rawcoact, ...
    binends, thresh, combinedscore, highratio2] = ...
    dg_computeCoactivation(activations, varargin)
% Computes various representations of coactivation of a set of channels
% over a set of trials, where activation is boolean.  Most measures are
% normalized with respect to the same computation for a number of shuffles
% where no two channels are from the same trial.  Confidence limits for the
% shuffles are computed by doing multiple sets of shuffles and using the
% corresponding prctiles of the shuffles.
%INPUTS
% activations:  a 3D logical array in channels X timepoints X trials form.
%OUTPUTS
% highscore:  in each time bin, the fraction of actual coactivation values
%   that are at least as high as the threshold for the <plevel> tail of one
%   set of shuffles is computed, and divided by the actual size of the
%   <plevel> tail of that set of shuffles.  The computation is repeated
%   separately for each set of shuffles, and then median and CLs are
%   computed. Columns are time bins, row one = median, row two = upper CL,
%   row three = lower CL.
% lowscore:  same computation as <highscore>, but for the low-end tails.
% highratio:  the median value of the high-end <plevel> tail of the actual
%   coactivation counts divided by that for the shuffles.  Columns are
%   time bins, row one = median, row two = upper CL, row three = lower CL.
%   (There is no "lowratio" because that will frequently be 0/0.)
% rawhist:  at each time point, the histogram of the numbers of coactive
%   channels across all trials.  In hist-bins X timepoints form.  The
%   number of hist-bins is equal to the number of channels plus one, with
%   zero active channels as the first bin.
% shufhist:  same as <rawhist>, but for shuffles instead of trials.
%   Columns are timepoints, rows are numbers of simultaneously active
%   channels starting with zero, layer one = median, layer two = upper CL,
%   layer three = lower CL.
% rawcoact:  the number of channels active in each time bin, in trials X
%   bins format.
% binends: the last column in <activations> of each bin.
% thresh:  a row vector of median right-hand tail thresholds for each bin.
% combinedscore:  same computation as <highscore>, but combining both
%   tails.
% highratio2:  same as highratio, except that instead of using the high-end
%   <plevel> tail of the actual coactivation counts, all the coactivation
%   values that are above the threshold of the <plevel> tail of the
%   shuffles is used, i.e. the same threshold value is used to define the
%   tails of both distributions.
%OPTIONS
% 'binsize', binsize - bins results in bins of <binsize> time points,
%   throwing away any timepoints at the end that do not make up a whole
%   bin. Default = 5.
% 'clevel', p - sets the p level for computing confidence limits (so 90%
%   confidence implies <p> = .05, since there are two tails).  Default .05.
% 'numsets', n - sets the number of sets shuffles to <n>.  This is used for
%   calculating confidence limits on the "null hypothesis", i.e. that there
%   is no correlation between different channels.  Default n = 100.  Note
%   that is procedure of doing multiple sets of shuffles only addresses the
%   vagaries introduced by the shuffling process itself.  Generally, the
%   input data will comprise a finite set of observations (e.g. trials)
%   whose finiteness also introduces some vagaries, which are not touched
%   by doing multiple sets of shuffles.  If it is at all possible to
%   bootstrap over the input data set, using only the median values
%   returned by this function, the bootstrap result set will yield more
%   reliable confidence limits.
% 'numshuffles', n - sets the number of shuffles in each set to <n>.
%   Default n = 100.  CLs are computed based on multiple sets of shuffles,
%   but the normalization factors are the median across all sets of
%   shuffles.
% 'offset', tolerance - before starting the analysis, <activations> is
%   converted to mark only the activiation offsets.  Each bin that is
%   <false> following a bin that is <true> is set <true>, along with
%   <tolerance> bins before and after it.  All other bins are set <false>.
%   If 'onset' is combined with 'offset', the resulting <activations> is
%   the logical OR of applying 'onset' and 'offset' separately.
% 'onset', tolerance - same as 'offset', except that each bin that is
%   <true> following a bin that is <false> is set <true> with <tolerance>
%   bins on each side.
% 'plevel', p - sets the p level for computing histogram tails; default
%   p = .05.
% 'rawonly' - computes nothing more than <rawcoact> and <rawhist>.
%NOTES
% All kinds of havoc are wrought by the fact that the number of active
% channels can only have integer values, so it becomes necessary to do
% silly things like compute the actual size of the nominally plevel tail.

%$Rev: 151 $
%$Date: 2012-06-18 17:06:32 -0400 (Mon, 18 Jun 2012) $
%$Author: dgibson $

binsize = 5;  % number of timepoints in one time bin
clevel = .05;
numshuffles = 100;  % number of shuffles in one set pseudotrials
numsets = 100;  % number of sets of pseudotrials (for computing CLs)
offsetflag = false;
onsetflag = false;
plevel = .05;
rawonlyflag = false;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'binsize'
            argnum = argnum + 1;
            binsize = varargin{argnum};
        case 'clevel'
            argnum = argnum + 1;
            clevel = varargin{argnum};
        case 'numsets'
            argnum = argnum + 1;
            numsets = varargin{argnum};
        case 'numshuffles'
            argnum = argnum + 1;
            numshuffles = varargin{argnum};
        case 'offset'
            offsetflag = true;
            argnum = argnum + 1;
            tolerance = varargin{argnum};
        case 'onset'
            onsetflag = true;
            argnum = argnum + 1;
            tolerance = varargin{argnum};
        case 'plevel'
            argnum = argnum + 1;
            plevel = varargin{argnum};
        case 'rawonly'
            rawonlyflag = true;
        case 'verbose'
            verboseflag = true;
        otherwise
            error('dg_computeCoactivation:badoption', ...
                ['The option "' ...
                dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

numchans = size(activations,1);
numtrials = size(activations,3);
if numchans > numtrials
    error('dg_computeCoactivation:shuf1', ...
        'There are not enough trials to do a shuffle.');
else
    s = warning('off', 'MATLAB:nchoosek:LargeCoefficient');
    comboexplosion = nchoosek(numtrials, numchans);
    warning(s.state, 'MATLAB:nchoosek:LargeCoefficient');
    if numshuffles > comboexplosion
        warning('dg_computeCoactivation:shuf2', ...
            'There are not enough trials to do %d unique shuffles.', ...
            numshuffles );
    elseif numshuffles * 1000 > comboexplosion
        warning('dg_computeCoactivation:shuf3', ...
            '<numshuffles> is greater than one-thousandth the possible number of unique shuffles.');
    end
end
if ~islogical(activations)
    warning('dg_computeCoactivation:illogical', ...
        'Converting <activations> to logical.');
    activations = logical(activations);
end

% Do tha binning:
binends = binsize : binsize : size(activations,2);
binned = NaN(size(activations,1), length(binends), size(activations,3));
for binnum = 1:length(binends)
    binned(:, binnum, :) = any( ...
        activations(:, binends(binnum)-binsize+1:binends(binnum), :), 2 );
end
activations = binned;
numbins = size(activations,2);

% Do the onset/offset detection:
if offsetflag || onsetflag
    diffedactivations = false(size(activations));
    if offsetflag
        for trialidx = 1:size(activations,3)
            for chanidx = 1:size(activations,1)
                offsetidx = find( ~activations(chanidx, 2:end, trialidx) ...
                    & activations(chanidx, 1:end-1, trialidx) ) + 1;
                diffedactivations(chanidx, offsetidx, trialidx) = true;
                for tolidx = 1:tolerance
                    diffedactivations(chanidx, ...
                        max(1, offsetidx - tolidx), trialidx) = true;
                    diffedactivations(chanidx, ...
                        min( size(diffedactivations, 2), ...
                        offsetidx + tolidx ), trialidx) = true;
                end
            end
        end
    end
    if onsetflag
        for trialidx = 1:size(activations,3)
            for chanidx = 1:size(activations,1)
                onsetidx = find( activations(chanidx, 2:end, trialidx) ...
                    & ~activations(chanidx, 1:end-1, trialidx) ) + 1;
                diffedactivations(chanidx, onsetidx, trialidx) = true;
                for tolidx = 1:tolerance
                    diffedactivations(chanidx, ...
                        max(1, onsetidx - tolidx), trialidx) = true;
                    diffedactivations(chanidx, ...
                        min( size(diffedactivations, 2), ...
                        onsetidx + tolidx ), trialidx) = true;
                end
            end
        end
    end
    activations = diffedactivations;
end

% Compute raw coactivations scores:
rawcoact = permute(sum(activations, 1), [3 2 1]); % trials X time bins
rawhist = NaN(numchans+2, numbins);
for binidx = 1:numbins
    rawhist(:,binidx) = ...
        histc(rawcoact(:, binidx), -0.5 : (numchans + 0.5));
end
if any(rawhist(end, :))
    error('Burma!');
end
rawhist(end, :) = [];
if rawonlyflag
    highscore = [];
    lowscore = [];
    highratio = [];
    shufhist = [];
    return
end

shuftailthresh1 = NaN(numsets, numbins);
shuftailthresh2 = NaN(numsets, numbins);
shuftailMedianVal = NaN(numsets, numbins);
shuftailsize1 = NaN(numsets, numbins);
shuftailsize2 = NaN(numsets, numbins);
isBadShufTail = false(numbins, numsets);
% <shufhistos> initially has one extra bin because histc returns an
% extra bin for values that are exactly equal to the last binedge
% (always zero here).
shufhistos = NaN(numchans+2, numbins, numsets);
parfor setnum = 1:numsets
    if verboseflag
        fprintf('Starting shuffle set #%d\n', setnum);
    end
    shufcoact = NaN(numshuffles, numbins); % shuffles X time bins
    for shufnum = 1:numshuffles
        % Produce a random shuffle in which no two channels are from the
        % same trial.
        randomlist = rand(numtrials,1);
        [b, randtrialidx] = sort(randomlist); %#ok<*ASGLU>
        while any(diff(sort(randtrialidx)) == 0)
            randomlist = rand(numtrials,1);
            [b, randtrialidx] = sort(randomlist);
        end
        % <randtrialidx> now contains each trial listed once, in random
        % order.  The first <numchans> are just as random as any, so use
        % them as the random selection of different trials for each
        % channel.  Note that the same trial can appear any number of times
        % across different shuffles, both across channels and within
        % channels.
        shuffled = false(numchans, numbins);
        for chan = 1:numchans
            shuffled(chan,:) = activations(chan, :, randtrialidx(chan));  %#ok<PFBNS>
        end
        shufcoact(shufnum,:) = sum(shuffled, 1);
    end
    shufhistos(:, :, setnum) = ...
        histc(shufcoact(:, :), -0.5 : (numchans + 0.5), 1);
    % shuftailthresh<n>, shuftailsize<n> are in shuffle-set X bins format;
    % 1 = left tail, 2 = right tail.  The reason they are no longer
    % shuftailthresh(1,:,:), shuftailthresh(2,:,:), etc. was that it was
    % making the parfor indexing restrictions impossible to satisfy (for no
    % good reason, IMHO, but that's not under my control).
    shuftailthresh1(setnum,:) = prctile(shufcoact, 100*plevel);
    shuftailthresh2(setnum,:) = prctile(shufcoact, 100*(1-plevel));
    shuftailsize2(setnum,:) = sum(shufcoact >= ...
        repmat(shuftailthresh2(setnum,:), numshuffles, 1), 1) ...
        / numshuffles;  %#ok<PFOUS>
    shuftailsize1(setnum,:) = sum(shufcoact <= ...
        repmat(shuftailthresh1(setnum,:), numshuffles, 1), 1) ...
        / numshuffles;  %#ok<PFOUS>
    isBadShufTail(:,setnum) = ...
        shuftailthresh1(setnum,:) == shuftailthresh2(setnum,:);
    % Compute <shuftailMedianVal>, the median value of the right-hand tail
    % of the shuffle. To match <shuftailthresh> etc., <shuftailMedianVal>
    % is in shuffle-set X bins format.  Parallelization trick:
    % make copy of <shufcoact> and set non-tail values to NaN for prctile.
    shufcoacttailvals = shufcoact;
    shufcoacttailvals(shufcoacttailvals < ...
        repmat(shuftailthresh2(setnum,:), size(shufcoact,1), 1)) = NaN;
    shuftailMedianVal(setnum,:) = prctile(shufcoacttailvals, 50, 1);
end
if any(reshape(shufhistos(end, :, :), [], 1))
    error('Kuala Lumpur!');
end
shufhistos(end, :, :) = [];
shuftailsize1(isBadShufTail) = NaN;
shuftailsize2(isBadShufTail) = NaN;

%
% Compute output values
%
if verboseflag
    fprintf('Computing output values\n');
end

% Compute raw tail sizes from the threshold for each shuffle set.  Then
% compute median and CLs across shuffle sets.  <rawtailsize> and <score>
% are in [left-tail|right-tail] X bins X setnum format.  I partially
% de-parallelized this because it was too memory intensive.  rawcoact:
% trials X bins. shuftailthresh<n>, shuftailsize<n>: shuffle-set X bins.
rawtailsize = NaN(2, numbins, numsets);
for setnum = 1:numsets
    for binnum = 1:numbins
        rawtailsize(1,binnum,setnum) = ...
            sum( rawcoact(:, binnum) <= shuftailthresh1(setnum,binnum), ...
            1 ) / numtrials;
        rawtailsize(2,binnum,setnum) = ...
            sum( rawcoact(:, binnum) >= shuftailthresh2(setnum,binnum), ...
            1 ) / numtrials;
    end
end
score = rawtailsize ./ [ permute(shuftailsize1, [3 2 1])
    permute(shuftailsize2, [3 2 1]) ];
lowscore = permute( ...
    prctile(score(1,:,:), [50 100*(1-clevel) 100*clevel], 3), [3 2 1] );
highscore = permute( ...
    prctile(score(2,:,:), [50 100*(1-clevel) 100*clevel], 3), [3 2 1] );

summedscore = sum(rawtailsize,1) ./ ...
    (permute(shuftailsize1, [3 2 1]) + permute(shuftailsize2, [3 2 1]));
combinedscore = permute( prctile(summedscore(1,:,:), ...
    [50 100*(1-clevel) 100*clevel], 3), [3 2 1] );

% Compute ratios of median right tail values.  Parallelization trick: make
% a copy of <rawcoact> and set all the values that are NOT right tail
% values to NaN, so they will have no effect on the "prctile" computation.
rawtailthresh = prctile(rawcoact, 100*(1-plevel)); % row vec time bins
rawcoactTailVals = rawcoact;
rawcoactTailVals(rawcoact < repmat(rawtailthresh, numtrials, 1)) = NaN;
rawtailMedianVal = reshape(prctile(rawcoactTailVals, 50, 1 ), ...
    1, []); % row vector of timebins
ratio = repmat(rawtailMedianVal, [1 1 numsets]) ...
    ./ permute(shuftailMedianVal, [3 2 1]);
highratio = permute( ...
    prctile(ratio, [50 100*(1-clevel) 100*clevel], 3), [3 2 1] );

% Compute ratios of median right tail values as above, except using the
% shuffled tail threshes to compute both the raw and shuffled median tail
% values.
rawtailMedianVal2 = NaN(numsets, numbins);
for setnum = 1:numsets
    rawcoactTailVals2 = rawcoact;
    rawcoactTailVals2( rawcoact < ...
        repmat(shuftailthresh2(setnum,:), numtrials, 1) ) = NaN;
    rawtailMedianVal2(setnum,:) = reshape(prctile(rawcoactTailVals2, 50, 1 ), ...
        1, []); % shuffle-sets X timebins
end
ratio2 = rawtailMedianVal2 ./ shuftailMedianVal;
highratio2 = prctile(ratio2, [50 100*(1-clevel) 100*clevel], 1);


shufhist(:,:,3) = prctile(shufhistos, 100*clevel, 3);
shufhist(:,:,2) = prctile(shufhistos, 100*(1-clevel), 3);
shufhist(:,:,1) = prctile(shufhistos, 50, 3);
thresh = median(shuftailthresh2, 1);

end


