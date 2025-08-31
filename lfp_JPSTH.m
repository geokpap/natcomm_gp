function JPSTHdata = lfp_JPSTH(trials, clustnums, ...
    window, varargin)
% Computes the Joint Peristimulus Time Histogram (JPSTH) according to the
% method of Mati Joshua, Avital Adler, Yifat Prut, Eilon Vaadia, Jeffery R.
% Wickens, and Hagai Bergman, Neuron 62, 695-704, June 11, 2009:  "For this
% analysis, we calculated the raw JPSTH matrix in which the (t1,t2)-th bin
% was the count of the number of times that a coincidence occurred, in
% which neuron #1 spiked in time bin t1 and neuron #2 spiked in time bin t2
% on the same trial (see examples in the first column of Figure S3). To
% correct for rate modulations we calculated the PSTH predictor (Aertsen et
% al., 1989). The predictor matrix is the product of the single-neuron
% PSTHs, i.e., the (t1,t2)-th bin is equal to PSTH1(t1)*PSTH2(t2)... The
% JPSTH was calculated as the subtraction of the number of coincident
% spikes expected by chance (PSTH predictor) from the raw matrix."
%   Time bins are constructed as when using 'newbinwidth' in
% lfp_spikeAnalysis, i.e. in integer arithmetic up to the last possible
% moment. Spike counts are divided by the square of the binwidth in seconds
% times the number of triggers to yield average joint firing rate in
% spikes^2/sec^2.  <window(1)> is a bin edge.  The case where one or both
% neurons emitted more than one spike in a time bin is handled by
% multiplying the numbers of spikes emitted by each neuron and incrementing
% the (t1,t2)-th bin by that product.
%INPUTS
% trials: as for lfp_disp, etc.
% clustnums: a two element vector of cluster numbers, where each cluster
%   represents a putative single neuron.  If empty, defaults to [1 2].
% window: as for lfp_disp, etc.
%OUTPUT
% JPSTHdata: a struct suitable for use with lfp_plotJPSTH.  Contains the
%       following fields, the last batch of which contain the values of the
%       like-named local options variables:
%   matrix - the JPSTH matrix in spikes^2/sec^2, corrected for rate
%       modulation
%   trialcountmatrix - same size as 'matrix', but instead of being
%       incremented by the product of spike counts on each trial, it is
%       only incremented by one.  This is useful to detect situations where
%       a large joint spike rate was produced by a single trial.
%   ntrigs - number of triggers that contributed to the matrix
%   trials - the list of trials actually contributing to the result (i.e.
%       excluding trials that lacked a reference event)
%   trialslabel - as computed when lfp_TrialStyle = 'rule'.
%   timepts - the bin center times for the matrix
%   win - the time window relative to lfp_AlignmentRef that was analyzed
%   spikenames - a two-element cell string array of spike channel names
%   sessionnames - value of lfp_SessionNames
%   align - value of lfp_AlignmentRef
%   evtmatrix - returned by lfp_getSpikes
%   trigfuncArgStr - a string representation of the 'trigfunc' function's
%       arguments (which could be arbitrarily large, and that would be bad)
%   countsflag
%   rawflag
%   evtavg2flag
%   evtbounds
%   multitrigflag
%   norefOKflag
%   trigfuncflag
%   trigfuncH
%OPTIONS
% 'binwidth', N - sets the bin width to <N> milliseconds.  Default is 50.
% 'counts' - does not divide the bin counts by the square of the binwidth
%   times number of triggers.
% 'evtavg', evts2avg - see lfp_getSpikes.
% 'evtavg2', evts2avg - see lfp_getSpikes.
% 'evtbounds', evtbounds - see lfp_getSpikes.
% 'multitrig' - see lfp_getSpikes.
% 'norefOK' - see lfp_getSpikes.
% 'offset' - makes <window(1)> a bin center.
% 'raw' - does not subtract the PSTH predictor matrix.
% 'trigfunc', funcHandle, args - see lfp_getSpikes.
%NOTES
%   The dirty little detail that they gloss over in the paper is how to
% account for the numbers of triggers, to wit:  the expected spikes per bin
% per trig would be PSTH/ntrigs, and so the expected spikes^2 per bin per
% trig between uncorrelated spike trains is PSTH1*PSTH2/ntrigs^2, and the
% expected accumulated number of spikes^2/bin over all trigs is:
%	ntrigs * PSTH1*PSTH2/ntrigs^2 = PSTH1*PSTH2/ntrigs
% i.e. they forgot to mention that factor of ntrigs in computing the
% predictor matrix.
%	"To test whether the population JPSTHs for two different events were
% significantly different, we performed a bin by bin paired t test. The
% surprise values were obtained by transforming the p value of this test by
% ln (p)."  What that means in Matlab terms is:
%       [h,p] = ttest(jpsth1.matrix(:), jpsth2.matrix(:))
% where jpsth1 and jpsth2 are the results returned by lfp_JPSTH.
%EXAMPLES
% Two examples that do the same thing:
%     >> lfp_plotJPSTH(lfp_JPSTH([],[1 1],[-1 2], 'raw', 'offset'), 'smooth', 1);
% 
%     >> thing = lfp_JPSTH([],[1 1],[-1 2], 'raw', 'offset');
%     >> lfp_plotJPSTH(thing, 'smooth', 1);

%$Rev: 324 $
%$Date: 2014-05-12 16:24:57 -0400 (Mon, 12 May 2014) $
%$Author: dgibson $

lfp_declareGlobals;

binwidth = 50;
countsflag = false;
evtavg2flag = false;
evtbounds = {};
multitrigflag = false;
norefOKflag = false;
offsetflag = false;
opts2delete = [];
rawflag = false;
trigfuncArgs = {};
trigfuncflag = false;
trigfuncH = [];

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'binwidth'
                opts2delete(end+1) = argnum; %#ok<*AGROW>
                argnum = argnum + 1;
                binwidth = varargin{argnum};
                if ~isnumeric(binwidth)
                    error('lfp_JPSTH:binwidth', ...
                        '<binwidth> must be numeric.');
                end
                opts2delete(end+1) = argnum; %#ok<*AGROW>
            case 'counts'
                countsflag = true;
                opts2delete(end+1) = argnum; %#ok<*AGROW>
            case 'evtavg'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
                if ~isnumeric(evts2avg)
                    error('lfp_JPSTH:evts2avg', ...
                        '<evts2avg> must be numeric.');
                end
            case 'evtavg2'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
                evtavg2flag = true;
                if ~isnumeric(evts2avg)
                    error('lfp_JPSTH:evts2avg2', ...
                        '<evts2avg> must be numeric.');
                end
            case 'evtbounds'
                argnum = argnum + 1;
                evtbounds = varargin{argnum};
                if ~iscell(evtbounds)
                    error('lfp_JPSTH:evtbounds', ...
                        '<evtbounds> must be a cell array.');
                end
            case 'multitrig'
                multitrigflag = true;
            case 'norefOK'
                norefOKflag = true;
            case 'offset'
                offsetflag = true;
                opts2delete(end+1) = argnum; %#ok<*AGROW>
            case 'raw'
                rawflag = true;
                opts2delete(end+1) = argnum; %#ok<*AGROW>
            case 'trigfunc'
                trigfuncflag = true;
                argnum = argnum +1;
                trigfuncH = varargin{argnum};
                if ~isa(trigfuncH, 'function_handle')
                    error('lfp_JPSTH:trigfuncH', ...
                        '<trigfuncH> must be a function handle.');
                end
                trigfuncArgs = varargin(argnum+1:end);
                argnum = length(varargin);
            otherwise
                error('lfp_JPSTH:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('lfp_JPSTH:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end
varargin(opts2delete) = [];

% lfp_lib boilerplate:
if nargin <3
    window = [];
elseif ~(strcmp(class(window), 'double') ...
        && isequal(size(window), [1 2]) ) ...
        && ~isempty(window)
    error('lfp_JPSTH:badwindow', ...
        '<window> must be 1x2 number array.' );
end
if nargin < 2 || isempty(clustnums)
    clustnums = 1:2;
end
if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
if ischar(trials)
    trials = lfp_parseTrialStr(trials);
end
% <trials> is now numeric, i.e. trialnums.
if any(trials > length(lfp_SelectedTrials))
    warning('lfp_JPSTH:trials', ...
        'Ignoring trials %s, which are beyond the last trial.', ...
        dg_canonicalSeries(trials(trials > length(lfp_SelectedTrials))) );
    trials(trials > length(lfp_SelectedTrials)) = [];
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

if length(clustnums) ~= 2
    error('lfp_JPSTH:clutnums', ...
        '<clustnums> must be a two-element vector.');
end
% Default window:
if isempty(window)
    window = lfp_XLimAll;
end

if isempty(window)
    window = lfp_findCommonTime(trials, 'continuous', 'norefOK');
    if isnan(window(1))
        error('lfp_JPSTH:noref', ...
            'No trials contain the reference event.');
    end
end

if offsetflag
    [binedges, binctrs] = dg_makeTimeBins( ...
        [window(1) - binwidth/2000, window(2)], binwidth);
else
    [binedges, binctrs] = dg_makeTimeBins(window, binwidth);
end

[spikes, triginfo, evtidx, evtmatrix, noreftrialidx] = ...
    lfp_getSpikes(trials, clustnums, window, varargin{:}); %#ok<ASGLU>
if ~isempty(noreftrialidx)
    warning('lfp_JPSTH:noref', ...
        'Trials skipped due to no reference event:\n%s', ...
        dg_canonicalSeries(trials(noreftrialidx)));
end
PSTH = zeros(2, length(binedges) - 1);  % clusts X bins
JPSTH = zeros(length(binedges) - 1, length(binedges) - 1);
trialcount = zeros(length(binedges) - 1, length(binedges) - 1);

% Both spike channels are from the same session, so both have identical
% sets of trials and triggers and both are contained in <spikes>.
for trialidx = 1:length(trials)
    for trignum = 1:size(spikes,2)
        if isempty(spikes{trialidx, trignum, 1}) && ...
                isempty(spikes{trialidx, trignum, 2})
            % no spikes in either channel, so nothing to do
            continue
        end
        % counts1 and counts2 are row vectors of single-trig spike counts
        % per bin for the two spike channels.  We must increment the
        % running totals in PSTH separately for each channel.
        if ~isempty(spikes{trialidx, trignum, 1})
            counts1 = reshape( ...
                histc(spikes{trialidx, trignum, 1}, binedges), 1, [] );
            PSTH(1,:) = PSTH(1,:) + counts1(1:end-1);
        end
        if ~isempty(spikes{trialidx, trignum, 2})
            counts2 = reshape( ...
                histc(spikes{trialidx, trignum, 2}, binedges), 1, [] );
            PSTH(2,:) = PSTH(2,:) + counts2(1:end-1);
        end
        if ~isempty(spikes{trialidx, trignum, 1}) && ...
                ~isempty(spikes{trialidx, trignum, 2})
            % Both channels had spikes, so we must increment the running
            % totals in JPSTH.
            for binidx1 = find(counts1)
                for binidx2 = find(counts2)
                    JPSTH(binidx1, binidx2) = JPSTH(binidx1, binidx2) ...
                        + counts1(binidx1) * counts2(binidx2);
                    trialcount(binidx1, binidx2) = ...
                        trialcount(binidx1, binidx2) + ...
                        (JPSTH(binidx1, binidx2) > 0);
                end
            end
        end
    end
end

predictor = [];
if rawflag
    JPSTHdata.matrix = JPSTH;
else
    predictor = PSTH(1,:)' * PSTH(2,:) / size(triginfo, 1);
    JPSTHdata.matrix = (JPSTH - predictor);
end
if ~countsflag
    JPSTHdata.matrix = JPSTHdata.matrix / ...
        (binwidth^2 * 1e-6 * size(triginfo, 1));
end
JPSTHdata.predictor = predictor;
JPSTHdata.trialcountmatrix = trialcount;
JPSTHdata.ntrigs = size(triginfo, 1);
JPSTHdata.trials = setdiff(trials, trials(noreftrialidx));
JPSTHdata.trialslabel = lfp_getTrialsLabel(JPSTHdata.trials, 'rule');
JPSTHdata.timepts = binctrs;
JPSTHdata.win = window;
JPSTHdata.spikenames = lfp_SpikeNames(clustnums);
JPSTHdata.sessionnames = lfp_SessionNames;
JPSTHdata.align = lfp_AlignmentRef;
JPSTHdata.evtmatrix = evtmatrix;
JPSTHdata.trigfuncArgStr = dg_thing2str(trigfuncArgs);
JPSTHdata.countsflag = countsflag;
JPSTHdata.rawflag = rawflag;
JPSTHdata.evtavg2flag = evtavg2flag;
JPSTHdata.evtbounds = evtbounds;
JPSTHdata.multitrigflag = multitrigflag;
JPSTHdata.norefOKflag = norefOKflag;
JPSTHdata.trigfuncflag = trigfuncflag;
JPSTHdata.trigfuncH = trigfuncH;

