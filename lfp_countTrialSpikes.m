function [totspikes, allgoodtrials] = lfp_countTrialSpikes(...
    aligns, wins, clustnums)
%INPUTS
% aligns: cell array of lfp_AlignmentRef values.
% wins: a cell array of time window values specified in seconds, with one
%   window for each lfp_AlignmentRef value in <aligns>.  Values are only
%   precise to around .001 s.
%OUTPUTS
% totspikes: in trials X clusters format; the total spike counts for each
%   cluster in each trial, summed over all <wins> relative to their
%   respective <aligns>.  Note that overlapping windows will cause some
%   spikes to be counted twice.  Spikes that do not fall in any window are
%   not counted at all.
% allgoodtrials: trialnum for each trial analyzed; same size and in same
%   order as <totspikes>.
%NOTES
% To construct a list of <groups> for Matlab's 'classify' function, run
%   lfp_countTrialSpikes and then do something like:
%       othertrials = setdiff(allgoodtrials, delibtrials);
%       groups = cell(size(allgoodtrials));
%       groups(ismember(allgoodtrials, othertrials)) = {'other'};
%       groups(ismember(allgoodtrials, delibtrials)) = {'delib'};
% Then you can call:
%   [hitrate, pred, grpidx, plevel] = dg_leaveOneOut(totspikes, groups);

%$Rev: 291 $
%$Date: 2012-12-21 13:59:34 -0500 (Fri, 21 Dec 2012) $
%$Author: dgibson $

%lfp_declareGlobals
global lfp_AlignmentRef

lfp_selectByRule('true');
for alignidx = 1:length(aligns)
    lfp_selectByRule( ...
        sprintf('HasEvent(%s)', mat2str(aligns{alignidx})), 'and' );
end
allgoodtrials = lfp_enabledTrials;
binwidth = 1;
totspikes = NaN(length(allgoodtrials), length(clustnums));
for trialidx = 1:length(allgoodtrials)
    trial = allgoodtrials(trialidx);
    for clustidx = 1:length(clustnums)
        clust = clustnums(clustidx);
        for alignidx = 1:length(aligns)
            spikecounts = NaN(size(aligns));
            lfp_AlignmentRef = aligns{alignidx};
            [result, spikes] = lfp_spikeAnalysis('his', trial, clust, ...
                wins{alignidx}, 'newbinwidth', binwidth);
            if length(result) ~= round((wins{alignidx}(2) - wins{alignidx}(1)) ...
                    / (binwidth/1000))
                warning('lfp_countTrialSpikes:roundoff', ...
                    'Wrong # bins returned for trial %d align %s', ...
                    trial, mat2str(lfp_AlignmentRef));
            end
            spikecounts(alignidx) = length(spikes{1});
        end
        totspikes(trialidx, clustidx) = sum(spikecounts);
    end
end

