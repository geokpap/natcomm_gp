function [reportstr, reportary] = dg_runClustTests(matfile, retrig, ...
    usepeakvalley)
% Runs a suite of cluster quality tests, returns tab-delimited report
% suitable for tabulation.
%INPUTS
% matfile: *.mat file whose first variable contains per-waveform data with
%   timestamps in seconds in col. 1, cluster IDs in col. 2, and waveform
%   data in the remaining columns.  One row per waveform.
% retrig: value for 'retrig' option to 'dg_shortISI', i.e. the "dead time"
%   after a spike trigger, during which there cannot be another trigger.
% usepeakvalley: optional arg, default is <true>.  If <false>, then only
%   the first 3 principle components are used as features in the call to
%   'dg_rateClusters'.
%OUTPUTS
% reportstr: returns a string containing one line for each cluster in
%   <matfile>, with the following tab-delimited fields:
%       1. Clust Num: from dg_rateClusters.
%       2. numspikes: from dg_rateClusters.
%       3. numwavs: total number of waveforms in file.
%       4. LRatio: from dg_rateClusters; <NaN> indicates inability to get a
%           result from 'mahal'.
%       5. LRatioThresh: from dg_rateClusters; <NaN> indicates inability to
%           get a result from 'mahal'.
%       6. isworthkeeping: rule of thumb based on 'LRatio' and
%           'LRatioThresh' from from dg_rateClustersto identify clusters
%           potentially worth analyzing further.
%       7. SNR: from 'dg_SE_snr'.
%       8. K: from 'dg_shortISI'; <NaN> indicates too few spikes to do the
%           computation; <Inf> indicates bursty MUA with a
%           more-than-Poissonian number of short ISIs.
%       9. minplev: from 'dg_shortISI'; minplev >= 0.05 indicates too few
%           spikes to reach the p<0.05 level for being a single unit (i.e.
%           not MUA).
%       10. avgrate: from <etc> from 'dg_shortISI'; <Inf> or negative value
%           indicates a data error.
%       11. fracshort: <numshort> from 'dg_shortISI' divided by number of
%           points in cluster.
%       12. trunc: logical, <h> returned from 'dg_is_trunc_clust'.
%       13. clustfrac: <frac> returned from 'dg_is_trunc_clust'.
%       14. category: string (possible values: see "Categorization
%           criteria").
% reportary: same content as <reportstr> but formatted as a cell array with
%   one row per cluster.
%NOTES
% LRatio, numspikes, LRatioThresh, and spikeratio as returned from
% 'dg_rateClusters' can be empty iff there is no cluster ID other than "0"
% in the input.  In this case, the empty string <''> is returned as the
% value of <reportstr> without doing any analysis.  Note that this implies
% that <numspikes> cannot be zero for any cluster.
%   The test for warning 'dg_runClustTests:bad' is based on the assumption
% that the spike snapshots are 1 ms at 32 kHz, i.e. 32 samples.  The
% timestamp and cluster number add another 2 columns to make 34.  If you
% use a different number of samples in your spike snapshots, please change
% the number "34" to the right value for your data.

%$Rev: 307 $
%$Date: 2023-10-06 13:06:15 -0400 (Fri, 06 Oct 2023) $
%$Author: dgibson $

if nargin < 3
    usepeakvalley = true;
end
fprintf('dg_runClustTests reading %s\n', matfile);
reportstr = '';
reportary = {};
s = load(matfile);
flds = fieldnames(s);
x = s.(flds{1}); % all multicluster data from the <matfile>
if size(x, 2) ~= 34
    warning('dg_runClustTests:bad', ...
        '%s does not appear to be a cluster file, skipping.', matfile);
    return
end
hasclustID = x(:,2) >= 0;
if ~any(hasclustID)
    warning('dg_runClustTests:empty', ...
        '%s does not contain any clusters, skipping.', matfile);
    return
end
[~, score] = pca(x(:, 3:end));
if usepeakvalley
    pkval = [max(x(:,3:end), [], 2) min(x(:,3:end), [], 2)];
    features = [score(hasclustID, 1:3) pkval(hasclustID, :)];
else
    features = score(hasclustID, 1:3);
end
numwavs = size(features, 1);
[LRatio, numspikes, clustIDs, LRatioThresh] = ...
    dg_rateClusters([x(hasclustID,2) features]);
criterion1 = LRatio < 15 * LRatioThresh;
criterion2 = LRatio < 0.5;
criterion3 = LRatio < 25 * LRatioThresh;
criterion4 = LRatio < 0.2;
criterion5 = LRatio < LRatioThresh;
criterion6 = LRatio < 1;
isworthkeeping = (criterion1 & criterion2) ...
| (criterion3 & criterion4) ...
| (criterion5 & criterion6);
reportary = cell(sum(clustIDs ~= 0), 14);
rownum = 0;
for clustidx = 1:length(clustIDs)
    clustIDnum = clustIDs(clustidx);
    % Do not rate noise cluster:
    if clustIDnum == 0
        continue
    end
    rownum = rownum + 1; % in output and in returns from 'dg_rateClusters'
    [K, ~, minplev, numshort, ~, etc] = ...
        dg_shortISI(x(x(hasclustID,2)==clustIDnum, 1), ...
        'ref', 1.5e-3, 'retrig', retrig);
    % First half: clustID, numspikes, LRatio, LRatioThresh, isworthkeeping,
    % SNR:
    SNR = mean(dg_RMS_snr(x(x(hasclustID,2)==clustIDnum, 3:end), 8));
    reportstr = sprintf( '%s%d\t%d\t%d\t%3.3g\t%3.3g\t%d\t%.2f', ...
        reportstr, clustIDnum, numspikes(clustIDnum), numwavs, ...
        LRatio(clustIDnum), ...
        LRatioThresh(clustIDnum), isworthkeeping(clustIDnum), SNR );
    fracshort = numshort / sum(x(hasclustID,2)==clustIDnum);
    pk_height = max(x(x(hasclustID,2)==clustIDnum, 3:end), [], 2);
    [trunc, clustfrac] = dg_is_trunc_clust(pk_height);
    if ~trunc
        warning('dg_runClustTests:trunc', ...
            'Truncation model does not fit clust %d in %s', ...
            clustidx, matfile);
    end
    % 2nd half: K, minplev, avgrate, fracshort, trunc, clustfrac:
    reportstr = sprintf( '%s\t%.2f\t%.3g\t%.3g\t%.3g\t%d\t%.3g', ...
        reportstr, K, minplev, etc.avgrate, fracshort, trunc, clustfrac );
    cat = dg_catClust( numwavs, numspikes(clustIDnum), LRatio(clustIDnum), ...
        LRatioThresh(clustIDnum), SNR, K, trunc, clustfrac );
    reportstr = sprintf( '%s\t%s\n', reportstr, cat);
    reportary(rownum, :) = { clustIDnum, numspikes(clustIDnum), ...
        numwavs, ...
        LRatio(clustIDnum), LRatioThresh(clustIDnum), ...
        isworthkeeping(clustIDnum), ...
        mean(dg_RMS_snr(x(x(hasclustID,2)==clustIDnum, 3:end), 8)), ...
        K, minplev, etc.avgrate, fracshort, trunc, clustfrac, cat };
end


