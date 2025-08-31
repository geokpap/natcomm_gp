function cat = dg_catClust(numwavs, numspikes, LRatio, ...
    LRatioThresh, isworthkeeping, SNR, K, trunc, clustfrac)
% Categorize one cluster based on test results.
%INPUTS
% numwavs: the total number of waveforms in the channel.
% numspikes: number of spikes in cluster. 
% LRatio, LRatioThresh: individual scalar elements of values returned by
%   'dg_rateClusters'.
% isworthkeeping: individual scalar element from <isworthkeeping> as
%   calculated in 'dg_runClustTests' code based on values returned by
%   'dg_rateClusters'.
% SNR: average result of 'dg_RMS_snr' over all waveforms in cluster
%   <clustIDnum>; scalar.
% K: scalar, as returned by 'dg_shortISI'.
% trunc: scalar, as returned by 'dg_is_trunc_clust'.
% clustfrac: scalar, as returned by 'dg_is_trunc_clust'.
%OUTPUT
% cat: a string designating the category to which the cluster was assigned.
%   Possible values:
%     'MUA'
%     'OKunit'
%     'goodunit'
%     'inspect'
%     'poorunit1'
%     'poorunit2'
%     'poorunit3'
%     'poorunit4'
%     'truncated'
%   ...plus any of the above prepended by 'noK' to indicate that there
%   weren't enough spikes to calculate <K>, and/or prepended by 'trunc' to
%   indicate that the cluster was probably truncated by the triggering
%   threshold (except there is no 'noKMUA', but there is 'noKbad').
%NOTES
% The arguments numspikes, LRatio, isworthkeeping, SNR, K, are in the same
% order as the corresponding consecutive columns in the <reportstr>
% returned by 'dg_runClustTests', as are trunc, clustfrac (but there are
% columns missing in between).

%$Rev: 292 $
%$Date: 2022-04-28 17:28:01 -0400 (Thu, 28 Apr 2022) $
%$Author: dgibson $

% Certain conditions don't require complicated analysis to know something's
% fishy:
ratioratio = LRatio / LRatioThresh;
if numspikes / numwavs > 0.9 ...
        || ratioratio > 40 || K < 2 && ratioratio > 30
    cat = 'inspect';
    return
end

% Categorization criteria:
if isnan(LRatio)
    % Could not compute Mahalanobis distances, so there is no <LRatio>
    % value (i.e. no isolation score, and <isworthkeeping> is
    % a nonsense value) and no way to evaluate the ostensible cluster's
    % quality.  If <numspikes> is large, this could indicate a software
    % error, so this raises a warning:
    warning('dg_runClustTests:noLRatio', ...
        'noLRatio, numspikes = %d', numspikes);
    cat = 'noLRatio';
else
    if isworthkeeping
        if isnan(K)
            % Potentially "well isolated", but undecidable as to how
            % many distinct units are in the cluster
            cat = [ 'noK' dg_catClustLRatio(LRatio, ...
                LRatioThresh, numspikes, ...
                numwavs, SNR) ];
        elseif isinf(K)
            % bursty MUA
            cat = 'MUA';
        elseif K>=2
            % not ridiculously non-poissonian, but still MUA
            cat = 'MUA';
        else
            % Cluster boundary should be apparent on visual inspection
            cat = dg_catClustLRatio( LRatio, ...
                LRatioThresh, numspikes, ...
                numwavs, SNR );
        end
        if trunc && clustfrac<0.75
            cat = ['trunc' cat];
        end
    else
        % Not likely to qualify as a true cluster by visual inspection
        if isnan(K)
            % Poorly isolated, and undecidable as to how many distinct
            % units are in the "cluster"
            cat = 'noKbad';
        elseif isinf(K)
            % bursty MUA
            cat = 'MUA';
        elseif K>=2
            % not ridiculously non-poissonian, but still MUA
            cat = 'MUA';
        elseif trunc && clustfrac<0.75
            cat = 'truncpoor';
        else
            % Probably a single unit, but there is no way to know how
            % complete it is because there is no clear cluster boundary
            cat = 'poorunit1';
        end
    end
end
end


function cat = dg_catClustLRatio(LRatio, LRatioThresh, numspikes, ...
    numwaves, SNR)
%INPUTS
% LRatio, LRatioThresh, numspikes, numwaves, SNR: scalar.
%OUTPUT
% cat: string, one of
%   'OKunit'
%   'badunit'
%   'goodunit'
%   'inspect'
%   'poorunit'
spikefrac = numspikes/numwaves;
ratioratio = LRatio/LRatioThresh;
goodiso = LRatio < 0.12 && ratioratio < 4;
OKiso = LRatio < 0.25 - 0.002 * spikefrac ...
    && ratioratio >= 3 && ratioratio < 10 - 5 * spikefrac;
questioniso = LRatio >= 0.05 && ratioratio < 5;
% Putative unit; isolation quality determined from
% LRatio and LRatioThresh
if goodiso
    cat = 'goodunit';
elseif OKiso
    cat = 'OKunit';
elseif questioniso
    if ratioratio > 3.5
        cat = 'poorunit2';
    elseif ratioratio > 1
        cat = 'OKunit';
    elseif spikefrac < 0.01
        % Clusters that contain very few spikes compared to the size of the
        % file tend to receive misleadingly high LRatio values, so we
        % correct that here:
        cat = 'goodunit';
    else
        cat = 'OKunit';
    end
else
    if LRatio < 0.2 && ratioratio < 5
        cat = 'poorunit3';
    elseif SNR > 2 && ( LRatio > 1 && ratioratio > 1 || ...
            SNR * LRatio > 0.1 && SNR * ratioratio > 30 )
        cat = 'inspect';
    else
        cat = 'poorunit4';
    end
end
% ~60% of 'OKunit's turn out to be poorly cut, and in test cases they
% satisfied this condition 5/6 times, whereas 5/5 test cases that didn't
% satisfy the condition didn't need to berecut.  Therefore we change these
% to 'inspect' (although recutting them only changed the category to
% 'goodunit' for 2/6 of them):
if isequal(cat, 'OKunit') && ratioratio > 2.2 && SNR > 3.4
    cat = 'inspect';
end
end
