function cat = dg_catClust(numwavs, numspikes, LRatio, ...
    LRatioThresh, SNR, K, trunc, clustfrac)
% Categorize one cluster based on test results.
%INPUTS
% numwavs: the total number of waveforms in the channel.
% numspikes: number of spikes in cluster. 
% LRatio, LRatioThresh: individual scalar elements of values returned by
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
% The arguments {numspikes, LRatio} are consecutive columns in the
% <reportstr> returned by 'dg_runClustTests', as are {SNR, K}, and also
% {trunc, clustfrac}.

%$Rev: 293 $
%$Date: 2022-06-01 18:51:22 -0400 (Wed, 01 Jun 2022) $
%$Author: dgibson $

noKflag = false; % no estimate of multi/unitary
truncflag = false; % cluster is substantially truncated
inspectflag = false; % cluster may be poorly cut
possflag = false; % it's possible that cluster is actually a 'goodunit'
% Certain conditions don't require complicated analysis to know something's
% fishy:
ratioratio = LRatio / LRatioThresh;
rspkwav = numspikes / numwavs;
if rspkwav > 0.9 ...
        || ratioratio > 40 || K < 2 && ratioratio > 30
    inspectflag = true;
end

% Categorization criteria:
if isnan(LRatio)
    % Could not compute Mahalanobis distances, so there is no <LRatio>
    % value (i.e. no isolation score) and no way to evaluate the ostensible
    % cluster's quality.  If <numspikes> is large, this could indicate a
    % software error, so this raises a warning:
    warning('dg_runClustTests:noLRatio', ...
        'noLRatio, numspikes = %d', numspikes);
    cat = 'noLRatio';
else
    if isnan(K)
        % Potentially "well isolated", but undecidable as to how
        % many distinct units are in the cluster
        noKflag = true;
        cat = dg_catClustLRatio( LRatio, ...
            LRatioThresh, rspkwav );
    elseif isinf(K)
        % bursty MUA
        cat = 'MUA';
    elseif K>=2
        % not ridiculously non-poissonian, but still MUA
        cat = 'MUA';
    else
        % Cluster boundary should be apparent on visual inspection
        cat = dg_catClustLRatio( LRatio, ...
            LRatioThresh, rspkwav );
    end
    if trunc && clustfrac<0.75
        truncflag = true;
    elseif isequal(cat, 'goodiso')
        % good isolation
        if (K < 1.2) && (SNR > 3)
            cat = 'goodunit';
        else
            if (K < 1.2) && (SNR > 2) || (K < 1.4) && (SNR > 3)
                % 'inspect-goodiso':
                inspectflag = true;
            else
                if (K < 1.2) || (SNR > 3)
                    % OR instead of AND
                    cat = 'OKunit';
                else
                    cat = 'poorunit';
                end
            end
        end
    elseif K < 1.2 && SNR > 1.5
        possflag = true;
    end
end
if ~isempty(strfind(cat, 'iso'))
    warning('dg_catClust:unresolvediso', ...
        'cat=%s.', cat);
end
% Check for conditions that are considered cause for visual inspection:
if ismember(cat, {'OKunit' 'OKiso'})
    if ratioratio > 2.2 && SNR > 3.5
        % With that good an SNR, and that bad a <ratioratio>, it's probably
        % mis-cut.
        inspectflag = true;
    else
        % If it just barely failed the LRatio criteria, and has pretty good
        % <K> and <SNR>, we will reclassify it as 'goodunit'.
        if ( (rspkwav < 0.05 && ratioratio < 4) || ...
                (rspkwav >= 0.05 && LRatio < 0.08) ) ...
                && K < 1.1 && SNR > 1.8
            warning('dg_catClust:reclassified', ...
                'Reclassified rspkwav=%d, ratioratio=%d, LRatio=%d, K=%d, SNR=%d\n', ...
                rspkwav, ratioratio, LRatio, K, SNR);
            cat = 'goodunit';
        else
            % At this point it's not obviously mis-cut and not arbitrarily
            % thresholded too tightly in spite of being good on <K> and
            % <SNR>, so it's either 'noKOKiso' or plain old 'OKunit'.
            if ~noKflag
                cat = 'OKunit';
            end
        end
    end
elseif ~isequal(cat, 'goodunit')
    
    inspectflag = true;
end
if noKflag
    cat = ['noK' cat];
end
if possflag
    cat = ['poss-' cat];
end
if truncflag
    cat = ['trunc' cat];
end
if inspectflag
    % some of these can be eliminated from inspection per notes of 05/18/22
    % 07:18:01 PM:
    if ~(LRatio > 0.1 && SNR < 3.2 && clustfrac < 0.92)
        cat = ['inspect-' cat];
    end
end
end


function cat = dg_catClustLRatio(LRatio, LRatioThresh, rspkwav)
%INPUTS
% LRatio, LRatioThresh, rspkwav: scalar.
%OUTPUT
% cat: string.
ratioratio = LRatio/LRatioThresh;
goodiso = (rspkwav < 0.05 && ratioratio < 3) || ...
                (rspkwav >= 0.05 && LRatio < 0.05);
OKiso = LRatio < 0.25 - 0.002 * rspkwav ...
    && ratioratio >= 3 && ratioratio < 10 - 5 * rspkwav;
questioniso = LRatio >= 0.05 && ratioratio < 5;
% Putative unit; isolation quality determined from
% LRatio and LRatioThresh
if goodiso
    cat = 'goodiso';
elseif OKiso
    cat = 'OKiso';
elseif questioniso
    if ratioratio > 3.5
        cat = 'pooriso2';
    elseif ratioratio > 1
        cat = 'OKiso';
    elseif rspkwav < 0.01
        % Clusters that contain very few spikes compared to the size of the
        % file tend to receive misleadingly high LRatio values, so we
        % correct that here:
        cat = 'goodiso';
    else
        cat = 'OKiso';
    end
else
    if LRatio < 0.2 && ratioratio < 5
        cat = 'pooriso3';
    else
        cat = 'MUA';
    end
end
end
