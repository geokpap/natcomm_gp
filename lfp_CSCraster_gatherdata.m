function [data, interval, trialevents, reftime, window] = ...
    lfp_CSCraster_gatherdata(filenums, trials, window)
global lfp_Samples lfp_XLimAll lfp_TrialIndex lfp_Events lfp_AlignmentRef
global lfp_SamplePeriod
trialevents = [];
reftime = [];

if isempty(lfp_XLimAll) && isempty(window)
    % Use whole trials and pad short ones with NaNs as needed
    interval = [];
    rawtrialinfo = [];
    for trial = reshape(trials, 1, [])
        trialevents = lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), ...
            : ); 
        reftime = trialevents( ...
            ismember(trialevents(:,2), lfp_AlignmentRef), 1 );
        if isempty(reftime)
            refpoint = 0;
        else
            refpoint = lfp_time2index(reftime(1));
        end
        rawtrialinfo(end+1,:) = ...
            [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4) refpoint]; 
    end
else
    if isempty(window)
        window = lfp_XLimAll;
    end
    [interval, rawtrialinfo] = lfp_findCommonTime(trials, 'recseg');
    xlimpoints = round(window/lfp_SamplePeriod);
    interval(1) = max(xlimpoints(1), interval(1));
    interval(2) = min(xlimpoints(2), interval(2));
    if numel(trials) == 1
        trialevents = lfp_Events( ...
            lfp_TrialIndex(trials,1) : lfp_TrialIndex(trials,2), ...
            : ); 
        reftime = trialevents( ...
            ismember(trialevents(:,2), lfp_AlignmentRef), 1 );
        trialevents = lfp_Events( ...
            lfp_Events(:,1) >= reftime + window(1) ...
            & lfp_Events(:,1) <= reftime + window(2), ...
            : );
    end
end
if any(rawtrialinfo(:,3)==0)
    error('lfp_CSCraster:noref', ...
        'The reference event is missing from trial(s) %s', ...
        dg_canonicalSeries(trials(rawtrialinfo(:,3)==0)) );
end
if numel(trials) > 1
    numrows = length(trials);
else
    numrows = length(filenums);
end
if isempty(interval) || numel(filenums) > 1
    % For whole trials or multiple channels, construct the <data> array row
    % by row.
    if isempty(interval)
        ptsBeforeRef = rawtrialinfo(:,3) - rawtrialinfo(:,1);
        ptsAfterRef = rawtrialinfo(:,2) - rawtrialinfo(:,3);
    else
        ptsBeforeRef = min( ...
            rawtrialinfo(:,3) - rawtrialinfo(:,1), -interval(1) );
        ptsAfterRef = min( ...
            rawtrialinfo(:,2) - rawtrialinfo(:,3), interval(2) );
    end        
    numpts = max(ptsBeforeRef) + max(ptsAfterRef) + 1;
    data = NaN(numrows, numpts);
    refpt = max(ptsBeforeRef) + 1;
    startpt = refpt - ptsBeforeRef;
    endpt = refpt + ptsAfterRef;
    if numel(trials) > 1
        for trialidx = 1:length(trials)
            data(trialidx, startpt(trialidx):endpt(trialidx)) = ...
                lfp_Samples{filenums}(rawtrialinfo(trialidx,1) : ...
                rawtrialinfo(trialidx,2)); 
        end
    else
        for filenumidx = 1:length(filenums)
            if isnan(filenums(filenumidx))
                data(filenumidx, startpt:endpt) = NaN;
            else
                data(filenumidx, startpt:endpt) = lfp_Samples{ ...
                    filenums(filenumidx) }( ...
                    rawtrialinfo(3) - ptsBeforeRef : ...
                    rawtrialinfo(3) + ptsAfterRef );
            end
        end
    end
    interval = [-max(ptsBeforeRef) max(ptsAfterRef)];
else
    % This should run faster for large numbers of samples
    idxrange = interval(1):interval(2);
    indices = (repmat(idxrange, size(rawtrialinfo(:,3))) ...
        + repmat(rawtrialinfo(:,3), size(idxrange)) );
    if isempty(indices)
        error('lfp_CSCraster:nodata', ...
            ['No samples were selected; note that if lfp_XLimAll is\n' ...
            'empty, <window> is clipped to start and end of trial.'])
    end
    data = reshape(lfp_Samples{filenums}(indices), size(indices)); 
end
end
