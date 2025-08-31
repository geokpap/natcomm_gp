function [saccades2, enroute] = lfp_saccTarg(xnum, ynum, saccades, ...
    targets, startEvts, stopEvts, RT)
%LFP_SACCTARG adds target ID columns to a saccades table.
%saccades2 = lfp_saccTarg(xnum, ynum, saccades, targets, ...
%       startEvts, stopEvts, RT)
%   lfp_SelectedTrials and lfp_BadTrials must be unchanged since
%   <saccades> was generated.  Unfortunately, all that we can do to
%   enforce that condition is check that there is still the same number of
%   trials selected.

% <xnum>, <ynum> are filenums for the x coordinate and y coordinate
% respectively, which should derive from the same channels that were
% submitted to lfp_EyeTabulation when <saccades> was computed.
% <startEvts>, <stopEvts>, <RT> are as for lfp_TargSeq.
% <saccades> is as returned by lfp_EyeTabulation.
% <saccades2> is a copy of <saccades> with two new columns appended on
% the right containing the target ID from which the saccade started,
% and that on which it ended.  In the case of a trial where no saccades are
% included, saccades2{trial} is empty (this differs from lfp_TargSeq).
% <enroute> is a cell row vector of cell column vectors, each of which
% serves as a third new column for that trial's saccade table; each element
% is the list of targets through which the saccade passed en route.  In
% other words, as the code says, enroute{trialidx}{saccidx} = enroutelist,
% where <enroutelist> is the numeric list of target IDs that were passed
% through.  A target is added to <enroutelist> whenever a new target is
% entered, except that "0" is always ignored for this purpose, and the
% target on which the saccade ended is also not included in <enroutelist>.
% IMPORTANT NOTES:
%   There is no time averaging done here, unlike lfp_TargSeq where the
% target identification is based on the (x, y) coordinates given in the
% fixations table, which are averaged over the duration of the fixation.
% For this reason, it may be desirable to apply additional smoothing to the
% x and y traces before calling this function.  It should also be kept in
% mind that this difference will cause discrepancies (usually small)
% between the saccade start and end positions that you would infer from the
% fixations table and the start and end positions that are used here to
% identify the targets.  Also, since fixation positions are not looked up,
% it is possible here to identify and report "start" and "end" targets for
% "undefined" saccades that have NaN listed as the values for direction and
% amplitude.  That could be a useful feature for some purposes, but for
% other purposes these target identifications may be meaningless, in which
% case you should explicitly check for NaN direction and amplitude before
% using the target ID data.
%   Also, there is no attempt to interpolate between samples.  This means
% that the results returned will match the eye plots, but there could still
% be instances where a target is briefly entered but not reported in
% between two successive time points that were both outside the target.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if length(trials) ~= length(saccades)
    error('lfp_saccTarg:badSelect', ...
        'Trial selection state has changed since generating <saccades>' );
end

startoffset = 0;    % ...for consistency with lfp_TargSeq
stopoffset = 0;     % ditto
saccades2 = saccades;
enroute = cell(size(saccades));
for trialidx = 1:length(saccades)
    if ~isempty(saccades2{trialidx})
        saccades2{trialidx}(:, end+1:end+2) = NaN;
        enroute{trialidx} = cell(size(saccades{trialidx}, 1), 1);
        evtrange = lfp_TrialIndex( ...
            trials(trialidx), 1) : lfp_TrialIndex(trials(trialidx), 2);
        startTS = startoffset + lfp_Events(find(...
            ismember(lfp_Events(evtrange, 2), startEvts) ...
            ) + evtrange(1) - 1, 1);
        stopTS = stopoffset + lfp_Events(find(...
            ismember(lfp_Events(evtrange, 2), stopEvts) ...
            ) + evtrange(1) - 1, 1);
        if isempty(startTS)
            warning('lfp_saccTarg:nostart', ...
                'No start event in trial %d', trials(trialidx) );
            continue
        end
        if isempty(stopTS)
            warning('lfp_saccTarg:nostop', ...
                'No stop event in trial %d', trials(trialidx) );
            continue
        end
        for saccidx = 1:size(saccades{trialidx}, 1)
            saccTS = saccades{trialidx}(saccidx,4);
            saccend = saccTS + saccades{trialidx}(saccidx,5);
            if (saccTS > startTS + RT) && (saccTS < stopTS)
                saccTSidx = lfp_time2index(saccTS);
                saccendidx = lfp_time2index(saccend);
                saccades2{trialidx}(saccidx, end-1) = lfp_xy2targID( ...
                    lfp_Samples{xnum}(saccTSidx), ...
                    lfp_Samples{ynum}(saccTSidx), ...
                    targets);
                saccades2{trialidx}(saccidx, end) = lfp_xy2targID( ...
                    lfp_Samples{xnum}(saccendidx), ...
                    lfp_Samples{ynum}(saccendidx), ...
                    targets);
                enroutelist = [];
                prevsampletarg = saccades2{trialidx}(saccidx, end-1);
                for sampleidx = (saccTSidx + 1) : (saccendidx - 1)
                    targ = lfp_xy2targID( ...
                        lfp_Samples{xnum}(sampleidx), ...
                        lfp_Samples{ynum}(sampleidx), ...
                        targets);
                    if targ && (targ ~= prevsampletarg)
                        enroutelist(end+1) = targ;
                        prevsampletarg = targ;
                    end
                end
                if ~isempty(enroutelist)&& ...
                        enroutelist(end) == saccades2{trialidx}(saccidx, end)
                    if length(enroutelist) == 1
                        % This eliminates unsightly "[1x0 double]" empties:
                        enroutelist = [];
                    else
                        enroutelist(end) = [];
                    end
                end
                enroute{trialidx}{saccidx} = enroutelist;
            elseif saccTS >= stopTS
                % remaining saccades are also after stopTS
                break
            end
        end
    end
end
