function lfp_cleanEyeEvents3(minsaccdur, minfixdur, minsaccint, minfixint)
%lfp_cleanEyeEvents3(minsaccdur, minfixdur, minsaccint, minfixint)
%   Edits lfp_Events to eliminate non-paired starts and ends of saccades,
%   and likewise for fixations (the closer starts and ends are kept).
%   Then applies the following algorithm:
%       1.	start with fixations and saccades that meet minimimum duration
%       requirements.
%       2.	for each fixation > mindur, attempt to extend it to the left by
%       fusing it with preceding fixations (of any duration) or saccades
%       (short ones only) that end within minfixint of the current fix
%       start.
%       3.	When you reach the big (i.e. always meets mindur) saccade that
%       started it all, "then you've gone too fah, ya gotta go back!"
%       4.  Delete any remaining saccades that do not meet their mindur.
%       5.  Delete any remaining fixations that do not meet their mindur.
%       6.  Fuse any pairs of saccades that are separated by < minsaccint.
%   minsaccdur:  min saccade duration in seconds
%   minfixdur:  min fixation duration in seconds
%   minsaccint:  min interval between successive saccades in seconds
%   minfixint:  min interval between successive fixations in seconds

% NOTE:  This code depends on these eye event constants being defined in
% the appropriate lfp_getEvtIDs_* script:  SaccStart, SaccEnd, FixStart,
% FixEnd, EyeEvents.

% Modified by DG 10-Apr-2008 from lfp_cleanEyeEvents, lfp_cleanEyeEvents2.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;

% Since they derive from lfp_Events, all vectors are column vectors unless
% specified otherwise:
saccstarts = find(lfp_Events(:,2) == SaccStart);
saccends = find(lfp_Events(:,2) == SaccEnd);

% remove unpaired ("extra") events:
[pairs, extrastart, extraend] = dg_zip(saccstarts, saccends);
disp(sprintf(...
    'Removing %d extra saccade starts, %d extra saccade ends', ...
    length(extrastart), length(extraend) ));
lfp_Events([saccstarts(extrastart); saccends(extraend)],:) = [];
fixstarts = find(lfp_Events(:,2) == FixStart);
fixends = find(lfp_Events(:,2) == FixEnd);
[pairs, extrastart, extraend] = dg_zip(fixstarts, fixends);
disp(sprintf(...
    'Removing %d extra fixation starts, %d extra fixation ends', ...
    length(extrastart), length(extraend) ));
lfp_Events([fixstarts(extrastart); fixends(extraend)],:) = [];
saccstarts = find(lfp_Events(:,2) == SaccStart);
saccends = find(lfp_Events(:,2) == SaccEnd);
fixstarts = find(lfp_Events(:,2) == FixStart);
fixends = find(lfp_Events(:,2) == FixEnd);

% Calculate durations:
if minsaccdur > 0
    saccdurations = lfp_Events(saccends,1) - lfp_Events(saccstarts,1);
    isbigsacc = (saccdurations > minsaccdur);
else
    isbigsacc = true(length(saccstarts));
end
if minfixdur > 0
    fixdurations = lfp_Events(fixends,1) - lfp_Events(fixstarts,1);
    isbigfix = (fixdurations > minsaccdur);
    bigfixidx = find(isbigfix);
else
    bigfixidx = (1:length(fixstarts))';
end

% 2. For each long fixation, going in reverse order, attempt to extend it
% to the left by fusing it with preceding fixations or saccades (of any
% duration) that end within minfixint of the current fix start.
commit2delete = false(size(lfp_Events,1),1); % ... rows from lfp_Events
psaccidx = length(saccends);
disp(sprintf('There are %.0f fixations > minfixdur', length(bigfixidx)));
for k = length(bigfixidx):-1:1
    if psaccidx < 1
        psaccidx = 1;
    end
    fixidx = bigfixidx(k);
    pfixidx = fixidx - 1;
    new2delete = [];
    markedbigsaccidx = 0;
    foundTBSTSIA = false; % "found The Big Saccade That Started It All"
    fusionstart = fixstarts(fixidx);
    fusionstartTS = lfp_Events(fusionstart, 1);
    
    % show progress:
    if ~mod(k,100) 
        disp(sprintf('fix %.0f', k));
    end
    
    % Find the last saccade before fusionstart:
    while psaccidx <= length(saccends) && saccends(psaccidx) < fusionstart
        psaccidx = psaccidx + 1;
    end
    while psaccidx > 0 && saccends(psaccidx) >= fusionstart
        psaccidx = psaccidx - 1;
    end
    
    fusedsacc = true;   % arbitrary, just to start the "while" loop
    while fusedsacc || fusedfix
        pfixendTS = 0;
        psaccendTS = 0;
        fusedfix = false;
        fusedsacc = false;
        if pfixidx ~= 0
            pfixendTS = lfp_Events(fixends(pfixidx),1);
        end
        if psaccidx ~= 0
            psaccendTS = lfp_Events(saccends(psaccidx),1);
        end
        if pfixendTS > psaccendTS
            % preceding eye event is fixation; if it is a long fixation and
            % fusionstart is a long saccade (and therefore is marked as
            % <markedbigsaccidx>), then we've found the Big Saccade That
            % Started It All, so no fusion.
            if markedbigsaccidx && isbigfix(pfixidx) && ...
                    fusionstart == saccstarts(markedbigsaccidx)
                foundTBSTSIA = true;
            end
            if fusionstartTS - pfixendTS < minfixint
                % fuse it
                new2delete(end+1,1) = fusionstart;
                new2delete(end+1,1) = fixends(pfixidx);
                fusionstart = fixstarts(pfixidx);
                fusedfix = true;
                pfixidx = pfixidx - 1;
            end
        else
            % preceding eye event (if any) is saccade
            if psaccidx ~= 0 && (fusionstartTS - psaccendTS < minfixint)
                % fuse it
                if isbigsacc(psaccidx)
                    % Preceding eye event is long sacc; commit deletions up
                    % to but not including the big saccade fusion, and mark
                    % the big saccade as <markedbigsaccidx>.  Then fuse the
                    % saccade tentatively.
                    commit2delete(new2delete) = true;
                    new2delete = [];
                    markedbigsaccidx = psaccidx;
                    new2delete(end+1,1) = fusionstart;
                    new2delete(end+1,1) = saccends(psaccidx);
                    fusionstart = saccstarts(psaccidx);
                    fusedsacc = true;
                    psaccidx = psaccidx - 1;
                else
                    % Preceding eye event is short sacc; if the event that
                    % starts the fusion is a long sacc, DON'T fuse, and
                    % treat fusionstart as Big Saccade That Started It All.
                    % If it is anything else, DO fuse.  If it was a long
                    % sacc, then it should be marked as <markedbigsaccidx>.
                    if lfp_Events(fusionstart, 2) == SaccStart && ...
                            fusionstart == saccstarts(markedbigsaccidx)
                        foundTBSTSIA = true;
                    else
                        new2delete(end+1,1) = fusionstart;
                        new2delete(end+1,1) = saccends(psaccidx);
                        fusionstart = saccstarts(psaccidx);
                        fusedsacc = true;
                        psaccidx = psaccidx - 1;
                    end
                end
            end
        end
        fusionstartTS = lfp_Events(fusionstart, 1);
        if foundTBSTSIA
            % We have found the Big Saccade That Started It All.  Dump
            % those last deletions, mark the next start event as a fixation
            % start, and move on.  The next start event will always be the
            % first item in <new2delete>.
            lfp_Events(new2delete(1), 2) = FixStart;
        end
    end
end
% Now actually execute the deletions and recalculate the events lists:
lfp_Events(commit2delete,:) = [];
saccstarts = find(lfp_Events(:,2) == SaccStart);
saccends = find(lfp_Events(:,2) == SaccEnd);
fixstarts = find(lfp_Events(:,2) == FixStart);
fixends = find(lfp_Events(:,2) == FixEnd);

% 4. Remove remaining saccades that are too short to be real:
if minsaccdur > 0
    saccdurations = lfp_Events(saccends,1) - lfp_Events(saccstarts,1);
    remove = find(saccdurations < minsaccdur);
    disp(sprintf(...
        'Removing %d short saccades', ...
        length(remove) ));
    commit2delete = [saccstarts(remove); saccends(remove)];
end

% 5. Remove remaining fixations that are too short to be real:
if minfixdur > 0
    fixdurations = lfp_Events(fixends,1) - lfp_Events(fixstarts,1);
    remove = find(fixdurations < minfixdur);
    disp(sprintf(...
        'Removing %d short fixations', ...
        length(remove) ));
    commit2delete = [ commit2delete
        fixstarts(remove)
        fixends(remove) ];
end
lfp_Events(commit2delete,:) = [];

% 6. Fuse remaining saccades that are separated by too little time to be
% distinct:
if minsaccint > 0
    saccstarts = find(lfp_Events(:,2) == SaccStart);
    saccends = find(lfp_Events(:,2) == SaccEnd);
    intervals = lfp_Events(saccstarts(2:end),1) ...
        - lfp_Events(saccends(1:end-1),1);
    remove = find(intervals < minsaccint);
    % But don't do it if any other eye events intervene!  (Non-eye events
    % are fine, because they could be absolutely anything.)
    if any(saccstarts(remove+1) - saccends(remove) > 1)
        % Yes, there some intervening events; must examine
        dontremove = [];
        for rmidx = 1:length(remove)
            eventrange = ...
                saccends(remove(rmidx)) + 1 : saccstarts(remove(rmidx)+1) - 1;
            if any(ismember(lfp_Events(eventrange, 2), EyeEvents))
                dontremove(end+1) = rmidx;
            end
        end
        disp(sprintf(...
            'Of %d short inter-saccade intervals, %d were interrupted', ...
            length(remove), length(dontremove) ));
        remove(dontremove) = [];
    end
    disp(sprintf(...
        'Removing %d short inter-saccade intervals', ...
        length(remove) ));
    lfp_Events([saccends(remove); saccstarts(remove+1)],:) = [];
end

lfp_Events = sortrows(lfp_Events);
disp('Recalculating lfp_TrialIndex...');
lfp_TrialIndex = lfp_createTrialIndex(lfp_Events);
[startSampleIndex, endSampleIndex] = ...
    lfp_addTIcols3n4(lfp_TrialIndex(:,1), lfp_TrialIndex(:,2), lfp_Events);
lfp_TrialIndex = [lfp_TrialIndex startSampleIndex endSampleIndex];

lfp_log(sprintf('lfp_cleanEyeEvents3 %d %d %d %d', ...
    minsaccdur, minfixdur, minsaccint, minfixint ));
