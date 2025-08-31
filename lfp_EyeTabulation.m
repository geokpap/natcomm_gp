function [fixations, saccades] = lfp_EyeTabulation(xnum, ynum, vnum, ...
    threshes, varargin)
%LFP_EYETABULATION produces a table of eye movements for enabled trials.
%[fixations, saccades] = lfp_EyeTabulation(xnum, ynum, vnum, threshes)
%lfp_EyeTabulation(xnum, ynum, vnum, threshes)
%lfp_EyeTabulation(xnum, ynum, vnum, threshes, <options>)

% <xnum>, <ynum>, <vnum> are filenums for the x coordinate, y coordinate,
% and velocity traces respectively.  The coordinate traces are assumed to
% be already calibrated, and so observe the "Cartesian" conventions, x<0 =
% left, y<0 = down.  <threshes> is as for lfp_parseEyeTraces.
%
% <fixations> and <saccades> are cell row vectors with one element for each
% selected trial, in the same order as in lfp_TrialIndex.  Each of those
% elements is an array with columns x, y, start time, duration, extra for
% fixations; and direction (in degrees), amplitude, peak velocity, start
% time, duration, extra for saccades.  The column "extra" is 1 if the event
% was considered extra when pairing fixations and saccades, 0 otherwise.
% If saccades are not processed, then "extra" is always 0.  The tables are
% sorted by start time.  Spatial units are the same as they are for the
% traces in <xnum> and <ynum>, and velocity is per sample (so to convert
% from per sample to per second, divide by lfp_SamplePeriod).  Saccades are
% identified strictly by means of SaccStart and SaccEnd events.  Fixations
% are considered to start at either FixStart or threshes.blinkTO s after
% BlinkEnd, and to end at either of FixEnd or threshes.blinkTO s before
% BlinkStart.  For the saccade-fixation pairing, blinks are also treated as
% saccades, but they are removed from the saccades table before further
% processing.
%
% Saccades and fixations are included only if they are contained completely
% within the boundaries of the trial, e.g. FixStart/BlinkEnd is after
% lfp_TrialIndex(trial,1) and FixEnd/BlinkStart is before
% lfp_TrialIndex(trial,2).
%
% If neither of the output arguments is used, lfp_EyeTabulation opens GUIs
% to assign file names to the outputs, and they are saved as files with
% Unique Trial IDs in the first column.
%
%OPTIONS
%   'alldata' - include fixations that are NOT contained
%       completely within the boundaries of the trial, and put NaN for
%       duration (and also for start time for fixations whose ends have
%       been found without their beginnings).  Such fixations
%       are always flagged arbitrarily as NOT being "extra".  The x and y
%       coordinates of the fixation are computed from the portion of the
%       fixation that IS contained withing the trial boundaries.  This
%       option does not directly affect saccades; additional saccades will
%       be included if they are preceded or followed by part of a fixation.
%   'extend' - 
%       WARNING:  THIS FEATURE IS NOT YET FULLY IMPLEMENTED.  The partial
%       implementation shouldn't hurt anything, but invoking it might cause
%       mysterious errors and will certainly not work as specified.  See
%       "lfp_lib maintenance vol 2.doc" for status.
%   'minT', n - is used only when processing saccades, to
%       eliminate trivial false events.  It specifies a minimum time
%       interval in seconds that serves both as the minimum duration of a
%       saccade (otherwise the saccade will be ignored) and as the minimum
%       time between successive fixations (otherwise they will be fused
%       into one long fixation).  'minT' must be followed by a number.
%   'preset', files - unconditionally creates output files, and bypasses
%       GUI for selecting output filenames, using the 2-element cell string
%       array <files> instead; note that if <files> does not contain
%       absolute pathnames, then the output files will be written to the
%       current working directory.  <files>{1} is the fixations table
%       filename, and <files>{2} is saccades.
%   'startwindow', [e1 e2] - same as 'window', except the condition is that
%       the saccade or fixation must START between the two events.
%   'window', [e1 e2] - applies the additional condition that the saccade
%       or fixation must be at least partially contained within a time
%       interval defined by two events.  'window' must be followed by a
%       1x2 cell array as in lfp_disp.  If an eye movement event has
%       exactly the same timestamp as e1 or e2, it is considered to be
%       contained within [e1 e2].
%
    % 10/19/04 DG changed errors to warning...continue in 380-390.
%
% This program uses the following symbolic constants that should be
% assigned in the user's lfp_getEvtIDs_<name> file:
%   BlinkStart
%   BlinkEnd
%   SaccStart
%   SaccEnd
%   FixStart
%   FixEnd

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;
if ~isa(threshes, 'struct')
    error('lfp_EyeTabulation:threshes', ...
        '<threshes> must be a struct as specified for lfp_parseEyeTraces.');
end

argnum = 1;
alldataflag = false;
extendflag = false;
minTflag = false;
presetflag = false;
startwinflag = false;
winflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'alldata'
            alldataflag = true;
        case 'extend'
            extendflag = true;
        case 'minT'
            minTflag = true;
            minT = varargin{argnum+1};
            argnum = argnum + 1;
        case 'preset'
            presetflag = true;
            fixOutputFileName = varargin{argnum + 1}{1};
            saccOutputFileName = varargin{argnum + 1}{2};
            OutputPathName = '';
            argnum = argnum + 1;
        case 'startwindow'
            startwinflag = true;
            window = varargin{argnum+1};
            argnum = argnum + 1;
            if ~(isequal(class(window), 'cell') && isequal(size(window), [1 2]))
                error('lfp_EyeTabulation:badwindow2', ...
                    'window must be a 1x2 cell array' );
            end
        case 'window'
            winflag = true;
            window = varargin{argnum+1};
            argnum = argnum + 1;
            if ~(isequal(class(window), 'cell') && isequal(size(window), [1 2]))
                error('lfp_EyeTabulation:badwindow', ...
                    'window must be a 1x2 cell array' );
            end
        otherwise
            error('lfp_fragmentFiles:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if nargout == 1
    do_saccades = false;
else
    do_saccades = true;
end

if nargout == 0 && ~presetflag
    fixOutputFileName = 'fixations.xls';
    [fixOutputFileName, OutputPathName] = uiputfile( ...
        fullfile(lfp_DataDir, fixOutputFileName), ...
        'Save fixation table as:' );
    if isequal(fixOutputFileName, 0)
        return
    end
    saccOutputFileName = 'saccades.xls';
    [saccOutputFileName, OutputPathName] = uiputfile( ...
        fullfile(lfp_DataDir, saccOutputFileName), ...
        'Save saccade table as:' );
    if isequal(saccOutputFileName, 0)
        return
    end
end
if nargout == 0 || presetflag
    fixoutfid = fopen(fullfile(OutputPathName, fixOutputFileName), 'w');
    if fixoutfid == -1
        error('lfp_EyeTabulation:noOutfid1', ...
            'Could not open fixations output file');
    end
    saccoutfid = fopen(fullfile(OutputPathName, saccOutputFileName), 'w');
    if saccoutfid == -1
        error('lfp_EyeTabulation:noOutfid2', ...
            'Could not open saccades output file');
    end
end

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));

fixations = cell(1, length(trials));
if do_saccades
    saccades = cell(1, length(trials));
end

extraFixCount = 0;
extraSaccCount = 0;
numSkipped = 0;
for trialidx = 1:length(trials)
    fixFrag1 = [];  % for 'alldata'
    fixFrag2 = [];  % for 'alldata'
    fixFrag2Start = []; % for 'alldata'
    trial = trials(trialidx);
    trialevents = lfp_Events( ...
        lfp_TrialIndex(trial,1) : ...
        lfp_TrialIndex(trial,2), : );
    % fixstarts and fixends are indices into trialevents:
    fixstarts = find(ismember(trialevents(:,2), [FixStart BlinkEnd]));
    fixends = find(ismember(trialevents(:,2), [FixEnd BlinkStart]));
    if isempty(fixstarts) || isempty(fixends)
        warning('lfp_EyeTabulation:nofix', ...
            'No fixations in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    
    % Remove unmatched fixation events at beginning and end:
    if trialevents(fixends(1),1) < trialevents(fixstarts(1),1)
        if alldataflag
            % create a dummy record to hold the fixation's coords:
            fixFrag1(3:4) = [NaN NaN];  % start time, duration
            fixFrag1(5) = 0;    % extraFix
            % ...and save the fixation's end time:
            fixFrag1End = trialevents(fixends(1),1);
        end
        fixends(1) = [];
    end
    if isempty(fixstarts) || isempty(fixends)
        warning('lfp_EyeTabulation:nofix2', ...
            'No fixations in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    if trialevents(fixends(end),1) < trialevents(fixstarts(end),1)
        if extendflag
            % search for the matching FixEnd:
            foundEnd = false;
            for evtidx = lfp_TrialIndex(trial,2) + 1 : size(lfp_Events, 1)
                if ismember(lfp_Events(evtidx, 2), [FixEnd BlinkStart])
                    % add it to both <trialevents> and <fixends>:
                    foundEnd = true;
                    fixends(end+1) = size(trialevents, 1) + 1;
                    trialevents(end+1, :) = lfp_Events(evtidx, :);
                    fixdur = lfp_Events(evtidx, 1) - ...
                        trialevents(fixstarts(end),1);
                    if fixdur > 1
                        warning('lfp_EyeTabulation:longExtFix', ...
                            '%d s long fixation found starting at timestamp %d', ...
                            fixdur, trialevents(fixstarts(end),1) );
                    end
                    break
                end
            end
            if ~foundEnd
                warning('lfp_EyeTabulation:noExtFE', ...
                    'Could not find fixation end starting at timestamp %d', ...
                    trialevents(fixstarts(end),1) );
            end
        else
            if alldataflag
                % Save the unmatched FixStart with NaN duration:
                fixFrag2(3) = trialevents(fixstarts(end),1);  % start time
                fixFrag2(4) = NaN;  % duration
                fixFrag2(5) = 0;    % extraFix
                fixFrag2Start = fixstarts(end);
            end
            % delete the unmatched FixStart:
            fixstarts(end) = [];
        end
    end
    if isempty(fixstarts) || isempty(fixends)
        warning('lfp_EyeTabulation:nofix3', ...
            'No fixations in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    if (length(fixstarts) ~= length(fixends)) || ...
            any(fixstarts > fixends)
        warning('lfp_EyeTabulation:mismatchF', ...
            'Fixation fragments in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    % End of pre-processing for fixations
    
    if do_saccades
        saccstarts = find(ismember(trialevents(:,2), [SaccStart BlinkStart]));
        saccends = find(ismember(trialevents(:,2), [SaccEnd BlinkEnd]));
        if isempty(saccstarts) || isempty(saccends)
            warning('lfp_EyeTabulation:nosacc1', ...
                'No saccades in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        % Remove unmatched saccade events at beginning and end:
        if trialevents(saccends(1),1) < trialevents(saccstarts(1),1)
            saccends(1) = [];
        end
        if isempty(saccstarts) || isempty(saccends)
            warning('lfp_EyeTabulation:nosacc4', ...
                'No saccades in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        if trialevents(saccends(end),1) < trialevents(saccstarts(end),1)
            saccstarts(end) = [];
        end
        if isempty(saccstarts) || isempty(saccends)
            warning('lfp_EyeTabulation:nosacc2', ...
                'No saccades in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        % If anything is still ill-formed about saccstarts and saccends,
        % give up on this trial:
        if (length(saccstarts) ~= length(saccends)) || ...
                any(saccstarts>saccends)
            warning('lfp_EyeTabulation:mismatchS', ...
                'Saccade fragments in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        % saccstarts and saccends contain the same numbers of elements from
        % here on.
        
        % Remove saccades that are before the first fixation or after the
        % last one:
        remove1 = [];
        remove2 = [];
        if isempty(fixFrag1)
            firstfixstart = trialevents(fixstarts(1),1);
        else
            firstfixstart = lfp_Events(lfp_TrialIndex(trial,1), 1);
        end
        if trialevents(saccstarts(1),1) < firstfixstart
            remove1 = find( trialevents(saccstarts,1) ...
                < firstfixstart );
        end
        if isempty(fixFrag2)
            if trialevents(saccstarts(end),1) > trialevents(fixstarts(end),1)
                remove2 = find( trialevents(saccstarts,1) ...
                    > trialevents(fixstarts(end),1) );
            end
        end
        remove = union(remove1, remove2);
        saccstarts(remove) = [];
        saccends(remove) = [];
        if isempty(saccstarts) || isempty(saccends)
            warning('lfp_EyeTabulation:nosacc3', ...
                'No saccades in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        
        if minTflag
            % Remove saccades that are too short to be real:
            saccdurations = trialevents(saccends,1) - trialevents(saccstarts,1);
            remove = find(saccdurations < minT);
            saccstarts(remove) = [];
            saccends(remove) = [];
            
            % Fuse fixations that are separated by too little time to be
            % distinct:
            fixspacings = trialevents(fixstarts(2:end),1) ...
                - trialevents(fixends(1:end-1),1);
            remove = find(fixspacings < minT);
            fixends(remove) = [];
            fixstarts(remove+1) = [];
        end

        % Include fixFrag1 and fixFrag2 in a new fixations list if they exist
        % (which will only be when running 'alldata').  The only requirements
        % here are that fixFrag1 sorts before the first saccade, and fixFrag2
        % sorts after the last saccade.  The starting time of fixFrag1 is
        % before the first event in trialevents, so entering its trialevents
        % index as "0" meets the sorting criterion.  fixFrag2Start was
        % previously set along with the other fixFrag2 values.
        if ~isempty(fixFrag1)
            fixFrag1Start = 0;
        else
            fixFrag1Start = [];
        end
        % <fixstartsplus> is a copy of <fixstarts> that includes
        % fragmentary fixations; a fixation tail at the beginning of the
        % trial is represented as 0, so in this case fixstartsplus cannot
        % be used directly as an index into trialevents (unlike fixstarts).
        fixstartsplus = [fixFrag1Start
            fixstarts
            fixFrag2Start];

        % Find & handle missing and extra fixations & saccades; use the
        % starting times to put them in alternating order.  Every saccade
        % will end up either in the pairs table or the extraSacc list,
        % but not both.
        [pairs, extraFix, extraSacc] = dg_zip(fixstartsplus, saccstarts);
        fixdurations = trialevents(fixends,1) - trialevents(fixstarts,1);
        extraFixCount = extraFixCount + length(extraFix);
        extraSaccCount = extraSaccCount + length(extraSacc);
        undefSacc = extraSacc;  % indices into saccstarts
        if length(extraSacc) > 0
            warning('lfp_EyeTabulation:extraSacc', ...
                'Trial %d has extra saccade(s)', trial );
            % Not only the extra saccades but also the preceding saccades
            % have undefined direction and amplitude; remove the preceding
            % saccades from the pairs list and add them to the "undefined" list.
            preceding = extraSacc - 1;
            if preceding(1) == 0
                preceding(1) = [];
            end
            undefSacc(end+1:end+length(preceding)) = preceding;
            rows2delete = find(ismember(pairs(:,2), preceding));
            pairs(rows2delete,:) = [];
        end
        if isempty(extraFix)
            noFinalFix = true;
            extraFixInMiddle = false;
        else
            noFinalFix = (extraFix(end) ~= length(fixstartsplus));
            extraFixInMiddle = (extraFix(1) ~= length(fixstartsplus));
        end
        if noFinalFix
            % If there is no final fixation, the final saccade must also be
            % considered undefined
            undefSacc(end+1) = length(saccstarts);
            if pairs(end,2) == length(saccstarts)
                pairs(end,:) = [];
            end
            warning('lfp_EyeTabulation:noFinalFix', ...
                'Trial %d lacks a final fixation', trial );
        end
        if extraFixInMiddle
            warning('lfp_EyeTabulation:extraFixInMiddle', ...
                'Trial %d has extra fixation(s)', trial );
        end
    end
    % End of pre-processing for saccades
    
    % Generate fixations table.
    
    % Find start and end samples, adjusted to eliminate blink artifact for
    % fixations that start or end with a blink:
    startsamples = lfp_time2index(trialevents(fixstarts,1));
    endsamples = lfp_time2index(trialevents(fixends,1));
    startisblink = trialevents(fixstarts,2) == BlinkEnd;
    endisblink = trialevents(fixends,2) == BlinkStart;
    startsamples(startisblink) = ...
        lfp_time2index( trialevents(fixstarts(startisblink),1) ...
        + threshes.blinkTO );
    endsamples(endisblink) = ...
        lfp_time2index( trialevents(fixends(endisblink),1) ...
        - threshes.blinkTO );

    % Find X coordinates:
    xcoords = zeros(numel(startsamples), 1);
    for idx = 1 : numel(startsamples)
        xcoords(idx) = ...
            mean(lfp_Samples{xnum}(startsamples(idx):endsamples(idx)));
    end
    % Find Y coordinates:
    ycoords = zeros(numel(startsamples), 1);
    for idx = 1 : numel(startsamples)
        ycoords(idx) = ...
            mean(lfp_Samples{ynum}(startsamples(idx):endsamples(idx)));
    end
    % Find X and Y coordinates for fragmentary movements (for 'alldata'):
    if ~isempty(fixFrag1)
        fixEnd1idx = lfp_time2index(fixFrag1End);
        fixFrag1(1) = mean(lfp_Samples{xnum}(...
            lfp_TrialIndex(trial,3):fixEnd1idx ));
        fixFrag1(2) = mean(lfp_Samples{ynum}(...
            lfp_TrialIndex(trial,3):fixEnd1idx ));
    end
    if ~isempty(fixFrag2)
        fixStart2idx = lfp_time2index(fixFrag2(3));
        fixFrag2(1) = mean(lfp_Samples{xnum}(...
            fixStart2idx:lfp_TrialIndex(trial,4) ));
        fixFrag2(2) = mean(lfp_Samples{ynum}(...
            fixStart2idx:lfp_TrialIndex(trial,4) ));
    end
    % Find durations col vector:
    durations = trialevents(fixends,1) - trialevents(fixstarts,1);
    % Save table for this trial; include fixFrag1 and fixFrag2 if they
    % exist (which will only be when running 'alldata').  
    % Note that fixations{trialidx} is here constructed in the same form as
    % fixstartsplus was previously, so an index into fixstartsplus also
    % works for fixations{trialidx}.
    fixations{trialidx} = [ fixFrag1
        [ xcoords ycoords trialevents(fixstarts,1) ...
        durations ismember((1:length(fixstarts))', extraFix) ]
        fixFrag2 ];

        
    if do_saccades
        % Generate saccades table.
        %
        % The <pairs> table (which contains pairs of indices into
        % fixstartsplus, saccstarts) now contains only saccades and blinks
        % with well defined fixations both before and after. pairs(:,1)
        % points at the fixation that precedes the saccade pairs(:,2).
        %
        % The elements of pairs are indices into fixstartsplus and
        % saccstarts, after they may have had elements removed pursuant to
        % minT. The values in fixstarts, fixends, saccstarts, saccends, and
        % fixstartsplus are indices into trialevents, except that
        % fixstartsplus may also contain the value 0.
        
        % Allocate for all "saccades" (based on saccstarts), including
        % blinks: 
        directions = zeros(length(saccstarts),1);
        amplitudes = zeros(length(saccstarts),1);
        velocities = zeros(length(saccstarts),1);
        durations = zeros(length(saccstarts),1);

        pairedsaccidx = pairs(:, 2);  % pairs only!  Indexes saccstarts
        pairedfixidx = pairs(:, 1);   % pairs only!  Indexes fixstartsplus
        
        % Calculate everything possible before removing blinks
        % (otherwise pairedsaccidx and pairedfixidx will no longer be
        % paired)
        deltaX = fixations{trialidx}(pairedfixidx+1, 1) - ...
            fixations{trialidx}(pairedfixidx, 1);
        deltaY = fixations{trialidx}(pairedfixidx+1, 2) - ...
            fixations{trialidx}(pairedfixidx, 2);
        

        % For each well defined saccade, find...
        
        % direction
        s = warning('off', 'MATLAB:divideByZero');
        directions(pairedsaccidx) = 360*atan(deltaY./deltaX)/(2*pi);
        warning(s);
        % put direction in the correct quadrant, using -180 to +180:
        nxpy = find(deltaX<0 & deltaY>=0);
        nxny = find(deltaX<0 & deltaY<0);
        directions(pairedsaccidx(nxpy)) = directions(pairedsaccidx(nxpy)) + 180;
        directions(pairedsaccidx(nxny)) = directions(pairedsaccidx(nxny)) - 180;
        % amplitude
        amplitudes(pairedsaccidx) = sqrt(deltaX.^2 + deltaY.^2);
        % duration
        durations(pairedsaccidx) = trialevents(saccends(pairedsaccidx),1) ...
            - trialevents(saccstarts(pairedsaccidx),1);
        
        % For each undefined saccade, report direction and amplitude as "NaN"
        directions(undefSacc) = NaN;
        % amplitude
        amplitudes(undefSacc) = NaN;
        % peak velocity
        velocities(undefSacc) = lfp_Samples{vnum}( lfp_findmax( ...
            vnum, ...
            lfp_time2index(trialevents(saccstarts(undefSacc),1)), ...
            lfp_time2index(trialevents(saccends(undefSacc),1)) ));
        % duration
        durations(undefSacc) = trialevents(saccends(undefSacc)) ...
            - trialevents(saccstarts(undefSacc));
        
        % Find the entries in saccstarts, saccends that point to
        % blinks rather than true saccades.  In case of disagreement
        % between starts and ends, raise a warning and treat all the
        % disagreeing entries as blinks.  isblinkstart, isblinkend and
        % isblink are all logical arrays that refer to saccstarts and
        % saccends, which are guaranteed to be the same length at this
        % point.
        isblinkstart = trialevents(saccstarts, 2) == BlinkStart;
        isblinkend = trialevents(saccends, 2) == BlinkEnd;
        if ~isequal(isblinkstart, isblinkend)
            badstartidx = isblinkstart & ~isblinkend;
            badendidx = isblinkend & ~isblinkstart;
            if isempty(badstartidx)
                msg1 = '';
            else
                msg1 = sprintf([ 'Putative saccades %s start w/ BlinkStart and end w/ ' ...
                    'SaccEnd;\nstartTS = %s\nendTS = %s\n' ], ...
                    dg_canonicalSeries(badstartidx), ...
                    mat2str(trialevents(saccstarts(badstartidx)), 1), ...
                    mat2str(trialevents(saccends(badstartidx)), 1) );
            end
            if isempty(badendidx)
                msg2 = '';
            else
                msg2 = sprintf([ 'Putative saccades %s start w/ SaccStart and end w/ ' ...
                    'BlinkEnd;\nstartTS = %s\nendTS = %s\n' ], ...
                    dg_canonicalSeries(badendidx), ...
                    mat2str(trialevents(saccstarts(badendidx)), 1), ...
                    mat2str(trialevents(saccends(badendidx)), 1) );
            end
            msg = [ sprintf('Trial %d:\n', trials(trialidx)) msg1 msg2 ...
                'Treating them all as blinks.' ];
            warning('lfp_EyeTabulation:mismatch', msg );
        end
        isblink = (isblinkstart | isblinkend);
        
        % Now remove the blinks from pairedsaccidx:
        ispairedblink = isblink(pairedsaccidx);
        pairedsaccidx(ispairedblink) = [];
       
        % peak velocity must be computed only for paired saccades, and this
        % must be done after removing blinks, because blinks raise
        % misleading 'lfp_findmax:endpoint' warnings:
        velocities(pairedsaccidx) = lfp_Samples{vnum}( lfp_findmax( ...
            vnum, ...
            lfp_time2index(trialevents(saccstarts(pairedsaccidx),1)), ...
            lfp_time2index(trialevents(saccends(pairedsaccidx),1)) ));

        % Construct saccades table for this trial, and remove blinks:
        saccades{trialidx} = [ directions amplitudes velocities ...
            trialevents(saccstarts,1) ...
            durations ismember((1:length(durations))', extraSacc) ];
        saccades{trialidx}(isblink, :) = [];
            
        if winflag || startwinflag
            % Remove fixations that are not contained within
            % the window
            startevts = find(ismember(trialevents(:,2), window{1}));
            if isempty(startevts)
                warning('lfp_EyeTabulation:nostart', ...
                    'There is no start event in trial %d', trial );
                continue
            end
            stopevts = find(ismember(trialevents(:,2), window{2}));
            if isempty(stopevts)
                warning('lfp_EyeTabulation:nostop', ...
                    'There is no stop event in trial %d', trial );
                continue
            end
            startTS = trialevents(startevts(1),1);
            stopTS = trialevents(stopevts(1),1);
            fixstartTS = fixations{trialidx}(:,3);
            fixendTS = fixstartTS + fixations{trialidx}(:,4);
            remove1 = find(fixendTS < startTS);
            if startwinflag
                remove2 = find(fixstartTS > stopTS | fixstartTS < startTS);
            else    % winflag
                remove2 = find(fixstartTS > stopTS);
            end
            fixations{trialidx}([remove1; remove2], :) = [];
            % Remove saccades that are not contained within
            % the window
            saccstartTS = saccades{trialidx}(:,4);
            saccendTS = saccstartTS + saccades{trialidx}(:,5);
            remove1 = find(saccendTS < startTS);
            if startwinflag
                remove2 = find(saccstartTS > stopTS | saccstartTS < startTS);
            else    % winflag
                remove2 = find(saccstartTS > stopTS);
            end
            saccades{trialidx}([remove1; remove2], :) = [];
        end
    end
end

if nargout == 0 || presetflag
    fprintf(fixoutfid, 'Trial\tTrialID\tX\tY\tStart\tDuration\tExtra\n');
    for trialidx = 1:length(trials)
        for row = 1:size(fixations{trialidx},1)
            fprintf(fixoutfid, '%d\t%s\t%d\t%d\t%.6f\t%.6f\t%d\n', ...
                trials(trialidx), lfp_getTrialID(trials(trialidx)), ...
                fixations{trialidx}(row,:) );
        end
    end
    fprintf(fixoutfid, '\n');
    fclose(fixoutfid);
    
    fprintf(saccoutfid, ...
        'Trial\tTrialID\tDirection\tAmplitude\tVelocity\tStart\tDuration\tExtra\n');
    for trialidx = 1:length(trials)
        for row = 1:size(saccades{trialidx},1)
            fprintf(saccoutfid, '%d\t%s\t%d\t%d\t%d\t%.6f\t%.6f\t%d\n', ...
                trials(trialidx), lfp_getTrialID(trials(trialidx)), ...
                saccades{trialidx}(row,:) );
        end
    end
    fprintf(saccoutfid, '\n');
    fclose(saccoutfid);
end

disp(sprintf(...
    '\nTotal extra fixations: %d\nTotal extra saccades: %d', ...
    extraFixCount, extraSaccCount ));
disp(sprintf(...
    'Total number of trials skipped: %d', numSkipped ));
