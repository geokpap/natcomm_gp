function [fixations, saccades, enroute] = lfp_EyeTabulation2(xnum, ynum, vnum, ...
    fixEvtIDs, saccEvtIDs, varargin)
%LFP_EYETABULATION2 produces a table of eye movements for enabled trials.
%[fixations, saccades, enroute] = lfp_EyeTabulation2(xnum, ynum, vnum,
%   fixEvtIDs, saccEvtIDs, varargin)

% <xnum>, <ynum>, <vnum> are filenums for the x coordinate, y coordinate,
%   and velocity traces respectively.  The coordinate traces are assumed to
%   be already calibrated, and so observe the "Cartesian" conventions, x<0 =
%   left, y<0 = down.
% <fixEvtIDs> and <saccEvtIDs> are 2 column arrays that contain starts in
%   the first column and ends in the second column.  Fixations and Saccades
%   are determined strictly by the starts and ends given here.
% <fixations> and <saccades> are cell row vectors with one element for each
%   selected trial, in the same order as in lfp_TrialIndex.
%   Fixations columns:
%   1) x (same units as xnum, usually pix)
%   2) y (same units as ynum, usually pix)
%   3) start time (sec)
%   4) duration (sec)
%   5) extra - 1 if the fix start and fix end are not in the time window,
%       otherwise value is 0
%   6) targID (optional)
%   Saccade columns:
%   1) direction (in degrees)
%   2) amplitude (same units as xnum & ynum, usually pix)
%   3) peak velocity (same units as vnum, usually pix/sec)
%   4) start time (sec)
%   5) duration (sec)
%   6) extra - currently nothing written to make this other than 0
%   7) fromID (optional)
%   8) toID (optional)
%   Saccades and fixations are included only if they are contained completely
%       within the boundaries of the trial, e.g. FixStart/BlinkEnd is after
%       lfp_TrialIndex(trial,1) and FixEnd/BlinkStart is before
%       lfp_TrialIndex(trial,2).
%   If neither of the output arguments is used, lfp_EyeTabulation2 opens GUIs
%       to assign file names to the outputs, and they are saved as files with
%       Unique Trial IDs in the first column.
% <enroute> is a cell row vector created only if the 'targIDs' option is
%   used and like fixations and saccades contains one element for each
%   selected trial, in the same order as in lfp_TrialIndex. Each element
%   also contains a cell vector with one element for each saccade in that
%   trial. That secondary cell vetor contains a list, in order, of all the
%   targets the saccade passed through, excluding the first (fromID) and
%   the last (toID) target
% lfp_EyeTabulation2(..., 'alldata') - include fixations that are NOT contained
%   completely within the boundaries of the trial, and put NaN for duration
%   (and also for start time for fixations whose ends have been found
%   without their beginnings).  Such fixations are always flagged as "extra".
%   The x and y coordinates of the fixation are computed from the portion
%   of the fixation that IS contained withing the trial boundaries.  This
%   option does not directly affect saccades.
% lfp_EyeTabulation2(..., 'enddata') - take data from the entire trial but
%   include any fixation that STARTS in the desired period.
% lfp_EyeTabulation2(..., 'extend') - WARNING:  THIS FEATURE IS NOT YET
%   FULLY IMPLEMENTED.  The partial implementation shouldn't hurt anything,
%   but invoking it might cause mysterious errors and will certainly not
%   work as specified.  See "lfp_lib maintenance vol 2.doc" for status.
% lfp_EyeTabulation2(..., 'preset', {files}) - unconditionally creates output
%   files, and bypasses GUI for selecting output filenames, using the
%   2-element cell string array <files> instead; note that if <files> does
%   not contain absolute pathnames, then the output files will be written
%   to the current working directory.  <files>{1} is the fixations table
%   filename, and <files>{2} is saccades.
% lfp_EyeTabulation2(..., 'startwindow', [e1 e2]) - same as 'window',
%   except the condition is that the saccade or fixation must START between
%   the two events.
% lfp_EyeTabulation2(..., 'window', [e1 e2]) - applies the additional
%   condition that the saccade or fixation must be at least partially
%   contained within a time interval defined by two events.  'window' must
%   be followed by a 1x2 cell array as in lfp_disp.  If an eye movement
%   event has exactly the same timestamp as e1 or e2, it is considered to
%   be contained within [e1 e2].
% lfp_EyeTabulation2(..., 'targIDs') - sets the option to add a target ID
%   containing column to the output for <fixations> and fromID and toID
%   columns to the <saccades> table in addition to the <enroute> cell array

% TMD modified from lfp_EyeTabulation 9/9/08
% Last modified 11/21/08 by TMD

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_getEvtIDs;

argnum = 1;
alldataflag = false;
enddataflag = false;
extendflag = false;
minTflag = false;
presetflag = false;
startwinflag = false;
winflag = false;
targIDsflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'alldata'
            alldataflag = true;
        case 'enddata'
            enddataflag = true;
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
                error('lfp_EyeTabulation2:badwindow2', ...
                    'window must be a 1x2 cell array' );
            end
        case 'window'
            winflag = true;
            window = varargin{argnum+1};
            argnum = argnum + 1;
            if ~(isequal(class(window), 'cell') && isequal(size(window), [1 2]))
                error('lfp_EyeTabulation2:badwindow', ...
                    'window must be a 1x2 cell array' );
            end
        case 'targIDs'
            targIDsflag = true;
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

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));

fixations = cell(1, length(trials));
if do_saccades
    saccades = cell(1, length(trials));
    if targIDsflag
        enroute = cell(1, length(trials));
    end
end

taskperiods = lfp_getTaskPeriods;
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
    fixstarts = find(ismember( trialevents(:,2), fixEvtIDs(:,1) ));
    fixends = find(ismember( trialevents(:,2), fixEvtIDs(:,2) ));
    if isempty(fixstarts) || isempty(fixends)
        warning('lfp_EyeTabulation2:nofix', ...
            'No fixations in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    
    % *********************
    % Pre-process fixations
    % *********************
    
    % Remove unmatched fixation events at beginning and end
    % Starts:
    if trialevents(fixends(1),1) < trialevents(fixstarts(1),1)
        if alldataflag
            % create a dummy record to hold the fixation's coords:
            fixFrag1(3:4) = [NaN NaN];  % start time, duration
            fixFrag1(5) = 1;    % extraFix
            extraFixCount = extraFixCount + 1;
            % ...and save the fixation's end time:
            fixFrag1End = trialevents(fixends(1),1);
        end
        fixends(1) = [];
    end
    if isempty(fixstarts) || isempty(fixends)
        warning('lfp_EyeTabulation2:nofix2', ...
            'No fixations in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    % Ends:
    if trialevents(fixends(end),1) < trialevents(fixstarts(end),1)
        if extendflag
            % search for the matching FixEnd:
            foundEnd = false;
            for evtidx = lfp_TrialIndex(trial,2) + 1 : size(lfp_Events, 1)
                if ismember(lfp_Events(evtidx, 2), fixEvtIDs(:,2))
                    % add it to both <trialevents> and <fixends>:
                    foundEnd = true;
                    fixends(end+1) = size(trialevents, 1) + 1;
                    trialevents(end+1, :) = lfp_Events(evtidx, :);
                    fixdur = lfp_Events(evtidx, 1) - ...
                        trialevents(fixstarts(end),1);
                    if fixdur > 1
                        warning('lfp_EyeTabulation2:longExtFix', ...
                            '%d s long fixation found starting at timestamp %d', ...
                            fixdur, trialevents(fixstarts(end),1) );
                    end
                    break
                end
            end
            if ~foundEnd
                warning('lfp_EyeTabulation2:noExtFE', ...
                    'Could not find fixation end starting at timestamp %d', ...
                    trialevents(fixstarts(end),1) );
            end
        else
            if alldataflag || enddataflag
                % Save the unmatched FixStart with NaN duration:
                fixFrag2(3) = trialevents(fixstarts(end),1);  % start time
                fixFrag2(4) = NaN;  % duration
                fixFrag2(5) = 1;    % extraFix
                extraFixCount = extraFixCount + 1;
                fixFrag2Start = fixstarts(end);
            end
            % delete the unmatched FixStart:
            fixstarts(end) = [];
        end
    end
    if isempty(fixstarts) || isempty(fixends)
        warning('lfp_EyeTabulation2:nofix3', ...
            'No fixations in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    if (length(fixstarts) ~= length(fixends)) || ...
            any(fixstarts > fixends)
        warning('lfp_EyeTabulation2:mismatchF', ...
            'Fixation fragments in trial %d, skipping', ...
            trial );
        numSkipped = numSkipped + 1;
        continue
    end
    % End of pre-processing for fixations
    
    % *********************
    % Pre-process saccades
    % *********************    
    
    if do_saccades
        saccstarts = find(ismember( trialevents(:,2), saccEvtIDs(:,1) ));
        saccends = find(ismember( trialevents(:,2), saccEvtIDs(:,2) ));
        if isempty(saccstarts) || isempty(saccends)
            warning('lfp_EyeTabulation2:nosacc1', ...
                'No saccades in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        % Remove unmatched saccade events at beginning and end:
        % Note that there probably shouldn't be any of these because of the
        % way we do parsing now.
        if trialevents(saccends(1),1) < trialevents(saccstarts(1),1)
            saccends(1) = [];
        end
        if isempty(saccstarts) || isempty(saccends)
            warning('lfp_EyeTabulation2:nosacc4', ...
                'No saccades in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        if trialevents(saccends(end),1) < trialevents(saccstarts(end),1)
            saccstarts(end) = [];
        end
        if isempty(saccstarts) || isempty(saccends)
            warning('lfp_EyeTabulation2:nosacc2', ...
                'No saccades in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        % If anything is still ill-formed about saccstarts and saccends,
        % give up on this trial:
        if (length(saccstarts) ~= length(saccends)) || ...
                any(saccstarts>saccends)
            warning('lfp_EyeTabulation2:mismatchS', ...
                'Saccade fragments in trial %d, skipping', ...
                trial );
            numSkipped = numSkipped + 1;
            continue
        end
        % saccstarts and saccends contain the same numbers of elements from
        % here on.
        
        % Housekeeping on Fixations & Saccades
        
        % Include fixFrag1 and fixFrag2 in a new fixations list if they exist
        % The only requirements here are that fixFrag1 sorts before the
        % first saccade, and fixFrag2 sorts after the last saccade.  The
        % starting time of fixFrag1 is before the first event in
        % trialevents, so entering its trialevents index as "0" meets the
        % sorting criterion.  fixFrag2Start was previously set along with
        % the other fixFrag2 values.
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

        % ******************
        % This is where I think the phillosophy of the previous
        % instantiation of the code and my phillosphy now differ.  I don't
        % care if the fixations are paired to saccades...
        % ******************
        
        % Just want all the saccades and fixations that are not truncated
        % in some way to be marked as not extra
        extraFix = zeros(length(fixstartsplus),1);
        extraSacc = zeros(length(saccstarts),1);
    end
    % End of pre-processing for saccades
    
    % ************************
    % Generate fixations table
    % ************************
    
    % Find start and end samples
    startsamples = lfp_time2index(trialevents(fixstarts,1));
    endsamples = lfp_time2index(trialevents(fixends,1));

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
    trialfixations = [ fixFrag1
        [ xcoords ycoords trialevents(fixstarts,1) ...
        durations ismember((1:length(fixstarts))', extraFix) ]
        fixFrag2 ];
    
    % Find the target IDs
    if targIDsflag
        % First get the task period times for this trial (taken from
        % lfp_replaceEyeFixIDs
        % taskperiodtimes contains start (col 1) and end (col 2) timestamps in
        % this trial for each row of taskperiods
        taskperiodtimes = zeros(size(taskperiods,1),2);
        for pd = 1:size(taskperiods,1)
            startidx = find(trialevents(:,2) == taskperiods{pd,2});
            endidx = find(trialevents(:,2) == taskperiods{pd,3});
            if length(startidx) > 1 || length(endidx) > 1
                warning('lfp_replaceEyeFixIDs:multibound', ...
                    'Multiple boundary events for period %s in trial %d', ...
                    taskperiods{pd,1}, trial );
                startidx = startidx(1);
                endidx = endidx(1);
            end
            taskperiodtimes(pd, :) = [ ...
                trialevents(startidx,1) trialevents(endidx,1) ];
        end
        % Cycle through each trial fixation and find the target ID
        % according to the task period
        trialfixations_targID = [trialfixations, zeros(size(trialfixations,1),1)];
        for fixidx = 1:size(trialfixations,1)
            fixTS = trialfixations(fixidx, 3);
            taskperiodidx = find( fixTS >= taskperiodtimes(:,1) ...
                & fixTS < taskperiodtimes(:,2) );
            if isempty(taskperiodidx)
                warning('lfp_replaceEyeFixIDs:noperiod', ...
                    'Fixation at %d is not in any task period; skipping.', fixTS );
                continue
            end
            targID = lfp_xy2targID( ...
                trialfixations(fixidx,1), trialfixations(fixidx,2), ...
                lfp_getTargets(trial, taskperiods{taskperiodidx}) );
            trialfixations_targID(fixidx,6) = targID;
        end
        % Save the table now with Target IDs attached
        fixations{trialidx} = trialfixations_targID;
    else % just save the table
        fixations{trialidx} = trialfixations;
    end
    
    % ***********************
    % Generate Saccades table
    % ***********************
        
    if do_saccades       
        % Allocate for all "saccades" (based on saccstarts)
        directions = zeros(length(saccstarts),1);
        amplitudes = zeros(length(saccstarts),1);
        velocities = zeros(length(saccstarts),1);
        durations = zeros(length(saccstarts),1);
        startXY = zeros(length(saccstarts),2);
        endXY = zeros(length(saccstarts),2);
        if targIDsflag
            starttargs = zeros(length(saccstarts),1);
            endtargs = zeros(length(saccstarts),1);
            enroute{trialidx} = cell(length(saccstarts),1);
        end
        
        % Reminder: saccstarts and saccends are indicies into trialevents
        saccstartsamples = lfp_time2index(trialevents(saccstarts,1));
        saccendsamples = lfp_time2index(trialevents(saccends,1));
        
        % Cycle through all the saccades and find the start and end
        % positions and all targets if applicable according to task period
        for sidx = 1:length(saccstarts)
            % Average over 5 samples to get the start and the end position
            startXY(sidx,1) =...
                mean(lfp_Samples{xnum}(saccstartsamples(sidx):saccstartsamples(sidx)+4));
            startXY(sidx,2) =...
                mean(lfp_Samples{ynum}(saccstartsamples(sidx):saccstartsamples(sidx)+4));
            endXY(sidx,1) =...
                mean(lfp_Samples{xnum}(saccendsamples(sidx)-4:saccendsamples(sidx)));
            endXY(sidx,2) =...
                mean(lfp_Samples{ynum}(saccendsamples(sidx)-4:saccendsamples(sidx)));
            if targIDsflag
                saccTS = trialevents(saccstarts(sidx),1);
                % taskperiodtimes calculated above in fixations table
                taskperiodidx = find( saccTS >= taskperiodtimes(:,1) ...
                    & saccTS < taskperiodtimes(:,2) );
                if isempty(taskperiodidx)
                    warning('lfp_replaceEyeFixIDs:noperiod', ...
                        'Fixation at %d is not in any task period; skipping.', fixTS );
                    continue
                end
                % Important to note that this means the targets that will
                % be used to determine saccade targetIDs will be the
                % targets that are on the screen at saccade start
                targets = lfp_getTargets(trial, taskperiods{taskperiodidx});
                starttargs(sidx) = lfp_xy2targID(startXY(sidx,1), ...
                    startXY(sidx,2), targets);
                endtargs(sidx) = lfp_xy2targID(endXY(sidx,1), ...
                    endXY(sidx,2), targets);
                % Cycle through all of the samples in the current saccade
                % to get the list of targets
                enroutelist = [];
                prevsampletarg = starttargs(sidx);
                % averaged the first five and last five samples to get
                % start and end position so don't need to go through them
                % here
                for sampleidx = saccstartsamples(sidx)+5 : saccendsamples(sidx)-5
                    targ = lfp_xy2targID( ...
                        lfp_Samples{xnum}(sampleidx), ...
                        lfp_Samples{ynum}(sampleidx), ...
                        targets);
                    if targ && (targ ~= prevsampletarg)
                        enroutelist(end+1) = targ;
                        prevsampletarg = targ;
                    end
                end
                % take off the last target in the enroute list if it's the
                % same as the last targ
                if ~isempty(enroutelist)&& ...
                        enroutelist(end) == endtargs(sidx)
                    if length(enroutelist) == 1
                        % This eliminates unsightly "[1x0 double]" empties:
                        enroutelist = [];
                    else
                        enroutelist(end) = [];
                    end
                end
                enroute{trialidx}{sidx} = enroutelist;
            end
        end        

        % For each well defined saccade, find...
        
        % direction
        deltaX = endXY(:,1) - startXY(:,1);
        deltaY = endXY(:,2) - startXY(:,2);
        s = warning('off', 'MATLAB:divideByZero');
        directions = 360*atan(deltaY./deltaX)/(2*pi);
        warning(s);
        % put direction in the correct quadrant, using -180 to +180:
        nxpy = find(deltaX<0 & deltaY>=0);
        nxny = find(deltaX<0 & deltaY<0);
        directions(nxpy) = directions(nxpy) + 180;
        directions(nxny) = directions(nxny) - 180;
        % amplitude
        amplitudes = sqrt(deltaX.^2 + deltaY.^2);
        % duration
        durations = trialevents(saccends, 1) - trialevents(saccstarts, 1);
        % peak velocity
        velocities = lfp_Samples{vnum}( lfp_findmax( ...
            vnum, saccstartsamples, saccendsamples ));
        
        % Construct saccades table for this trial        
        if targIDsflag
            saccades{trialidx} = [ directions amplitudes velocities' ...
                trialevents(saccstarts,1) ...
                durations ismember((1:length(durations))', extraSacc) ...
                starttargs endtargs ];
        else
            saccades{trialidx} = [ directions amplitudes velocities' ...
                trialevents(saccstarts,1) ...
                durations ismember((1:length(durations))', extraSacc) ];
        end
            
        % *************
        % Apply Windows
        % *************

        if winflag || startwinflag
            % Remove fixations that are not contained within
            % the window
            startevts = find(ismember(trialevents(:,2), window{1}));
            if isempty(startevts)
                warning('lfp_EyeTabulation2:nostart', ...
                    'There is no start event in trial %d', trial );
                continue
            end
            stopevts = find(ismember(trialevents(:,2), window{2}));
            if isempty(stopevts)
                warning('lfp_EyeTabulation2:nostop', ...
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
            if targIDsflag
                enroute{trialidx}([remove1; remove2], :) = [];
            end
        end
    end
end

% *****************
% Output text files
% *****************

if nargout == 0 || presetflag
    if targIDsflag
        % [3 4] & [4 5] refer to the columns of fixations that include time
        % stamps as they must be saved with a different precision
        lfp_writeEyeTab(fixOutputFileName, ...
            sprintf('Trial\tTrialID\tX\tY\tStart\tDuration\tExtra\tTargID'), ...
            [3 4], fixations)
        lfp_writeEyeTab(saccOutputFileName, ...
            sprintf('Trial\tTrialID\tDirection\tAmplitude\tVelocity\tStart\tDuration\tExtra\tFromID\tToID\tEnroute'), ...
            [4 5], saccades, enroute)
    else
        lfp_writeEyeTab(fixOutputFileName, ...
            sprintf('Trial\tTrialID\tX\tY\tStart\tDuration\tExtra'), ...
            [3 4], fixations)
        lfp_writeEyeTab(saccOutputFileName, ...
            sprintf('Trial\tTrialID\tDirection\tAmplitude\tVelocity\tStart\tDuration\tExtra'), ...
            [4 5], saccades, enroute)
    end
end
    
disp(sprintf(...
    '\nTotal extra fixations: %d\nTotal extra saccades: %d', ...
    extraFixCount, extraSaccCount ));
disp(sprintf(...
    'Total number of trials skipped: %d', numSkipped ));
