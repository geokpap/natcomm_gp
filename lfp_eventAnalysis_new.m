function [missingevts, disordered, intervals, trialIDs] = ...
    lfp_eventAnalysis(sessiondir, eventlist, missingevts, disordered, ...
    repairevts, varargin)
% lfp_read is invoked on <sessiondir> with the events file chosen in the
% following order of preference: 'events.evtsav' 'events.nev' 'events.dat'.
% This makes it easy to perform manual repairs on an events file before
% invoking lfp_eventAnalysis. <eventlist> is a cell array where each entry
% is a row vector of possible alternative event IDs.  
%   <missingevts>, <disordered> are accumulators for counting missing
% events and instances of particular patterns of disordered events.  Both
% are cell arrays containing string labels in first col, counts in second
% col, trial ID list column vectors in third column.
%   <repairevts> is an optional argument which, if equal to 1, saves two 
% different versions of "repaired" events files ("Events_cleaned.evtsav"
% and "Events_intrp.evtsav") for <sessiondir>, where trials that have
% disordered events marked as "bad" in both files (i.e. event 99) and with
% events that are missing replaced by substitutes (subject to the condition
% that there be 4 or fewer instances of the missing event) interpolated in
% between the preceding and the following events at a time such that the
% ratio of the interval to the preceding event and the interval to the
% following event is the same as the ratio between the median preceding
% interval and the median following interval. 
%   When interpolating a missing event that is specified in <eventlist>
% by a list of alternative IDs, the ID that "matches" (i.e.
% occupies the same position in its own list) the last event in the trial
% is used.  If there are more alternatives in the last element of
% <eventlist> than there are for the missing event, the alternatives are
% cycled as many times as needed to find the right match.  E.g., if
% eventlist = {10 [15 16] [17 18 27 28]}, then 15 will be inserted for [15
% 16] in a trial that ends 17 or 27, and 16 will be inserted for [15 16] in
% a trial that ends with 18 or 28.  No attempt is made to replace the first
% or last entries in <eventlist>.  A warning 'lfp_eventAnalysis:hopeless4'
% is raised and no subsitution is made if there is a missing event with
% multiple IDs that is eligible for substitution on a trial that lacks the
% last event.
%   Also, if any trials were marked bad or any substitute events were
% inserted, then all the T files for <sessiondir> are renamed as raw*.t<n>
% (e.g. "rawacq05.t08") and replaced by versions modified to match
% "Events_intrp.evtsav".  
%   If <repairevts> is 2, then the T files are unconditionally re-written
% without modifying the events. This option is for the purpose of
% re-writing all the T files based on new versions of the VT events file
% and (if present) the "events.evtsav" file.
%OPTIONS
% 'delay', eventlistswitches - <eventlistswitches> is a Boolean vector of
%   the same length as <eventlist>.  For each element of
%   <eventlistswitches> that is true, the corresponding event in
%   <eventlist> is interpolated using the median delay time between the
%   preceding event and the missing event instead of the interval ratio
%   described above.  This is more suitable for interpolating e.g. missing
%   Gate events, where the relation to the Warning Click is more reliable
%   than the relation to the following event (Out of Start or Start
%   Locomotion).
% 'substIDs', substIDs - uses the event IDs in <substIDs> instead of the
%   ones in <eventlist> when substituting missing events.

%NOTES
%This is only partially coded in this version (abandoned as unnecessary):
% A log of all event interpolations is kept in the file
% 'lfp_eventAnalysis_log.mat' in the current working directory.  This file
% contains one variable, <lfp_eventAnalysis_log>, which is a cell array
% containing session ID strings in the first column and numeric arrays in
% the second column; those numeric arrays contain event IDs that were
% interpolated in the second column, and the repaired-T-file trial numbers
% where the listed event ID was interpolated in the first column.  If the
% 'lfp_eventAnalysis_log.mat' file does not exist, it is created.

%$Rev: 82 $
%$Date: 2009-09-14 20:27:33 -0400 (Mon, 14 Sep 2009) $
%$Author: dgibson $

eventlistswitches = [];
matchevtsflag = false;
substIDs = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'delay'
            argnum = argnum + 1;
            eventlistswitches = varargin{argnum};
        case 'substIDs'
            argnum = argnum + 1;
            substIDs = varargin{argnum};
        otherwise
            error('lfp_eventAnalysis:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
% There was briefly a 'matchevts' option, but it is now standard, not
% optional: 
matchevtsflag = true;
numtrialtypes = length(eventlist{end});

if isempty(substIDs)
    substIDs = eventlist;
else
    if length(substIDs) ~= length(eventlist)
        error('lfp_eventAnalysis:badsubstIDs', ...
            '<substIDs> must be same length as <eventlist>' );
    end
end

if ~isempty(eventlistswitches) ...
        && length(eventlistswitches) ~= length(eventlist)
    error('lfp_eventAnalysis:eventlistswitches', ...
        '<eventlistswitches> must be same length as <eventlist>' );
end

intervals = [];
trialIDs = {};
evtfilenames = {'Events_intrp.EVTSAV' 'Events_cleaned.EVTSAV' 'Events.EVTSAV' 'events.nev' 'events.dat'};
evtsfound = false;
for k = 1:length(evtfilenames)
    if exist(fullfile(sessiondir, evtfilenames{k}), 'file')
        lfp_read2('preset', sessiondir, evtfilenames(k),'alwaysuseVT','fixeventdup',eventlist,12);
        evtsfound = true;
        break
    end
end
if ~evtsfound
    warning('lfp_eventAnalysis:noevts', ...
        'Session %s contains no events file.', sessiondir);
    return
end
    
lfp_declareGlobals;

% Guarantee that eventlist is a column vector:
eventlist = reshape(eventlist, [], 1);
evtstrings = cell(size(eventlist));
for k = 1:length(eventlist)
    evtstrings{k,1} = mat2str(eventlist{k});
end
if nargin < 3 || isempty(missingevts)
    missingevctr = zeros(numel(eventlist), 1);
    missingevts = cell(numel(eventlist), 3);
    missingevts(:,1) = evtstrings;
    missingevts(:,3) = repmat({{}}, size(missingevts,1), 1);
else
    if ~isequal(missingevts(:,1), evtstrings)
        error('lfp_eventAnalysis:badevtlist', ...
            'The events in <missingevts> do not match those in <eventlist>.');
    end
    missingevctr = cell2mat(missingevts(:,2));
    missingevts(:,2) = cell(size(missingevts(:,2)));
end
if nargin < 4 || isempty(disordered)
    % this is a 2-column association list, key first
    disordered = cell(0, 3);
end
if nargin < 5
    repairevts = 0;
end
if repairevts
    missingevttrials = repmat({[]}, size(eventlist));
end

trials = 1:length(lfp_SelectedTrials);
% Compute inter-event intervals and disordered trials
intervals = NaN(length(trials), length(eventlist)-1);
for trial = 1:length(lfp_SelectedTrials)
    trialevents = lfp_Events( ...
        lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), : );
    thatTS = trialevents(ismember( ...
        trialevents(:,2), eventlist{1} ), 1);
    if isempty(thatTS)
        lfp_log(sprintf('Trial %d has no event %s.', ...
            trial, mat2str(eventlist{1}) ));
        missingevctr(1) = missingevctr(1) + 1;
        missingevts{1,3} = [ missingevts{1,3}
            {lfp_getTrialID(trial)} ];
        thatTS = NaN;
    else
        thatTS = thatTS(1);
    end
    allTS = NaN(length(eventlist),1);
    for evtidx = 1:length(eventlist)
        allTS(evtidx) = thatTS;
        thisTS = thatTS;
        if evtidx < length(eventlist)
            thatTS = ...
                trialevents(ismember(trialevents(:,2), eventlist{evtidx+1}), 1);
            if isempty(thatTS)
                lfp_log(sprintf('Trial %d has no event %s.', ...
                    trial, mat2str(eventlist{evtidx+1}) ));
                missingevctr(evtidx+1) = missingevctr(evtidx+1) + 1;
                missingevts{evtidx+1,3} = [ missingevts{evtidx+1,3}
                    {lfp_getTrialID(trial)} ];
                if repairevts == 1
                    missingevttrials{evtidx+1}(end+1,1) = trial;
                end
                thatTS = NaN;
            else
                thatTS = thatTS(1);
            end
            if ~isempty(thisTS) && ~isempty(thatTS)
                intervals(trial, evtidx) = thatTS - thisTS;
                trialIDs{trial,1} = lfp_getTrialID(trial);
            end
        end
    end
    % NaNs represent missing events; find the order of the remaining
    % events, explicitly marking events that have the same timestamp as
    % disordered:
    remainevts = eventlist(~isnan(allTS));
    [sordid, evtorder] = sort(allTS(~isnan(allTS)));
    consecutive = evtorder(2:end) == evtorder(1:end-1) + 1;
    consecutive = consecutive & ~(sordid(2:end) == sordid(1:end-1));
    if ~any(~consecutive)
        continue
    end
    msg = sprintf('Trial %.0f has disordered events: ', trial);
    for k = 1:length(evtorder)
        msg = sprintf('%s %s', msg, mat2str(remainevts{evtorder(k)}));
    end
    lfp_log(msg);
    singlet = false(size(consecutive));
    singlet(2:end-1) = consecutive(1:end-2) & ~consecutive(2:end-1) & ...
        consecutive(3:end);
    singlet(1) = ~consecutive(1) & consecutive(2);
    singlet(end) = consecutive(end-1) & ~consecutive(end);
    multiplet = ~consecutive & ~singlet;
    % So the events that we need to report as being disordered are the ones
    % pointed to by find(multiplet), plus the ones pointed to by
    % find(singlet) and find(singlet)+1.  For the purpose of computing
    % summary tallies, we will count the number of times each particular
    % run of non-consecutive events occurred, so that means each singlet
    % and each multiplet run.  The key to <disordered> must be a string,
    % because ismember will only work on numbers and strings.
    if any(singlet)
        for evtorderidx = find(singlet')
            key = [ mat2str(remainevts{evtorder(evtorderidx)}) ...
                ' ', mat2str(remainevts{evtorder(evtorderidx+1)}) ];
            keyidx = find(ismember(disordered(:,1), key));
            if isempty(keyidx)
                disordered(end+1,:) = {key 1 {lfp_getTrialID(trial)}};
            else
                disordered{keyidx,2} = disordered{keyidx,2} + 1;
                disordered{keyidx,3} = [ disordered{keyidx,3}
                    {lfp_getTrialID(trial)} ];
            end
        end
    end
    if any(multiplet)
        firstMultipletOfRun = find([ multiplet(1)
            ~multiplet(1:end-1) & multiplet(2:end) ]);
        lastMultipletOfRun = find([ multiplet(1:end-1) & ~multiplet(2:end)
            multiplet(end) ]);
        for multidx = 1:length(firstMultipletOfRun)
            key = sprintf('%s', ...
                    mat2str(remainevts{evtorder( ...
                    firstMultipletOfRun(multidx))}) );
            for offset = 1 : ( lastMultipletOfRun(multidx) - ...
                    firstMultipletOfRun(multidx) )
                key = sprintf('%s %s', key, ...
                    mat2str(remainevts{evtorder( ...
                    firstMultipletOfRun(multidx) + offset )}) );
            end
            keyidx = find(ismember(disordered(:,1), key));
            if isempty(keyidx)
                disordered(end+1,:) = {key 1 {lfp_getTrialID(trial)}};
            else
                disordered{keyidx,2} = disordered{keyidx,2} + 1;
                disordered{keyidx,3} = [ disordered{keyidx,3}
                    {lfp_getTrialID(trial)} ];
            end
        end
    end
end
missingevts(:,2) = mat2cell(missingevctr, ones(size(missingevctr)), 1);

if repairevts
    oldevents = lfp_Events;
    oldtrialindex = lfp_TrialIndex;
    if repairevts == 1
        % Mark disordered trials bad
        badtrialoutcomeidx = [];
        for disidx = 1:size(disordered,1)
            for trialidx = 1:length(disordered{disidx,3})
                trialID = disordered{disidx,3}{trialidx};
                if ~isequal(lfp_SessionNames{1}, ...
                        trialID(1:length(lfp_SessionNames{1})))
                    continue
                end
                trial = lfp_getTrialNum(trialID);
                if trial == size(lfp_TrialIndex, 1)
                    endevt = size(lfp_Events,1);
                else
                    endevt = lfp_TrialIndex(trial+1, 1);
                end
                posttrialevts = lfp_Events(lfp_TrialIndex(trial,2):endevt, :);
                outcomeidx = find(ismember(posttrialevts(:,2), 90:92));
                if length(outcomeidx) ~= 1
                    error('lfp_eventAnalysis:badoutcomeidx', ...
                        'There is not a unique good trial marker in trial %.0f.', ...
                        trial );
                end
                badtrialoutcomeidx = union(badtrialoutcomeidx, ...
                    outcomeidx + lfp_TrialIndex(trial,2) - 1 );
                % for benefit of lfp_save(..., 'rod', 'selectedtrials', ...):
                lfp_SelectedTrials(trial) = false;
            end
        end
        if ~isempty(badtrialoutcomeidx)
            lfp_Events(badtrialoutcomeidx, 2) = 99;
            lfp_save('preset', 'Events_cleaned', 'evt');
        end
        
        % Insert estimates for missing events
        substTS = cell(size(eventlist));
        substType = cell(size(eventlist));
        for evtidx = 2:(length(eventlist)-1)
            if length(missingevttrials{evtidx}) > 4
                msg = sprintf('Too many missing events %s to interpolate', ...
                    evtstrings{evtidx} );
                label = 'lfp_eventAnalysis:hopeless3';
                warning(label, msg);
                lfp_log([label ' ' msg]);
                continue
            end
            premedian = prctile(intervals(:,evtidx-1), 50);
            if isempty(eventlistswitches) || ...
                    ~eventlistswitches(evtidx)
                postmedian = prctile(intervals(:,evtidx), 50);
            else
                postmedian = [];
            end
            if isnan(premedian) || ...
                    (~isempty(postmedian) && isnan(postmedian))
                label = 'lfp_eventAnalysis:hopeless1';
                msg = sprintf('Insufficient data to replace any events %s', ...
                    evtstrings{evtidx} );
                warning(label, msg);
                lfp_log([label ' ' msg]);
                continue
            end
            for trial = missingevttrials{evtidx}'
                trialevents = lfp_Events( ...
                    lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), : );
                preTS = trialevents(ismember( ...
                    trialevents(:,2), eventlist{evtidx-1} ), 1);
                if isempty(eventlistswitches) || ...
                        ~eventlistswitches(evtidx)
                    postTS = trialevents(ismember( ...
                        trialevents(:,2), eventlist{evtidx+1} ), 1);
                else
                    postTS = NaN;
                end
                if isempty(preTS) || isempty(postTS)
                    label = 'lfp_eventAnalysis:hopeless2';
                    msg = sprintf('The missing event %s in trial %.0f cannot be repaired', ...
                        evtstrings{evtidx}, trial );
                    warning(label, msg);
                    lfp_log([label ' ' msg]);
                else
                    if isempty(eventlistswitches) || ...
                            ~eventlistswitches(evtidx)
                        surroundinterval = postTS(1) - preTS(1);
                        substTS{evtidx}(end+1) = preTS(1) + ...
                            (premedian/(premedian+postmedian)) ...
                            * surroundinterval;
                    else
                        substTS{evtidx}(end+1) = preTS(1) + premedian;
                    end
                    if matchevtsflag
                        if length(substIDs{evtidx}) > 1
                            lastevt = find(ismember( ...
                                trialevents(:,2), eventlist{end} ));
                            if isempty(lastevt)
                                label = 'lfp_eventAnalysis:hopeless4';
                                msg = sprintf('The missing event %s in trial %.0f cannot be repaired', ...
                                    evtstrings{evtidx}, trial );
                                warning(label, msg);
                                lfp_log([label ' ' msg]);
                                substType{evtidx}(1, end+1) = NaN;
                            else
                                substType{evtidx}(1, end+1) = mod( ...
                                    find(ismember(eventlist{end}, ...
                                    trialevents(lastevt(1), 2) )) - 1, ...
                                    length(substIDs{evtidx}) ) + 1;
                            end
                        else
                            substType{evtidx}(1, end+1) = 1;
                        end
                    end
                end
            end
        end
        for evtidx = 2:(length(eventlist)-1)
            % This still needs code inserted here to log the missing event
            % insertion to 'lfp_eventAnalysis_log.mat'.
            if matchevtsflag
                for trialtype = 1:length(substIDs{evtidx})
                    lfp_createEvents(@(a) a, ...
                        substIDs{evtidx}(trialtype), ...
                        substTS{evtidx}(substType{evtidx} == trialtype));
                end
            else
                lfp_createEvents(@(a) a, substIDs{evtidx}(1), ...
                    substTS{evtidx});
            end
        end
        if any(~cell2mat(dg_mapfunc(@isempty, substTS)))
            lfp_save('preset', 'Events_intrp', 'evt');
        end
    end
    
    % Make new T files if necessary
    if repairevts == 2 || ~isempty(badtrialoutcomeidx) || ...
            any(~cell2mat(dg_mapfunc(@isempty, substTS)))
        if exist('lfp_eventAnalysis_log.mat')
            load('lfp_eventAnalysis_log.mat');
        else
            lfp_eventAnalysis_log = cell(0,2);
        end
        newevents = lfp_Events;
        newtrialindex = lfp_TrialIndex;
        [animaldir, sessionID] = fileparts(sessiondir);
        files = dir(animaldir);
        tfileidx = zeros(1,0);
        % Check that there is a non-'raw' file for every 'raw' T file
        for fnum = (1:length(files))
            if (~files(fnum).isdir)
                [pathstr,name,ext] = fileparts(files(fnum).name);
                if isequal(['RAW' upper(sessionID)], upper(name)) && ( ...
                        ~isempty(regexpi(ext, '^\.TT[0-9]$')) || ...
                        ~isempty(regexpi(ext, '^\.T[0-9][0-9]$')) )
                    if ~exist(fullfile(animaldir, [upper(sessionID) ext]), 'file')
                        warning('lfp_eventAnalysis:copyingraw', ...
                            'Copying ''raw'' T file to normally-named T file' );
                        copyfile(fullfile(animaldir, files(fnum).name), ...
                            fullfile(animaldir, [upper(sessionID) ext]) );
                    end
                end
            end
        end
        % Rewrite the T files
        for fnum = (1:length(files))
            if (~files(fnum).isdir)
                [pathstr,name,ext] = fileparts(files(fnum).name);
                if isequal(upper(sessionID), upper(name)) && ( ...
                        ~isempty(regexpi(ext, '^\.TT[0-9]$')) || ...
                        ~isempty(regexpi(ext, '^\.T[0-9][0-9]$')) )
                    % Load, move, and re-write the T file
                    srcname = files(fnum).name;
                    newname = sprintf('raw%s%s', name, ext);
                    if exist(fullfile(animaldir, newname), 'file')
                        warning('lfp_eventAnalysis:tfile', ...
                            'The file %s already exists; using it as source for raw data.',...
                            fullfile(animaldir, newname) );
                        copyfile(fullfile(animaldir, newname), ...
                            fullfile(animaldir, srcname), 'f');
                    end
                    startclustidx = length(lfp_Spikes) + 1;
                    lfp_Events = oldevents;
                    lfp_TrialIndex = oldtrialindex;
                    lfp_add('preset', animaldir, {srcname}, ...
                        'Rodent Clusters (*.Tnn, *.TTn)', false );
                    lfp_Events = newevents;
                    lfp_TrialIndex = newtrialindex;
                    if ~exist(fullfile(animaldir, newname), 'file')
                        movefile(fullfile(animaldir, files(fnum).name), ...
                            fullfile(animaldir, newname) );
                    end
                    FileHeader = dg_ReadRodentFormat( ...
                        fullfile(animaldir, newname), 'header' );
                    if ~ismember('LeftCS', fieldnames(FileHeader))
                        if FileHeader.RightCS == 1
                            FileHeader.LeftCS = 8;
                        else
                            FileHeader.RightCS = 1;
                        end
                    end
                     if ~ismember('NoGoCS', fieldnames(FileHeader))
                        FileHeader.NoGoCS = 0;
                     end
                    lfp_save('preset', lfp_SpikeNames(startclustidx:end), ...
                        'rod', 'selectedtrials', 'noclick', ...
                        [ FileHeader.ProcType FileHeader.RightCS ...
                        FileHeader.LeftCS FileHeader.NoGoCS ]);
                end
            end
        end
    end
end




