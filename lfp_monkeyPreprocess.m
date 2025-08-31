function [params, events] = lfp_monkeyPreprocess(events)

%$Rev: 221 $
%$Date: 2011-04-27 15:32:35 -0400 (Wed, 27 Apr 2011) $
%$Author: dgibson $

lfp_declareGlobals;

% Move contents of Behavioral Data sections from events to
% params, leaving only the BDOutputStart and BDOutputStop events.
% At the same time, search for incomplete BD sections.
%
% The constant caveat here is that data values in the BD section could
% potentially be anything, including the same values as BDOutputStart and
% BDOutputStop.  The occasional caveat is that to allow for unfinished BD
% sections, it is always necessary to check to be sure that the location in
% events that we want to read next actually exists.
%
% In the case of a BD section that is unfinished, it will be detected by
% reading the BD section and failing to find the BDOutputStop at the
% expected location.  Delete all events from and including the previous
% BDOutputStop (which is lfp_NominalTrialStart) to but not including the
% next BDOutputStop.  This will effectively delete the entire trial that
% had the incomplete BD section together with whatever tail end of a trial
% precedes the following BD section, while preserving the correct timing
% for the beginning of the next trial.
%
% For a tail end of a BD section that has no beginning, this will be
% detected by finding a BDOutputStop before finding a BDOutputStart. Delete
% all events from and including the prev BDOutputStop to but not including
% the BDOutputStop of the fragmentary BD section.  In case the unexpected
% BDOutputStop is really a BD value, there will be yet another unexpected
% BDOutputStop; applying the same deletion procedure will delete the BD
% value event that was previously misidentified as the stop signal, without
% requiring any backtracking or other changes to the algorithm.
%
% BD sections that are too short will look like unfinished sections to this
% algorithm, so they will be completely deleted at once, together with
% their trial data and the entire following trial and its BD section.
% Sections that are too long will also look like unfinished sections, but
% in this case, the deletion algorithm results in deletion only of the
% too-long BD section and its own trial data.
%
% The number of parameters in the BD sections, not counting BDOutputStart
% or BDOutputStop, is denoted BDlength in this file. If the first value
% after BDOutputStart is 7 or greater, it is taken as a format version
% number, whose corresponding format is known authoritatively by the
% scientist that uses it.  If it is less than 7, then the data are in the
% "old format", wherein it represents block #, and BDlength is 19.

lfp_getEvtIDs;  % BDOutputStart, BDOutputStop

bdstarts = find(events(:,2) == BDOutputStart);
bdstops = find(events(:,2) == BDOutputStop);
if isempty(bdstops) || isempty(bdstarts)
    error('lfp_monkeyPreprocess:noBD', ...
        'The session is devoid of either BDOutputStarts or BDOutputStops' );
end
if length(bdstarts) == length(bdstops) ...
        && all(bdstops > bdstarts) ...
        && all((bdstops - bdstarts) == (bdstops(1) - bdstarts(1)))
    [params, events] = moveTrialParams(events, bdstarts, bdstops);
    events = cleanupEnd(events, BDOutputStop);
else
    params = {};
    [params, events] = cleanupBrokenTrials(events);
    events = cleanupEnd(events, BDOutputStop);
end
switch lfp_SetupType
    case 'naotaka'
        % Delete any events that have more than one bit set in the high
        % byte plus the high bit of the low byte, as these are almost
        % certainly noise.
        noiseidx = find(events(:,2) > 127);
        numnoise = length(noiseidx);
        numhighbits = zeros(numnoise, 1);
        for k = 7:14
            numhighbits = numhighbits + logical( ...
                bitand(events(noiseidx,2), repmat(2^k, numnoise, 1)) );
        end
        evts2delete = noiseidx(numhighbits>1);
        if ~isempty(evts2delete)
            events(evts2delete, :) = [];
            warning('lfp_monkeyPreprocess:noise', ...
                'Deleted %d noise events', length(evts2delete));
        end
        % Delete all events up to but not including the first event 1, and
        % replace the last BDOutputStop before the first event 1 with dummy
        % event
        starttask = find(events(:,2) == 1);
        if isempty(starttask)
            error('lfp_monkeyPreprocess:nostarttask', ...
                'There are no "start task" events in this session' );
        end
        bdstops = find(events(1:starttask(1),2) == BDOutputStop);
        if isempty(bdstops)
            if starttask(1) > 1
                warning('lfp_monkeyPreprocess:extraevents2', ...
                    'There are events before first task start and no BDOutputStop' );
            end
            warning('lfp_monkeyPreprocess:insertevent', ...
                'Inserting nominal trial start before first event' );
            % Assuming that timestamps are in microseconds here:
            events = [ events(1,1) - 1, 95
                events ];
        else
            if bdstops(end) < (starttask(1) - 1)
                warning('lfp_monkeyPreprocess:extraevents', ...
                    'There are events between first task start and preceding BDOutputStop' );
            end
            if any(events(1:bdstops(end), 2) == 96)
                % There was an entire BD section before the first trial.
                % That means we need to delete it from <params>.
                params(1) = [];
            end
            events(bdstops(end), 2) = 95;
            events(1:(bdstops(end)-1), :) = [];
        end

        % Delete all events starting with the event after the last bdstop
        % through the end of the file
        bdstops = find(events(:,2) == BDOutputStop);
        events(bdstops(end)+1:end, :) = [];
        
        % unpack the BD parameters; this is Naotaka-only, so BDlength = 5.
        % The unpacked parameters are stored back in <params> in this
        % order:
        %   1-6. saccade directions #1 through #6
        %   7. yellow
        %   8. target off timing
        %   9. interval
        %   10. task condition
        %   11. number of saccades
        %   12. grid on
        for trial = 1:length(params)
            unpacked = zeros(1,12);
            sacc12data = params{trial}(1);
            if sacc12data < 98
                unpacked([1 2 7]) = -1;
            else
                % the '-2' and '-98' both represent the offset of 0x62 that
                % Naotaka added to everything, but the bitshifts actually
                % make the upper or lower digit irrelevant depending on
                % which half-byte you're unpacking
                unpacked(1:2) = getSaccadeCode(bitand(sacc12data-2,15));
                unpacked(7) = bitshift(sacc12data-98,-4);
            end
            
            sacc34data = params{trial}(2);
            if sacc34data < 2
                unpacked([3 4 8]) = -1;
            else
                unpacked(3:4) = getSaccadeCode(bitand(sacc34data-2,15));
                targetOffCode = bitshift(sacc34data,-4);
                if targetOffCode == 0 | targetOffCode == 1
                    unpacked(8) = targetOffCode + 2;
                elseif targetOffCode == 6 | targetOffCode == 7
                    unpacked(8) = targetOffCode - 6;
                else
                    unpacked(8) = -1; % something is wrong.
                end
            end
            
            sacc56data = params{trial}(3);
            if sacc56data < 2
                unpacked([5 6 9]) = -1;
            else
                unpacked(5:6) = getSaccadeCode(bitand(sacc56data-2,15));
                interval0 = bitshift(sacc56data,-4);
                if  interval0 == 6		% interval 0
                    unpacked(9) = 4444;
                elseif interval0 == 7	% interval 1
                    unpacked(9) = 4844;
                elseif interval0 == 0	% interval 2
                    unpacked(9) = 8444;
                elseif interval0 == 1	% interval 3
                    unpacked(9) = 6666;
                elseif interval0 == 2	% interval 4
                    unpacked(9) = 8888;
                else
                    unpacked(9) = -1;   % error
                end
            end
            
            unpacked(10) = params{trial}(4) - hex2dec('62');
            
            numsaccdata = params{trial}(5) ...
                - hex2dec('62') - params{trial}(4) + hex2dec('80');
            if numsaccdata < 0
                unpacked([11 12])= -1;
            else
                unpacked(11) = bitand(numsaccdata, 15);
                unpacked(12) = bitshift(numsaccdata, -4);
            end
            params{trial} = unpacked;
        end
end



function events = cleanupEnd(events, stopcode)
% Delete any events after the end of the last BD section
lfp_declareGlobals;
bdstops = find(events(:,2) == stopcode);
events(bdstops(end)+1:end, :) = [];



function [params, newEvents] = moveTrialParams(events, bdstarts, bdstops)
% The quick and simple way, does not correct errors or accommodate to
% variable BD section lengths.  bdstarts and bdstops must be of equal
% length to use this function.  To avoid excessive copying of events,
% the edited version is constructed in newEvents and then returned.
lfp_declareGlobals;
BDlength = bdstops(1) - bdstarts(1) - 1;
hWaitBar = waitbar(0, '', 'Name', 'Moving Simple BD Sections');
newEvents = zeros(length(events) - length(bdstarts)*BDlength, 2);
new_lastrow = 0;    % last row written in newEvents
old_lastrow = 0;    % last row read from events
for trial = 1:length(bdstarts)
    waitbar(trial/length(bdstarts), hWaitBar);
    params{trial} = ...
        events(bdstarts(trial)+1 : bdstops(trial)-1, 2)';
    new_lastrow2 = new_lastrow + (bdstarts(trial) - old_lastrow);
    newEvents(new_lastrow + 1 : new_lastrow2, :) ...
        = events(old_lastrow + 1 : bdstarts(trial), :);
    new_lastrow = new_lastrow2;
    old_lastrow = bdstops(trial) - 1;
end
newEvents(new_lastrow + 1 : end, :) ...
    = events(old_lastrow + 1 : end, :);
close(hWaitBar);


function [params, events] = moveOneTrial(params, events, trial, start, stop)
% Moves trial params for one trial
lfp_declareGlobals;
params{trial} = ...
    events(start+1 : stop-1, 2)';
events(start+1 : stop-1, :) = [];


function [params, events] = cleanupBrokenTrials(events)
% The slow and complicated way, removes fragmentary BD sections
lfp_declareGlobals;
BDOutputStart = 96;
BDOutputStop = 97;
bdstarts = find(events(:,2) == BDOutputStart);
bdstops = find(events(:,2) == BDOutputStop);
lastgood = [];   % points to end of last certified complete trial
trial = 1;
params = {};

% Eliminate initial tail ends with no beginnings.
% If the very first event is a bdstart, treat the dangling BD section as a
% tail fragment.
while bdstops(1) < bdstarts(1) || bdstarts(1) == 1
    warning('lfp_monkeyPreprocess:deletetail1', ...
        'Deleting initial tail fragment around timestamp %g.', ...
        events(bdstops(1), 1) );
    [startdel, events] = deleteTrialFragments(events, bdstops(1));
    bdstarts = find(events(:,2) == BDOutputStart);
    bdstops = find(events(:,2) == BDOutputStop);
    if isempty(bdstops) || isempty(bdstarts)
        error('lfp_monkeyPreprocess:noBD2', ...
            'The session has no intact BD sections' );
    end
end

hWaitBar = waitbar(0, '', 'Name', 'Processing Variable BD Sections');

while ~isempty(bdstarts)
    waitbar(bdstarts(1)/length(events), hWaitBar);
    % Get format code
    formatnum = 0;
    if bdstarts(1) < size(events, 1)
        formatmarker = events(bdstarts(1)+1, 2);
        formatmarker = bitand(formatmarker, 255);
    else
        formatmarker = 0;
    end
    if strcmp(lfp_SetupType, 'naotaka')    % Naotaka format
        BDlength = 5;
    elseif formatmarker < 7 % old Joey format
        BDlength = 19;
    elseif formatmarker == 7    % Joey format 7
        BDlength = 30;
    else
        % 2-byte format code.
        % Construct format code from this and the next param, treating
        % them respectively as the high-order and low-order pairs of hex
        % digits of a single number.
        %   NOTE: this entire function does not run if the conditions for
        %   running moveTrialParams instead are met; consequently, this is
        %   NOT a good place to check for unknown formats.
        if bdstarts(1) < size(events, 1) + 1
            formatmarker2 = events(bdstarts(1)+2, 2);
            formatmarker2 = bitand(formatmarker2, 255);
        else
            formatmarker2 = 0;
        end
        formatnum = 256*formatmarker + formatmarker2;
        switch formatnum
            case 2048   % 0x800: start of Joey's new (Feb. 2004) formats
                BDlength = 30;
            case 2304   % 0x900: start of Theresa's formats
                BDlength = 28;
            case 2305   % 0x901
                BDlength = 32;
            case 2306   % 0x902
                BDlength = 33;
            case 2307   % 0x903
                BDlength = 36;
            case 2308   % 0x904
                BDlength = 36;
            case 2309   % 0x905
                BDlength = 61;
            case 2310   % 0x906
                BDlength = 64;
            case 2311   % 0x907
                BDlength = 68;
            case 2312   % 0x908
                BDlength = 75;
            case 2313   % 0x909
                BDlength = 81;
            case 2314   % 0x90A
                BDlength = 84;
            case 2315   % 0x90B
                BDlength = 86;
            case 2316   % 0x90C
                BDlength = 88;
            case 2317   % 0x090D
                BDlength = 89;
            case 2318   % 0x090E
                BDlength = 91;
            case 2560   % 0xA00; P Tierney & A Kell
                BDlength = 11;
            case 2817   % 0xB01; start of Joey's newer (summer 2004) formats
                BDlength = 39;
            case 2818   % 0xB02;
                BDlength = 39;
            case 2819   % 0xB03;
                BDlength = 40;
            case 2820   % 0xB04;
                BDlength = 46;
            case 2821   % 0xB05;
                BDlength = 47;
            case 2822   % 0xB06;
                BDlength = 51;
            case 2823   % 0xB07;
                BDlength = 56;
            case 2824   % 0xB08;
                BDlength = 57;
            case 2825   % 0xB09;
                BDlength = 58;
            case 2826   % 0xBA;
                BDlength = 59;
            case 2827   % 0xBB;
                BDlength = 60;
            case 2828   % 0xBC;
                BDlength = 61;
       	    case 2829   % 0xBD;
                BDlength = 62;
            % this is TEMPORARY CHANGE ONLY for HH080106 and HH092906
            % necessary when fragmenting from original events.nev file,
            % which mistakenly has the wrong BDformat number ($BD instead
            % of $BE)
%             case 2829   % 0xBD;
%                 BDlength = 63;
            case 2830   % 0xBE;
                BDlength = 63;
    	    case {2831, 2832}   % 0xBF, 0xB10;
                BDlength = 64;
            otherwise
                error('lfp_monkeyPreprocess:badformat', ...
                    'Unknown format: %d (%F)', formatnum, formatnum );
        end
    end

    % Process this BDOutputStart - is BD section complete?
    stop = bdstarts(1) + BDlength + 1;
    if stop > size(events, 1) || events(stop, 2) ~= BDOutputStop
        % The section has no end.
        warning('lfp_monkeyPreprocess:deletehead', ...
            'Deleting head fragment around timestamp %.0f.\nformatnum=0x%X', ...
            events(bdstarts(1), 1), formatnum );
        [lastgood, events] = deleteTrialFragments(events, bdstarts(1));
    else
        % The section is complete.
        [params, events] = moveOneTrial(params, events, trial, bdstarts(1), stop);
        lastgood = bdstarts(1) + 1;
        trial = trial + 1;
    end

    % update start & stop lists
    if lastgood < size(events, 1)
        bdstarts = find(events(lastgood:end, 2) == BDOutputStart) ...
            + lastgood - 1;
        bdstops =  find(events(lastgood+1:end, 2) == BDOutputStop) ...
            + lastgood;
    elseif lastgood == size(events, 1)
        bdstarts = [];
        bdstops = [];
    else
        error('lfp_monkeyPreprocess:pointerError', 'oops');
    end

    % Check for a tail end with no beginning:
    if ~isempty(bdstops) && ( isempty(bdstarts) ...
            || bdstops(1) < bdstarts(1) )
        % Found tail with no beginning
        warning('lfp_monkeyPreprocess:deletetail2', ...
            'Deleting tail fragment around timestamp %g.', ...
            events(bdstops(1), 1) );
        [lastgood, events] = deleteTrialFragments(events, bdstops(1));
        % update start & stop lists
        if lastgood < size(events, 1)
            bdstarts = find(events(lastgood:end, 2) == BDOutputStart) ...
                + lastgood - 1;
            bdstops =  find(events(lastgood+1:end, 2) == BDOutputStop) ...
                + lastgood;
        elseif lastgood == size(events, 1)
            bdstarts = [];
            bdstops = [];
        else
            error('lfp_monkeyPreprocess:pointerError2', 'ups');
        end
    end
end
close(hWaitBar);



function [startdel, events] = deleteTrialFragments(events, refevent)
% <refevent> is a row index into events. In events, deletes all
% events from and including the last BDOutputStop before refevent, and to
% but not including the first BDOutputStop at or after refevent.  Returns
% the new index into events of the the BDOutputStop that was not
% deleted (which is the same as the old index of the BDOutputStop that was
% deleted, or 1 if there was no such).
%
% Special cases:
% 1. If there is no BDOutputStop before refevent, that means that
% events starts with a tail end fragment. In this case, the following
% BDOutputStop IS deleted, or else we get stuck with a tail end fragment
% that we can't get rid of.
% 2. If there is no BDOutputStop at or after refevent, then events ends
% with an unfinished BD section.  Since there is no BDOutputStop following
% the fragment, it is neither necessary nor possible to replace the
% preceding BDOutputStop with the following one (so we don't).
lfp_declareGlobals;
if refevent > size(events, 1)
    error('deleteTrialFragments:callError', 'oops');
end
BDOutputStop = 97;

precedingstops = find(events(1:refevent-1, 2) == BDOutputStop);
if isempty(precedingstops)
    startdel = 1;
else
    startdel = precedingstops(end);
end
followingstops = find(events(refevent:end, 2) == BDOutputStop) ...
    + refevent - 1;
if isempty(followingstops)
    events(startdel+1:end, :) = [];
elseif startdel == 1
    events(startdel:followingstops(1), :) = [];
else
    events(startdel:followingstops(1)-1, :) = [];
end


function saccCode = getSaccadeCode(codeNum)
%	function saccCode = getSaccadeCode(codeNum)
%	This function converts the codeNum to saccade directions.
%	The function is mainly used in readDataFile.m
%	0 -> up, 1-> right, 2-> down, 3-> left.
%	Parameters:
%	codeNum  -- code number for the direction two consecutive saccades.
%	saccCode -- 1D array of two numbers for each directions.
%	saccCodeString -- letters for the directions.
%
%	Written by Dezhe Jin, djin@mit.edu
%	Date 1/6/2003.
%   Hacked by Dan Gibson 2/2/10.

if codeNum == 0
	saccCode = [0 1];
	%saccCodeString = 'UR';
elseif codeNum == 1	
	saccCode = [1 2];
	%saccCodeString = 'RD';
elseif codeNum == 2	
	saccCode = [2 3];
	%saccCodeString = 'DL';
elseif codeNum == 3	
	saccCode = [3 0];
	%saccCodeString = 'LU';
elseif codeNum == 4	
	saccCode = [1 0];
	%saccCodeString = 'RU';
elseif codeNum == 5	
	saccCode = [2 1];
	%saccCodeString = 'DR';
elseif codeNum == 6	
	saccCode = [3 2];
	%saccCodeString = 'LD';
elseif codeNum == 7	
	saccCode = [0 3];
	%saccCodeString = 'UL';
else
	saccCode = [-1 -1];
	%saccCodeString = 'WW';	% wrong code number!
end
