function [params, events, fnameEV, badtrials] = lfp_readEvents(read_mode, ...
    Filename1, Filename2, alwaysuseVTflag, silentflag, usestrobe)
% <events> is in same format as lfp_Events, except timestamps are in
% microseconds.
% <badtrials> contains trials detected as bad while processing the events
% file.

% Sets lfp_DataDir as a side-effect if lfp_UseFileSelect.
% Note that Filename1 and Filename2 are only required if
% ~lfp_UseFileSelect.  In that case, if read_mode is 'preset', only
% Filename1 is attempted anyway.
% All returned values will be empty if the user cancels the operation via
% the GUI.

% In cases where lfp_SetupType is not one of
%           'monkey'
%           'theresa'
%           'naotaka'
%           'rodent'
%           'wheel'
% <usestrobe> determines whether the events file is read according to the
% old Neuralynx convention in which a particular bit of the TTL value is
% used to trigger the recording of an event on its transition from 0 to 1
% (usestrobe = true), or the new convention in which an event is recorded
% whenever any bit changes (usestrobe = false).  The new convention
% implicitly requires that there be a specific TTL value that is used as a
% neutral reference state from which at least one bit changes in order to
% record an event.  I hereby define that reference value to be 0, which is
% a problematic event code for Matlab purposes anyway.  When <usestrobe> is
% false, events with zero TTL value are removed silently.

%$Rev: 422 $
%$Date: 2023-09-08 14:25:38 -0400 (Fri, 08 Sep 2023) $
%$Author: dgibson $

lfp_declareGlobals;

events = [];
params = [];
fnameEV = [];
badtrials = [];

if nargin < 4
    alwaysuseVTflag = false;
end
if nargin < 5
    silentflag = false;
end
if nargin < 6
    usestrobe = true;
end

if lfp_UseFileSelect
    selectedfile = lfp_fileselect('evt', 'Please select an events file:');
    if isempty(selectedfile)
        return
    else
        if exist('lfp_fileselect.mat','file') == 2
            load('lfp_fileselect.mat');
        else
            error('lfp_readEvents:nofsfile', ...
                'lfp_fileselect.mat does not exist.' );
        end
        switch char(lfp_fileselect_FileTypesPopup_string)
            case 'Events (saved - *.EVTSAV)'
                ev_read_mode = 'evtsav';
            case 'Events (new - *.NEV)'
                ev_read_mode = 'nlx';
            case 'Events (old - *.DAT)'
                ev_read_mode = 'nlx';
            case 'Events (Matlab - *.MAT)'
                ev_read_mode = 'mat';
        end
    end
    fnameEV = fullfile(lfp_DataDir, char(selectedfile));
else
    % Automatically find Events file.
    switch read_mode
        case 'mat'
            % Filename1 is first preference, but change to .MAT
            fnameEV=fullfile(lfp_DataDir, Filename1);
            [pathstr,name] = fileparts(fnameEV);
            fnameEV = fullfile(pathstr,[name,'.MAT']);
        otherwise
            % Filename1 is first preference
            fnameEV=fullfile(lfp_DataDir, Filename1);
    end
    if ~(exist(fnameEV, 'file'))
        switch read_mode
            case 'preset'
                error('lfp_readEvents:badPreset1', ...
                    'Could not find events file "%s".\n', fnameEV );
            otherwise
                % In desperation, try Filename2
                fnameEV=fullfile(lfp_DataDir, Filename2);
                if strcmp(read_mode, 'mat')
                    [pathstr,name] = fileparts(fnameEV);
                    fnameEV = fullfile(pathstr,[name,'.MAT']);
                end
                if ~(exist(fnameEV, 'file'))
                    error('lfp_readEvents:notfound', ...
                        'Could not find alternate events file "%s".\n', ...
                        fnameEV );
                end
        end
    end
    ev_read_mode = read_mode;
end
if ~silentflag
    fprintf(1, 'Opening events file "%s".\n', fnameEV);
end
if strcmp(ev_read_mode, 'preset')
    [pathstr,name,ext] = fileparts(fnameEV); %#ok<ASGLU>
    switch lower(ext)
        case {'.nev' '.dat'}
            ev_read_mode = 'nlx';
        case '.mat'
            ev_read_mode = 'mat';
        case '.evtsav'
            ev_read_mode = 'evtsav';
        otherwise
            error('lfp_readEvents:badPreset', ...
                'Preset events files must be .NEV, .DAT, .EVTSAV, or .MAT' );
    end
end
if strcmp(ev_read_mode, 'evtsav')
    load(fnameEV, '-mat');
    events = lfp_save_events;
    clear lfp_save_events;
    events(:,1) = round(1e6 * events(:,1));
    params = lfp_save_params;
    clear lfp_save_params;
end
if ~strcmp(ev_read_mode, 'evtsav') || alwaysuseVTflag
    % Note that alwaysuseVTflag and ~strcmp(ev_read_mode, 'evtsav') may
    % BOTH be true at the same time.
    if ~strcmp(ev_read_mode, 'evtsav')
        % Process event data. In the case of an 'nlx' file,TS and TTL
        % values are used (case 'nlx' above). In the case of a 'mat' file,
        % string processing is necessary (case 'mat' above).
        switch lower(ev_read_mode)
            case 'nlx'
                try
                    [TS, TTL] = dg_readEvents(fnameEV);
                catch %#ok<CTCH>
                    error('lfp_read2:badEV', ...
                        'Failed to read file %s', fnameEV);
                end
                if isempty(TTL)
                    error('lfp_readEvents:noevents', ...
                        'There are no events in the event file.' );
                end
                % Check TTL for entries outside of the legitimate range
                % 0x8000 - 0xFFFF (i.e. defective events). Negative
                % numbers are due to sign bit propagation in Matlab.
                % Note that in JF data, zero is a legitimate event code
                % inside the BD sections, whereas in TMD data or rodent
                % data it is never legitimate.
                unstrobedEvents = TTL > -1;
                zeroEvents = TTL == -32768;
                
                if usestrobe && any(unstrobedEvents)
                    warning('lfp_readEvents:defectiveEvent',...
                        ['The events file contains %d unstrobed events.\n' ...
                        'See ''usestrobe'' in ''lfp_getEvtIDs''.\n' ...
                        'The first defective event is #%d.\n'], ...
                        sum(unstrobedEvents),find(unstrobedEvents,1));
                    % Remove unstrobed events (zero events are not
                    % removed).
                    TS = TS(~unstrobedEvents);
                    TTL = TTL(~unstrobedEvents);
                end
                
                if usestrobe && any(zeroEvents)
                    warning('lfp_readEvents:zeroEvent',...
                        ['The events file contains %d zero event(s). ', ...
                        'The first zero event is #%d.\n'], ...
                        sum(zeroEvents),find(zeroEvents,1));
                end
                
                if usestrobe
                    % TTL IDs are returned as signed integers, with 2^15
                    % bit (strobe/sync) propagated to the left.  To extract
                    % the lower 15 bits, add 2^15 to negative values to
                    % yield positive integers.
                    events = [TS' TTL'+2^15];
                else
                    % We are not using the strobe bit, which means that
                    % the entire TTL value should be interpreted as an
                    % unsigned integer.  That means we must undo the
                    % damage caused by the Nlx2Mat function's returning
                    % it as a signed integer.  Also, we must get rid of
                    % the TTL==0 events that we assume are intercalated
                    % between actual events markers, except in the case of
                    % Georgios data which uses neither strobe NOR
                    % interacalated 0s (and therefore drops repeated event
                    % IDs).
                    TTL(TTL<0) = TTL(TTL<0) + 2^16;
                    if isequal(lfp_SetupType, 'georgios_ApAv_rev_only')
                        events = [TS' TTL'];
                    else
                        unstrobedZeroEvents = TTL == 0;
                        events = [TS(~unstrobedZeroEvents)' ...
                            TTL(~unstrobedZeroEvents)'];
                    end
                end
            case 'mat' % Output of dg_Nlx2Mat (lfp_save saves .evtsav)
                load(fnameEV, '-mat');
                TimeStamps = dg_Nlx2Mat_Timestamps;
                clear dg_Nlx2Mat_Timestamps;
                if exist('dg_Nlx2Mat_TTL', 'var')
                    clear dg_Nlx2Mat_EventStrings;
                    if usestrobe
                        nostrobe = dg_Nlx2Mat_TTL >= 0; %#ok<NODEF>
                        if any(nostrobe)
                            warning('lfp_readEvents:nosync2', ...
                                'The events file contains %d event IDs with sync bit = 0; ignoring them.', ...
                                sum(nostrobe));
                        end
                        TimeStamps(nostrobe) = [];
                        dg_Nlx2Mat_TTL(nostrobe) = [];
                        dg_Nlx2Mat_TTL = dg_Nlx2Mat_TTL + 2^15;
                    end
                    events = [reshape(TimeStamps,[],1) reshape(dg_Nlx2Mat_TTL,[],1)];
                    clear dg_Nlx2Mat_TTL;
                else
                    EventStrings = dg_Nlx2Mat_EventStrings;
                    clear dg_Nlx2Mat_EventStrings;
                    % Do string processing stuff on TimeStamps and EventStrings here in
                    % the case of a .mat file
                    CharStrings = char(EventStrings);
                    if isempty(CharStrings)
                        error('lfp_readEvents:noevents', ...
                            'There are no events in the event file.' );
                    end
                    % Test for and delete any defective events, with warnings.
                    % First, check for non-hex characters (including blank, which appears
                    % when there are fewer than four hex digits):
                    putativedigits = CharStrings(:, end-3:end); 
                    [nondigitr, nondigitc] = find( ...
                        ~ismember(putativedigits, '0123456789ABCDEF') ); %#ok<NASGU>
                    if ~isempty(nondigitr)
                        rows = unique(nondigitr);
                        warning('lfp_readEvents:badID', ...
                            ['The events file contains %d corrupted ID strings in %d ' ...
                            'records total.\nThe first corrupted ID is record #%d.' ], ...
                            length(rows), length(TimeStamps), rows(1) );
                        TimeStamps(nondigitr) = [];
                        putativedigits(nondigitr,:) = [];
                    end
                    % Check for events that do not have the highest ("sync" or
                    % "strobe") bit set. Takashi's wheel setup is missing 2^7 bit, so
                    % the strobe bit ends up in 2^14 bit instead of 2^15 bit.
                    if isequal(lfp_SetupType, 'wheel')
                        firstdigits = hex2dec(putativedigits(:, 1));
                        nosync = find(firstdigits < 4);
                        if ~isempty(nosync)
                            warning('lfp_readEvents:nosync', ...
                                'The events file contains event IDs with sync bit = 0.');
                            TimeStamps(nosync) = [];
                            putativedigits(nosync,:) = [];
                        end
                    else
                        firstdigits = hex2dec(putativedigits(:, 1));
                        nosync = find(firstdigits < 8);
                        if ~isempty(nosync)
                            warning('lfp_readEvents:nosync', ...
                                'The events file contains event IDs with sync bit = 0.');
                            TimeStamps(nosync) = [];
                            putativedigits(nosync,:) = [];
                        end
                    end
                    disp('Converting event IDs to numeric form.');
                    Evt = hex2dec(putativedigits);
                    switch lfp_SetupType
                        case 'naotaka2'
                            % Use this setupType when NOT adding Naotaka's *.DWD files,
                            % because there are often false positives on all bits of
                            % his high order byte; however, better agreement with his
                            % *.DWD is obtained by letting lfp_monkeyPreprocess get rid
                            % of the BD fragments that result when you include the high
                            % order byte garbage.
                            Evt = bitand(Evt, hex2dec('7F')); % kill entire high byte
                        case 'wheel'
                            Evt = bitand(Evt, hex2dec('3FFF')); % kill strobe bit
                        otherwise
                            if usestrobe
                                Evt = bitand(Evt, hex2dec('7FFF')); % kill strobe bit
                            else
                                zeroEvents = Evt == 0;
                                Evt(zeroEvents) = [];
                                TimeStamps(zeroEvents) = [];
                            end
                    end
                    events = [TimeStamps' Evt];
                end
        end
    end
    switch lfp_SetupType
        case {'monkey' 'theresa' 'naotaka'}
            [params, events] = lfp_monkeyPreprocess(events);
        case {'ken' 'georgios_ApAv' 'georgios_ApAvApAp' ...
                'georgios_ApAp' 'georgios_ApAv_rev_only'}
            [params, events, badtrials] = lfp_otherPreprocess(events);
        case 'wheel'
            events = lfp_wheelPreprocess(events);
            params = {};
        case 'rodent'
            events = lfp_rodentPreProcess(events);
            params = {};
        otherwise % this includes 'simple'
            params = {};
    end
end
lfp_log(sprintf('Read events file "%s"', fnameEV));
