function lfp_save_noGUI(name, ftype, varargin)
%lfp_save_noGUI(name, ftype, varargin)
% Replacement for lfp_save('preset', ...).
%INPUTS
%   name: string specifying filename to save to.  <name> can be  anything
%       when saving event files, and the resulting file is named
%       "<name>.evtsav".  When saving other Matlab file types, the
%       resulting file is named "<name>.mat". For CSC files, <name> must be
%       one of the names of the CSC files currently in memory.   When
%       saving Rodent Cluster files, <name> must be a string or cell string
%       array containing names of spike channels currently in memory; if
%       more than one is given, the first one is used as the basis for the
%       filename being saved.
%   ftype: determines what type of file to save (the value of <ftype> is
%       not case sensitive). 
%       'evt' - saves events
%       'csc' - saves CSC data
%       'rod' - saves "Yasuo format" Rodent Clusters *.Tnn file (see
%           dg_WriteRodentFormat)
%       'mcc' - saves Multiple Cut Clusters
%OPTIONS
%   'destdir', destdir - <destdir> is used in place of lfp_DataDir as the
%       location to which to write the file when using 'preset' option.  If
%       <destdir> does not exist, it is created.
%   'double' - the opposite of 'int', i.e. does not round the values; this
%       is the default.
%   'int' - rounds CSC values to integers, potentially making the saved
%       file roughly up to eight times smaller (depending on the values).
%       Matlab is allowed to make its own decisions how to store the data,
%       so if the data stored as doubles in memory are already integer
%       values, specifying 'int' will not change anything.
%   'noclick' - requires an additional argument, <taskparams>.  Bypasses
%       the GUI for specifying task parameters when saving Rodent Clusters.
%       <taskparams> is the vector of Rodent Clusters file header values:
%       [ ProcType RightCS LeftCS NoGoCS ].
%   'selectedtrials' - save only trials that are not on the lfp_BadTrials
%       list, and for which lfp_SelectedTrials is true.  Works only for
%       saving Rodent Cluster files.
%   'single' - converts values to single-precision float before saving.
%NOTES
% In the distant mists of ancient history, this was a GUI-based function
% running only under Windows, hence the Windows-like use of uppercase
% filename extensions in the <filetype> labels.  However, the actual files
% saved have lowercase extensions.

%$Rev: 400 $
%$Date: 2019-07-16 14:14:13 -0400 (Tue, 16 Jul 2019) $
%$Author: dgibson $

lfp_declareGlobals;

handles.output = [];

argnum = 1;
handles.destdir = false;
handles.datatype = 'double';
handles.ftype = ftype;
handles.name = name;
handles.selectedtrialsflag = false;
handles.noclickflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'destdir'
            argnum = argnum + 1;
            handles.destdir = varargin{argnum};
        case 'double'
            handles.datatype = 'double';
        case 'int'
            handles.datatype = 'int';
        case 'single'
            handles.datatype = 'single';
        case 'selectedtrials'
            handles.selectedtrialsflag = true;
        case 'noclick'
            if length(varargin) < argnum + 1
                error('"noclick" option requires array input of form [proctype right_stim left_stim no_go_stim]');
            end
            handles.noclickflag=true;
            handles.inputtaskparams = varargin(argnum + 1);
            argnum = argnum + 1;
        otherwise
            error('lfp_save:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

%     SaveFileButton_Callback(handles.SaveFileButton, eventdata, handles);
switch lower(handles.ftype)
    case 'evt'
        filetype = 'Events (*.EVTSAV)';
    case 'csc'
        filetype = 'CSC Channel (*.MAT)';
    case 'rod'
        filetype = 'Rodent Clusters (*.Tnn)';
    case 'mcc'
        filetype = 'Multiple Cut Cluster (*.MAT)';
    otherwise
        error('lfp_save:badpresetext', ...
            'The ftype given with ''preset'' option is not recognized' );
end
if isequal(handles.destdir, false)
    OutputPathName = lfp_DataDir;
else
    OutputPathName = handles.destdir;
    if ~isempty(OutputPathName) && ~exist(OutputPathName, 'dir')
        mkdir(OutputPathName);
    end
end
switch filetype
    case 'CSC Channel (*.MAT)'
        channelnames = lfp_FileNames(lfp_ActiveFilenums); %#ok<*USENS>
        channels = find(ismember(channelnames, handles.name));
        if isempty(channels)
            error('lfp_save:badpresetname', ...
                'The name given with the preset option must match a loaded CSC channel' );
        end
        dg_Nlx2Mat_Timestamps = round(1e6*lfp_TimeStamps); %#ok<*NASGU>
        fprintf('Saving filenums %s\n', ...
            dg_thing2str(lfp_ActiveFilenums(channels)));
        for filenum = lfp_ActiveFilenums(reshape(channels, 1, []))
            OutputFileName = [lfp_FileNames{filenum} '.mat'];
            extrasamples = rem(numel(lfp_Samples{filenum}), lfp_SamplesPerFrame);
            padding = [];
            if extrasamples ~= 0
                msgbox(sprintf(...
                    'NOTE: padded "%s" with zeros to fill last frame', OutputFileName ));
                padding = zeros(1, lfp_SamplesPerFrame-extrasamples);
            end
            dg_Nlx2Mat_Samples = reshape(...
                [lfp_Samples{filenum}(1:end), padding], ...
                lfp_SamplesPerFrame, [] );
            switch handles.datatype
                case 'int'
                dg_Nlx2Mat_Samples = round(dg_Nlx2Mat_Samples);
                case 'single'
                dg_Nlx2Mat_Samples = single(dg_Nlx2Mat_Samples);
                case 'double'
                % do nothing
                otherwise
                    error('Doh!');
            end
            dg_Nlx2Mat_SamplesUnits = lfp_SamplesUnits{filenum};
            save(fullfile(OutputPathName, OutputFileName), ...
                'dg_Nlx2Mat_Timestamps', 'dg_Nlx2Mat_Samples', ...
                'dg_Nlx2Mat_SamplesUnits', '-v7.3' );
            lfp_log(sprintf('Saved CSC file %s', fullfile(OutputPathName, OutputFileName)));
        end
    case 'Events (*.EVTSAV)'
        OutputFileName = [handles.name '.evtsav'];
        lfp_save_events = lfp_Events;
        lfp_save_params = lfp_TrialParams;
        save(fullfile(OutputPathName, OutputFileName), ...
            'lfp_save_events', 'lfp_save_params');
        lfp_log(sprintf('Saved Events file %s', fullfile(OutputPathName, OutputFileName)));
    case 'Rodent Clusters (*.Tnn)'
        % Because all time stamps are saved relative to the estimated start
        % of recording for each trial, this "should work" equally well on
        % any session when multiple sessions are loaded.
        channelnames = lfp_SpikeNames;
        channels = find(ismember(channelnames, handles.name));
        if isempty(channels)
            error('lfp_save:badpresetname2', ...
                'The name given with the preset option must match a loaded spike channel' );
        end
        fprintf('Saving spike channels %s\n', dg_thing2str(channels));
        OutputFileName = makeClusterFilename('t', ...
            lfp_SpikeNames{channels(1)});
        if handles.selectedtrialsflag
            trials = lfp_enabledTrials(1:size(lfp_TrialIndex,1));
        else
            trials = 1:size(lfp_TrialIndex,1);
        end
        % Create file header:
        FileHeader.Format = hex2dec('FFF3');
        clockval = clock;
        FileHeader.Year     = clockval(1);
        FileHeader.Month    = clockval(2);
        FileHeader.Day      = clockval(3);
        FileHeader.SSize    = 0;
        channels = reshape(channels, 1, []);
        if handles.noclickflag
            taskparams = cell2mat(handles.inputtaskparams);
        else
            taskparams = lfp_getRodentTaskParams;
        end
        FileHeader.CSize    = length(channels);   % number of clusters in file
        FileHeader.TSize    = length(trials);   % number of trials
        FileHeader.ProcType = taskparams(1);
        FileHeader.RightCS  = taskparams(2);
        FileHeader.LeftCS   = taskparams(3);
        FileHeader.NoGoCS   = taskparams(4);
        % Re-package spikes & event data according to trial structure, and
        % save the file.  Note that trial numbers are strictly sequential,
        % and may thus disagree with the trial numbers that are saved by
        % Yasuo's programs.  Spike and event times are converted to
        % tenths of milliseconds relative to a fictitious starting time
        % that is arbitrarily set to coincide with the start of the
        % "pre-roll" period.  The actual nominal trial start time cannot be
        % used because then the start event would have a zero timestamp and
        % appear in the Yasuo format to be missing.
        TrialData = [];
        preRoll = 2.0012;
        postRoll = 0.5;
        for trialidx = 1:length(trials)
            trial = trials(trialidx);
            spikes = lfp_getTrialSpikes(trial, channels, preRoll, postRoll);
            for channelidx = 1:length(channels)
                FileHeader.SSize = ...
                    FileHeader.SSize + length(spikes{channelidx});
            end
            TrialData(trialidx).header.TNumber = trial; %#ok<*AGROW>
            % 12 is value of dg_WRF_MaxCluster from dg_WriteRodentFormat:
            TrialData(trialidx).header.SSize = zeros(1, 12);
            for cluster = 1:FileHeader.CSize
                TrialData(trialidx).header.SSize(cluster) = ...
                    numel(spikes{cluster});
            end
            TrialData(trialidx).header.TSSize = ...
                sum(TrialData(trialidx).header.SSize);
            TrialData(trialidx).header.Free = [1 1 1];
            paddedspikes = zeros(max(TrialData(trialidx).header.SSize), ...
                length(spikes));
            reftime = lfp_Events(lfp_TrialIndex(trial,1),1) - preRoll;
            for cluster = 1:FileHeader.CSize
                paddedspikes(1:length(spikes{cluster}), cluster) = ...
                    round(1e4*(spikes{cluster} - reftime));
            end
            TrialData(trialidx).spikes = paddedspikes;
            % Create Yasuo-style events array, assuming there is no more
            % than one of each event ID per trial.  Include the event
            % before lfp_NominalTrialStart if it is ID 1 (Record On).
            maxEvent = 50; % value of dg_WRF_MaxEvent from dg_WriteRodentFormat
            eventsarray = zeros(1, maxEvent);
            startevtidx = lfp_TrialIndex(trial,1);
            if startevtidx > 1 && lfp_Events(startevtidx - 1, 2) == 1
                startevtidx = startevtidx - 1;
            end
            for evtidx = startevtidx : lfp_TrialIndex(trial,2)
                eventsarray(lfp_Events(evtidx,2)) = round(...
                    1e4*(lfp_Events(evtidx,1) - reftime) );
            end
            if length(eventsarray) < maxEvent
                eventsarray(end+1:maxEvent) = 0;
            end
            TrialData(trialidx).events = eventsarray(1:maxEvent);
            % These are correct SType values when there is an event [31
            % 38 21 22], based on \\Matrisome2\RData\A75\acq02\EACQ02.DAT:
            if eventsarray(31) || eventsarray(21)
                TrialData(trialidx).header.SType = 1;
            elseif eventsarray(38) || eventsarray(22)
                TrialData(trialidx).header.SType = 8;
            else
                TrialData(trialidx).header.SType = 0;
            end
            % These RType values are correct regarding direction:
            if eventsarray(15)
                % Right Turn
                if TrialData(trialidx).header.SType == FileHeader.RightCS
                    % correct
                    TrialData(trialidx).header.RType = 1;
                else
                    % incorrect
                    TrialData(trialidx).header.RType = 4;
                end
            elseif eventsarray(16)
                % Left Turn
                if TrialData(trialidx).header.SType == FileHeader.LeftCS
                    % correct
                    TrialData(trialidx).header.RType = 2;
                else
                    % incorrect
                    TrialData(trialidx).header.RType = 5;
                end
            else
                % Score any other response as "incomplete"; note this does
                % not cover the no-go task:
                TrialData(trialidx).header.RType = 0;
            end
        end
        dg_WriteRodentFormat(fullfile(OutputPathName, OutputFileName), ...
            FileHeader, TrialData);
        lfp_log(sprintf('Saved Rodent Clusters file %s', ...
            fullfile(OutputPathName, OutputFileName)));
    case 'Multiple Cut Cluster (*.MAT)'
        % For compatibility with lfp_add, file must contain array with
        % spike timestamps in seconds in col 1 and cluster ID in col 2 as
        % first variable saved (does not care what its name is)
        channelnames = lfp_SpikeNames;
        channels = find(ismember(channelnames, handles.name));
        if isempty(channels)
            error('lfp_save:badpresetname3', ...
                'The name given with the preset option must match a loaded spike channel' );
        end
        fprintf('Saving spike channels %s\n', dg_thing2str(channels));
        OutputFileName = makeClusterFilename('m', lfp_SpikeNames{channels(1)});
        numrows = 0;
        for ch = reshape(channels, 1, [])
            numrows = numrows + numel(lfp_Spikes{ch});
        end
        lfp_save_spikes = zeros(numrows, 2);
        nextrow = 1;
        arbseqnum = 1000;
        for ch = reshape(channels, 1, [])
            idx = regexp(lfp_SpikeNames{ch}, 'C\d+$');
            if isempty(idx)
                arbseqnum = arbseqnum + 1;
                clust = arbseqnum;
            else
                clust = str2double(lfp_SpikeNames{ch}(idx+1:end));
            end
            numspikes = numel(lfp_Spikes{ch});
            lfp_save_spikes(nextrow : nextrow + numspikes - 1, :) = ...
                [ reshape(lfp_Spikes{ch}, [], 1) ...
                repmat(clust, numspikes, 1) ];
            nextrow = nextrow + numspikes;
        end
        save(fullfile(OutputPathName, OutputFileName), 'lfp_save_spikes');
        lfp_log(sprintf('Saved Multiple Cut Clusters file %s', ...
            fullfile(OutputPathName, OutputFileName)));
    otherwise
        error('lfp_save:internalError2', 'programming error');
end
end

