function dg_makeDummyCSC(evtfilepath, cscfilepath)
%dg_makeDummyCSC makes fake data files.
%dg_makeDummyCSC(evtfilepath, cscfilepath)
%   Reads the Neuralynx events file at <evtfilepath>, and creates a new CSC
%   file in the same directory containing frames full of zeros that fully
%   enclose all the timestamps in the events file.  If no <cscfilepath> is
%   specified or if it is empty, then the new file is created with name
%   'LFP99.dat' in the same directory as the events file.  If <evtfilepath>
%   has the extension of a Rodent Cluster file, then the events are read
%   from it, and if the session directory does not exist then it is
%   created; the fake CSC file is written to the session directory; if
%   the Events.Nev file does not exist then a fake one is created and a
%   README.TXT file is also created.  If the Events.Nev file does exist,
%   then an error is raised because you should have read the events from
%   there instead.

%$Rev: 207 $
%$Date: 2014-10-16 19:07:56 -0400 (Thu, 16 Oct 2014) $
%$Author: dgibson $

rodentclusterflag = false;
[pathstr,name,ext] = fileparts(evtfilepath);
if ~isempty(regexpi(ext, '^\.TT[0-9]$')) ...
        || ~isempty(regexpi(ext, '^\.T[0-9][0-9]$'))
    rodentclusterflag = true;
    sessiondir = fullfile(pathstr,name);
    if nargin < 2
        cscfilepath = fullfile(sessiondir, 'LFP99.dat');
    end
    fakeEvtFilePath = fullfile(sessiondir, 'Events.Nev');
    if exist(sessiondir) == 0
        mkdir(sessiondir);
    else
        if exist(fakeEvtFilePath)
            error('dg_makeDummyCSC:fileexists', ...
                'File %s already exists.', fakeEvtFilePath);
        end
    end
    [FileHeader, TrialData] = dg_ReadRodentFormat(evtfilepath);
    eventses = cell(size(TrialData));
    [eventses{:}] = deal(TrialData.events);
    eventses = 100 * cell2mat(eventses');  % convert from 0.1 ms to us
    eventses(:,50) = [];    % always contains garbage
    % Put an ITI of 10 sec between every pair of trials:
    EventTimeStamps(1,:) = eventses(1,:);
    lastevt = max(eventses, [], 2);
    if size(eventses,1) >= 2
        offset(1,1) = lastevt(1) + 5e6;
        for k = 2:(numel(lastevt) - 1)
            offset(k,1) = offset(k-1,1) + lastevt(k) + 10e6;
        end
        EventTimeStamps(2:size(eventses,1),:) = eventses(2:end,:) + ...
            repmat(offset, 1, size(eventses,2));
    else
        offset=[];
    end
    % make into vector and get rid of non-events:
    EventTimeStamps = reshape(EventTimeStamps,1,[]);
    EventTimeStamps(find(eventses == 0)) = [];
    [I, eventids] = ind2sub(size(eventses), find(eventses));
    % We now have EventTimeStamps and eventids in eventses-column order,
    % but we want them in chronological order; also for each trial we must
    % add another recorded segment containing an event 90 to mark it as a
    % "good trial", 5 sec after last event of trial:
    goodtrialstamps = [];
    for lastevtstamp = (lastevt + [0; offset])'
        goodtrialstamps = ...
            [ goodtrialstamps ([5.000e6 5.001e6 5.002e6] + lastevtstamp) ];
    end
    goodtrialevts = repmat([1 90 2],1,numel(lastevt));
    EventTimeStamps(end+1:end+numel(goodtrialstamps)) = goodtrialstamps;
    eventids(end+1:end+numel(goodtrialevts)) = goodtrialevts;
    [B,IX] = sort(EventTimeStamps);
    EventTimeStamps = reshape(EventTimeStamps(IX),1,[]);
    eventids = reshape(eventids(IX),1,[]);
    % There is some disagreement between Neuralynx and Matlab as to whether
    % these 16-bit integers are signed; subtract 2^15 to convert to
    % Neuralynx representation: 
    Nttls = eventids - 32768;
    Extras = zeros(8,numel(EventTimeStamps));
    EventStrings = cell(numel(EventTimeStamps),1);
    eventids = eventids + hex2dec('8000');  % set "strobe" bit
    for k = 1:numel(EventTimeStamps)
        EventStrings{k} = sprintf('0x%04X', eventids(k));
    end
    NlxEvtIDs = repmat(119, size(EventTimeStamps));
    if ~ispc
        error('dg_makeDummyCSC:notPC', ...
            'Neuralynx format files can only be created on Windows machines.');
    end
    Mat2NlxEV(fakeEvtFilePath, EventTimeStamps, NlxEvtIDs, Nttls, ...
        Extras, EventStrings, numel(EventTimeStamps));
    readmefid = fopen(fullfile(sessiondir, 'README.TXT'), 'a');
    fprintf(readmefid, 'The Events.Nev file in this directory is a fake created\r\n');
    fprintf(readmefid, 'by dg_makeDummyCSC from the Rodent Cluster file %s.\r\n', ...
        evtfilepath);
    fclose(readmefid);
else
    if nargin < 2
        cscfilepath = fullfile(pathstr, 'LFP99.dat');
    end
    [EventTimeStamps, EventStrings] = Nlx2MatEV(evtfilepath, ...
        [1, 0, 0, 0, 1], 0, 1);
end
duration = max(EventTimeStamps) - min(EventTimeStamps) + 1e6;
starttime = min(EventTimeStamps) - 0.5e6;
if starttime < 0
    error('oops... I assumed wrong!');
end
sampleperiod = 1000;    % 1 ms = 1000 us
numrecs = ceil(duration/(512*sampleperiod));
TimeStamps = (0:(numrecs-1)) * 512 * sampleperiod + starttime;
ChannelNumbers = zeros(1,numrecs);
SampleFrequencies = repmat(1000,1,numrecs);
NumberValidSamples = repmat(512,1,numrecs);
samples = zeros(512,numrecs);
if ~ispc
    error('dg_makeDummyCSC:notPC', ...
        'Neuralynx format files can only be created on Windows machines.');
end
Mat2NlxCSC(cscfilepath, TimeStamps, ChannelNumbers, ...
    SampleFrequencies, NumberValidSamples, samples, numrecs);
