function [FileHeader, TrialData] = dg_ReadRodentFormat(filename, varargin)
%DG_READRODENTFORMAT Read a file in one of Yasuo's formats 0xFFF1, 0xFFF2 or
% 0xFFF3.

% FileHeader = dg_ReadRodentFormat(..., 'header')
%   Reads only the FileHeader portion of the file.
% [FileHeader, TrialData] = dg_ReadRodentFormat(filename)
% FileHeader is a structure containing the fields of the Delphi Mouse
% Analysis "FileRec" structure.  This depends on the format; see comments
% about FileHeader in the code below.
% TrialData is a 1-D array of structures, where
% each 'trial' structure has the following form:
%   trial.header - a sructure that contains all the fields of the Mouse 
%                 Analysis TrialRec, except that SSize1 - SSize10 have been
%                 replaced by one array (SSize); ditto for Free1 - Free5.
%   trial.spikes - a 2-D array, where spikes(i, j) contains the ith
%                 timestamp of cluster j (i.e. each cluster's spike data is
%                 represented by a column vector); note that columns are
%                 padded with zeros at the end to make the array
%                 rectangular.
%   trial.events - a 1-D array, where events(i) contains the time stamp of
%                 the event whose TTL code is i; for standard rat T-maze
%                 data in 0xFFF3 format, the codes are:
%                     1: Start of Recording
%                     2: End of Recording
%                     3: Baseline On    (Event 10 - 700)
%                     4: Baseline Off   (Event 10 - 200)
%                     5: Trial On  (not used, use event 1 instead)
%                     6: Trial off  (actual end of trial)
%                    10: Click Warning Cue
%                    11: Gate Opening
%                    12: Locomotion Onset
%                    13: Start
%                    14: Turn On
%                    15: Right Turn Off
%                    16: Left Turn Off
%                    17: Right Goal
%                    18: Left Goal
%                    21: Rough On
%                    22: Smooth On
%                    30: White Noise On
%                    31: 1 kHz On
%                    38: 8 kHz On
%                    40: White Noise Off
%                    41: 1 kHz Off
%                    48: 8 kHz Off

% Similarly to the Delphi code, <position> points at the last read location.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

headerflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'header'
            headerflag = true;
        otherwise
            error('dg_ReadRodentFormat:badoption', ...
                ['The option "' varargin{argnum} ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~isequal(class(filename), 'char')
    error('dg_ReadRodentFormat:badfilename', ...
        '<filename> must be a string' );
end

shared.MaxEvent = 50;
shared.MaxCluster = 10;
shared.filename = filename;
fid = fopen(shared.filename, 'r');
if fid == -1
    error('dg_ReadRodentFormat:badfname', ...
        'Could not read file %s', shared.filename );
end
if headerflag
    [shared.rawdata, count] = fread(fid, 11, 'int32');
else
    [shared.rawdata, count] = fread(fid, inf, 'int32');
end
fclose(fid);

% The FileHeader field names are the same as in Mouse Analysis, except
% where noted.  Explanatory comments are mostly copied from Mouse Analysis
% and a few printed pages of Yasuo's data-writing code.
FileHeader.Format   = shared.rawdata(1);   % called 'Marker' in Mouse Analysis
FileHeader.Year     = shared.rawdata(2);
FileHeader.Month    = shared.rawdata(3);
FileHeader.Day      = shared.rawdata(4);
FileHeader.SSize    = shared.rawdata(5);   % number of spikes; might be unused
FileHeader.CSize    = double(shared.rawdata(6));   % number of clusters in file
FileHeader.TSize    = shared.rawdata(7);   % number of trials
if FileHeader.Format == hex2dec('FFF3')
    shared.MaxCluster = 12;
end
if FileHeader.Format == hex2dec('FFF2') || ...
        FileHeader.Format == hex2dec('FFF3')
	FileHeader.ProcType = shared.rawdata(8);
	FileHeader.RightCS  = shared.rawdata(9);   % Right Turn CS
	FileHeader.LeftCS   = shared.rawdata(10);  % Left Turn CS
	FileHeader.NoGoCS   = shared.rawdata(11);  % NoGo Turn CS
    position = 11;
elseif FileHeader.Format == hex2dec('FFF1')
    FileHeader.RightCS  = shared.rawdata(8);
    FileHeader.ProcType = 1;
    position = 8;
else
    error('dg_ReadRodentFormat:BadFileFormat', ...
        'Cannot read format 0x%s.', dec2hex(double(FileHeader.Format)) );
end
if headerflag
    return
end
TrialData = [];
NTrialsRead = 0;
while NTrialsRead < FileHeader.TSize
    [trial, position] = ReadOneTrial(shared, position, FileHeader);
    TrialData = [ TrialData trial ];
    NTrialsRead = NTrialsRead + 1;
end
if position ~= length(shared.rawdata)
    error('dg_ReadRodentFormat:FileLength', ...
        'File length discrepancy at %s position %d', ...
        shared.filename, position+1);
end


function [trial, position] = ReadOneTrial(shared, position, fheader, theader)
% Read one trial's worth of data (trial header, spikes) from the
% shared.rawdata array starting at position.  Update position to point at last
% read element in shared.rawdata.  fheader is the file header structure.
[trial.header, position] = ReadTrialHeader(shared, position);
[trial.spikes, position] = ReadTrialSpikes(shared, position, fheader, ...
    trial.header);
[trial.events, position] = ReadEvents(shared, position);


function [header, position] = ReadTrialHeader(shared, position)
if position + 20 > length(shared.rawdata)
    error('dg_ReadRodentFormat:TruncatedFile1', ...
        'File ends before end of trial header');
end
header.Marker  = shared.rawdata(position+1);   % should always be 0xEEEE = 61166
if header.Marker ~= 61166
    error('dg_ReadRodentFormat:BadTrialMarker', ...
        [ 'Bad trial header marker at %s position %d', shared.filename, ...
            int2str(position+1)]);
end
header.TNumber = shared.rawdata(position+2);       % Trial number
header.SType   = shared.rawdata(position+3);       % Stimulus Type: 1, 8, or 0
header.RType   = shared.rawdata(position+4);       % Response Type: 1 to 6
header.TSSize  = shared.rawdata(position+5);       % Total number of spikes
header.SSize(1:shared.MaxCluster) = double(shared.rawdata(position+6 : ...
    position+5+shared.MaxCluster));     % the numbers of spikes in each cluster
numfree = 15-shared.MaxCluster;
header.Free(1:numfree) = ...
    shared.rawdata(position+6+shared.MaxCluster : ...
    position+5+shared.MaxCluster+numfree);    % unused space
position = position + 20;


function [spikes, position] = ReadTrialSpikes(shared, position, fheader, theader)
spikes = [];
if position + shared.MaxCluster + sum(theader.SSize(1:fheader.CSize)) ...
        > length(shared.rawdata)
    error('dg_ReadRodentFormat:TruncatedFile2', ...
        'File ends before end of spike data');
end
for cluster = (1:fheader.CSize)
    position = position + 1;
    % The spike data "Marker" is 0xDDDD + cluster; verify it:
    if shared.rawdata(position) ~= 56797 + cluster
        error('dg_ReadRodentFormat:BadSpikeMarker', ...
            'Bad spike data marker at %s position %d', ...
            shared.filename, int2str(position));
    end
    spikesSize = size(spikes);
    spikes = [ spikes zeros(spikesSize(1), 1) ];
    if theader.SSize(cluster) > spikesSize(1)
        spikes = [ spikes
            zeros(theader.SSize(cluster) - spikesSize(1), cluster) ];
    end
    spikes(1:theader.SSize(cluster), cluster) = ...
        shared.rawdata(position+1 : position+theader.SSize(cluster));
    position = position + theader.SSize(cluster);
end
% Skip spike data markers for unused clusters
position = position + shared.MaxCluster - fheader.CSize;


function [events, position] = ReadEvents(shared, position)
if position + 1 + shared.MaxEvent > length(shared.rawdata)
    error('dg_ReadRodentFormat:TruncatedFile3', ...
        'File ends before end of event data');
end
position = position + 1;
if shared.rawdata(position) ~= 48059
    error(['Bad event data marker at ', shared.filename, ' position ', int2str(position+1)]);
end
events(1:shared.MaxEvent) = shared.rawdata(position + 1 : position + shared.MaxEvent);
position = position + shared.MaxEvent;
