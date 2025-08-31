function dg_WriteRodentFormat(filename, FileHeader, TrialData)
%dg_WriteRodentFormat is the inverse of dg_ReadRodentFormat.  Writes
%   Yasuo's format 0xFFF3.  The FileHeader.Format field is ignored.

%dg_WriteRodentFormat(filename, FileHeader, TrialData)
% FileHeader is a structure containing the fields of the Delphi Mouse
% Analysis "FileRec" structure.  See comments about FileHeader in the code
% below. TrialData is a 1-D array of structures, where each 'trial'
% structure has the following form:
%   trial.header - a sructure that contains all the fields of the Mouse 
%                 Analysis TrialRec, except that SSize1 - SSize10 have been
%                 replaced by one array (SSize); ditto for Free1 - Free5.
%   trial.spikes - a 2-D array, where spikes(i, j) contains the ith
%                 timestamp of cluster j (i.e. each cluster's spike data is
%                 represented by a column vector); note that columns are
%                 padded with zeros at the end to make the array
%                 rectangular.
%   trial.events - a 1-D array, where events(i) contains the time stamp of
%                 the event whose TTL code is i.

% Similarly to the Delphi code, <position> points at the last written location.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

global dg_WRF_data dg_WRF_MaxCluster dg_WRF_MaxEvent
dg_WRF_data = [];
dg_WRF_MaxEvent = 50;
dg_WRF_MaxCluster = 12;

if FileHeader.Format ~= 65523
    error('dg_WriteRodentFormat:badFormat', ...
        'Cannot write format 0x%s', dec2hex(FileHeader.Format) );
end
fid = fopen(filename, 'w');
if fid == -1
    error('dg_WriteRodentFormat:badfname', ...
        'Could not write file %s', filename );
end

% The FileHeader field names are the same as in Mouse Analysis, except
% where noted.  Explanatory comments are mostly copied from Mouse Analysis
% and a few printed pages of Yasuo's data-writing code.
dg_WRF_data(1) = hex2dec('FFF3');   % called 'Marker' in Mouse Analysis
dg_WRF_data(2) = FileHeader.Year;
dg_WRF_data(3) = FileHeader.Month;
dg_WRF_data(4) = FileHeader.Day;
dg_WRF_data(5) = FileHeader.SSize;   % number of spikes; might be unused
dg_WRF_data(6) = FileHeader.CSize;   % number of clusters in file
dg_WRF_data(7) = FileHeader.TSize;   % number of trials
dg_WRF_data(8) = FileHeader.ProcType;
dg_WRF_data(9) = FileHeader.RightCS;   % Right Turn CS
dg_WRF_data(10) = FileHeader.LeftCS;  % Left Turn CS
dg_WRF_data(11) = FileHeader.NoGoCS;  % NoGo Turn CS
position = 11;
NTrialsWritten = 0;
while NTrialsWritten < FileHeader.TSize
    position = WriteOneTrial(TrialData(NTrialsWritten + 1), ...
        position, FileHeader);
    NTrialsWritten = NTrialsWritten + 1;
end

count = fwrite(fid, dg_WRF_data, 'int32');
fclose(fid);


function position = WriteOneTrial(trial, position, FileHeader)
% Write one trial's worth of data (trial header, spikes) from the
% dg_WRF_data array starting at position.  Update position to
% point at last written element in dg_WRF_data.  fheader is the
% file header structure.
global dg_WRF_data
position = WriteTrialHeader(trial.header, position);
position = WriteTrialSpikes(trial.spikes, position, FileHeader, trial.header);
position = WriteEvents(trial.events, position);


function position = WriteTrialHeader(header, position);
global dg_WRF_data dg_WRF_MaxCluster
if length(header.SSize) > dg_WRF_MaxCluster
    error('dg_WriteRodentFormat:bigSSize', ...
        'T files cannot have more than %d clusters', dg_WRF_MaxCluster );
end
dg_WRF_data(position+1) = 61166;   % should always be 0xEEEE = 61166
dg_WRF_data(position+2) = header.TNumber;       % Trial number
dg_WRF_data(position+3) = header.SType;       % Stimulus Type: 1, 8, or 0
dg_WRF_data(position+4) = header.RType;       % Response Type: 1 to 6
dg_WRF_data(position+5)  = header.TSSize;       % Total number of spikes
dg_WRF_data(position+6 : ...
    position+5+dg_WRF_MaxCluster) = ...
    header.SSize(1:dg_WRF_MaxCluster);     % the numbers of spikes in each cluster
numfree = 15-dg_WRF_MaxCluster;
dg_WRF_data(position+6+dg_WRF_MaxCluster : ...
    position+5+dg_WRF_MaxCluster+numfree) = ...
    header.Free(1:numfree);    % unused space
position = position + 20;


function position = WriteTrialSpikes(spikes, position, fheader, theader)
global dg_WRF_data dg_WRF_MaxCluster
for cluster = (1:fheader.CSize)
    position = position + 1;
    % The spike data "Marker" is 0xDDDD + cluster:
    dg_WRF_data(position) = 56797 + cluster;
    dg_WRF_data(position+1 : position+theader.SSize(cluster)) = ...
        spikes(1:theader.SSize(cluster), cluster);
    position = position + theader.SSize(cluster);
end
% Write spike data markers for unused clusters (note local inversion of
% order of incrementing <position>, so <position> briefly points at the
% to-be-written location).
for cluster = (fheader.CSize + 1) : dg_WRF_MaxCluster
    position = position + 1;
    dg_WRF_data(position) = 56797 + cluster;
end


function position = WriteEvents(events, position)
global dg_WRF_data dg_WRF_MaxEvent
position = position + 1;
dg_WRF_data(position) = 48059;
dg_WRF_data(position + 1 : ...
    position + dg_WRF_MaxEvent) = ...
    events(1:dg_WRF_MaxEvent);
position = position + dg_WRF_MaxEvent;
