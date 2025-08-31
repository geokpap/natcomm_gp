function [TimeStamps, TTLIDs, trialends] = lfp_readYasuoEvents(filename)
%[TimeStamps, TTLIDs] = readYasuoEvents(filename)
% Returns the timestamp and TTL ID of every event in the file that has a
% nonzero timestamp.  <trialends> contains the index into TimeStamps and
% TTLIDs of the last event of each trial.  Return values are all row
% vectors.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

MaxEvent = 50;
fid = fopen(filename, 'r');

TimeStamps = [];
TTLIDs = [];
trialends = [];
NTrial = 0;

line = [];
while length(line) ~= 1 || line ~= -1
    line = fgetl(fid);
    if length(line) > 1 && isequal(line(1:2), 'E ')
        [token, rem] = strtok(line, 'E ');
        timestamp = str2num(token);
        TimeStamps = [ TimeStamps timestamp ];
        token = strtok(rem);
        eventID = str2num(token);
        TTLIDs = [ TTLIDs eventID ];
    elseif length(line) > 4 && isequal(line(1:5), 'Trial')
        trialends = [ trialends length(TTLIDs) ];
    end
end
trialends = [trialends(2:end) length(TTLIDs) ];

fclose(fid);
