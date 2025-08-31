function TrialData = dg_ReadMouseEventFormat(filename)
%DEPRECATED.  Use dg_ReadRodentFormat instead.

% [FileHeader, TrialData] = dg_ReadMouseEventFormat(filename)
% Return values are as in dg_ReadMouseFormat, except that everything except
% the event data and the trial header fields TNumber, Stype, Rtype is
% returned empty.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

MaxEvent = 50;
fid = fopen(filename, 'r');

TrialData = [];
NTrial = 0;

line = [];
while length(line) ~= 1 || line ~= -1   % while not EOF
    line = fgetl(fid);
    if length(line) == 0
        continue
    end
    if line(1) == 'T'
        % begin a new trial
        if NTrial > 0
            TrialData = [TrialData trial];
        end
        NTrial = NTrial + 1;
        trial.events(1:MaxEvent) = 0;
        tdata = sscanf(line(6:end), '%d', 3);
        trial.header.TNumber = tdata(1);
        trial.header.SType = tdata(2);
        trial.header.RType = tdata(3);
    elseif line(1) == 'E'
        % read one event
        [token, rem] = strtok(line, 'E ');
        timestamp = str2num(token);
        token = strtok(rem);
        eventID = str2num(token);
        trial.events(eventID) = timestamp;
    end
end
TrialData = [TrialData trial];

fclose(fid);
