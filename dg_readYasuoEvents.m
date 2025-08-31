function [header, trialdata] = dg_readYasuoEvents(filepath)
%   Written to format of RVTrack7 output, but contents of trialdata should
%   be valid regardless of format.
%   Returns entire contents of Yasuo format event file.  <header> is a 1x1
%   structure, trialdata is header.Numtrials x 1 structure.  <header>
%   fields that could not be converted are left empty.
%   header.AnimalID
%   header.Session
%   header.Procedure {Procedure type: 1 = R/L, 2 = R/L/NG, 3 = ET R/L}
%   header.Program {Program used to create events file}
%   header.Date = [Month, Day, Year]
%   header.Numtrials
%   header.Rstim
%   header.Lstim
%   header.NGstim {No-Go stimulus}
%   trialdata.trialnum
%   trialdata.stimtype
%   trialdata.response
%   trialdata.events = 2-column array, 1 row per event; col 1 = timestamp,
%       col 2 = eventID

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

fid = fopen(filepath, 'r');
if fid == -1
    error('dg_readYasuoEvents:open', ...
        'Could not open file "%s"', filepath );
end
line = fgetl(fid);
linenum = 1;
trialdata = [];
while ~isequal(line, -1)
    if length(line)>4 && isequal(line(1:5), 'Trial')
        A = sscanf(line, 'Trial %d %d %d');
        if ~isequal(size(A), [ 3 1 ])
            warning('dg_readYasuoEvents:badline', ...
                'Bad trial header at line %d', linenum);
        else
            trialdata(end+1).trialnum = A(1);
            trialdata(end).stimtype = A(2);
            trialdata(end).response = A(3);
            trialdata(end).events = zeros(0, 2);
        end
    elseif length(line)>1 && isequal(line(1:2), 'E ')
        A = sscanf(line, 'E %d %d %d');
        if ~isequal(size(A), [ 2 1 ])
            warning('dg_readYasuoEvents:badline2', ...
                'Bad event record at line %d', linenum);
        else
            trialdata(end).events(end+1,:) = A';
        end
    else
        switch linenum
            case 2
                header.AnimalID = strtok(line);
            case 3
                header.Session = strtok(line);
            case 4
                header.Procedure = sscanf(line, '%d');
            case 5
                header.Program = strtok(line);
            case 6
                header.Date = sscanf(line, '%d %d %d')';
            case 7
                header.Numtrials = sscanf(line, '%d');
            case 8
                header.Rstim = sscanf(line, '%d');
            case 9
                header.Lstim = sscanf(line, '%d');
            case 10
                header.NGstim = sscanf(line, '%d');
        end
    end
    line = fgetl(fid);
    linenum = linenum + 1;
end
fclose(fid);