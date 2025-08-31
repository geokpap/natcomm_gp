function [histos, sessions, mode, ntrials] = dg_getTrackedUnit( ...
    unitcatalog, unitID, eventID)
%dg_getTrackedUnit returns PETHs for all sessions for a single unit
%	tracked across sessions.  It is assumed that the animal directories are
%	in the same directory as <unitcatalog>.
%INPUTS
%   unitcatalog: a string containing the pathname to an Excel workbook
%       containing the catalog of tracked units, one per worksheet.
%   unitID: a string that matches the tab on one of the worksheets.
%   eventID: numerical event ID to use as time zero in the PETH.
%OUTPUTS
%   histos: a 40-column array with one row for each session.  If there are
%       no events <eventID> in a particular session, then the entire row
%       contains NaN.  Only those trials that contain an event <eventID>
%       are counted, unless <eventID> is 1, in which case all trials are
%       counted and the reference time is the start of recording.
%   sessions: the session number for each row of <histos>.
%   mode: a character vector (i.e. a string) containing the "A/T" column 
%       for the unit.
%   ntrials: the number of trials that were counted for the row.  If there
%       are no events <eventID> in a particular session, then the value is
%       0.
%ERRORS
%   "Specified worksheet was not found."  <unitID> must be an exact match
% for the tab label on one of the worksheets.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

histos = [];
mode = '';
[num, txt] = xlsread(unitcatalog, unitID);
histos = NaN(size(num,1), 40);
ntrials = zeros(size(num,1), 1);
numsessions = size(num,1);
if numsessions == 0
    warning('dg_getTrackedUnit:nosessions', ...
        'The sheet for %s contains no sessions.', ...
        unitID );
    return
end
sessions = num(:,1);
mode = upper(char(txt(4:end,end)));
dataroot = fileparts(unitcatalog);
animalID = txt{4,1};
tetrode = num(1,2);
binwidth = 200; % 0.1 ms units
binedges = (-20:20) * binwidth;
for sessionidx = 1:numsessions
    tfilename = fullfile(dataroot, animalID, ...
        sprintf('ACQ%02d.TT%d', num(sessionidx, 1), tetrode) );
    [FileHeader, TrialData] = dg_ReadRodentFormat(tfilename);
    clustnum = num(sessionidx, 3);
    spiketimes = NaN(FileHeader.SSize, 1);
    spikeidx = 1;
    for trial = 1:length(TrialData)
        reftime = TrialData(trial).events(eventID);
        if reftime ~= 0 || eventID == 1
            ntrials(sessionidx) = ntrials(sessionidx) + 1;
            trialspiketimes = TrialData(trial).spikes(:,clustnum);
            trialspiketimes(trialspiketimes==0) = [];
            spiketimes( ...
                (1:TrialData(trial).header.SSize(clustnum)) + spikeidx - 1 ) ...
                = trialspiketimes - reftime;
            spikeidx = spikeidx + TrialData(trial).header.SSize(clustnum);
        end
    end
    if ntrials(sessionidx)
        counts = histc(spiketimes(~isnan(spiketimes)), binedges);
        histos(sessionidx,:) = counts(1:end-1);
    end
end

        