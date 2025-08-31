function dg_opcham2nev(opchamfile, nlxfile)
% Converts a "sessiondata" file as saved by operant chamber code to a
% Neuralynx events (.NEV) file. Heedlessly blows away any file that may
% already exist at <nlxfile>.
%NOTES
% The "sessiondata" file should contain a variable that is a scalar struct
% with fields 'sessionDir', 'sessionName', and 'trials'. The first two
% fields are strings.  The last is a struct vector, one element per trial,
% with fields 'successful' (a Boolean scalar), and 'events'.  The 'events'
% field is a cell vector of struct arrays, but in my sample file each cell
% just contains a single element with fields:
%   TimeStampArray *
%   EventIDArray
%   TtlValueArray *
%   EventStringArray *
%   NumRecordsReturned
%   NumRecordsDropped
%   EventName
% Asterisks identify the fields actually used in this code. This function
% issues a warning if there is more than one value in any cell, and
% processes all values as if each were a distinct event, but I don't know
% if that's the correct interpretation.
%   This code can be made to run faster.  It is currently written for
% simplicity, not speed.

%$Rev: 161 $
%$Date: 2012-11-14 17:33:14 -0500 (Wed, 14 Nov 2012) $
%$Author: dgibson $

S = load(opchamfile);
fields = fieldnames(S);
sessiondata = S.(fields{1});
TS = [];
TTL = [];
ES = {};
Hdr = {
    '######## Neuralynx Data File Header'
    sprintf('## File Name: %s ', ...
    fullfile(sessiondata.sessionDir, sessiondata.sessionName))
    sprintf('## Time Opened: (m/d/y): %s  At Time: %s ', ...
    datestr(now, 2), datestr(now, 2))
    '## Converted by dg_opcham2nev '
    ' '
    };

for trialidx = 1:length(sessiondata.trials)
    for evtidx = 1:length(sessiondata.trials(trialidx).events)
        evts = sessiondata.trials(trialidx).events{evtidx};
        if length(evts) > 1
            warning('dg_opcham2nev:evts', ...
                'Trial %d event %d contains more than one value', ...
                trialidx, evtidx);
        end
        for tsidx = 1:length(evts.TimeStampArray)
            TS(end+1) = evts.TimeStampArray(tsidx); %#ok<*AGROW>
            if isempty(evts.TtlValueArray)
                TTL(end+1) = 0;
            else
                TTL(end+1) = evts.TtlValueArray(tsidx);
            end
            if isempty(evts.EventStringArray)
                ES{end+1} = '';
            else
                ES{end+1} = evts.EventStringArray{tsidx};
            end
        end
    end
end

dg_writeNlxEvents(nlxfile, TS, TTL, ES, Hdr);

