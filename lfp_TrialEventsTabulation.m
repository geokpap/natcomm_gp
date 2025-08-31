function lfp_TrialEventsTabulation(ltet_eventnames, varargin)
%lfp_TrialEventsTabulation(ltet_eventnames)
%lfp_TrialEventsTabulation(..., 'preset', filename)

%INPUTS
% ltet_eventnames - a cell string vector containing names of events that
% are
%   defined in the current setup's lfp_getEvtIDs_<setupname> script.
%DESCRIPTION
% Creates a tab-delimited text file containing a table with Unique Trial ID
% in column 1, trial start and end timestamps in columns 2 and 3 (from
% events listed in lfp_TrialIndex), and timestamp of the first instance of
% each event in <ltet_eventnames> in the following columns.  The first row
% is a header row that lists 'TrialID', 'start', 'end', and the names from
% <ltet_eventnames>.  Trials are included according to the usual
% combination of their selection state with lfp_BadTrials.
%OPTIONS
% 'preset' - bypasses GUI for selecting output filename and uses
% <filename>
%   instead; note that if <filename> is not an absolute pathname, then
%   the output file will be written to the current working directory.

%NOTES
% All the variables local to this function have names beginning with
% 'ltet_' to prevent conflicts with names defined in
% lfp_getEvtIDs_<setupname>, except for <varargin>, which should not be
% addressed at all after the call to lfp_getEvtIDs.  As for conflicts with
% function names, we cross our fingers and hope that Matlab has gotten
% smart enough to figure out which is which.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;

ltet_presetflag = false;
ltet_argnum = 1;
while ltet_argnum <= length(varargin)
    switch varargin{ltet_argnum}
        case 'preset'
            ltet_presetflag = true;
            ltet_OutputFileName = varargin{ltet_argnum + 1};
            ltet_OutputPathName = '';
            ltet_argnum = ltet_argnum + 1;
        otherwise
            error('lfp_TrialEventsTabulation:badoption', ...
                ['The option "' varargin{ltet_argnum} '" is not recognized.'] );
    end
    ltet_argnum = ltet_argnum + 1;
end

lfp_getEvtIDs;

ltet_trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));

if ~ltet_presetflag
    ltet_OutputFileName = 'trialevents.xls';
    [ltet_OutputFileName, ltet_OutputPathName] = uiputfile( ...
        fullfile(lfp_DataDir, ltet_OutputFileName), ...
        'Save trial events table as:' );
    if isequal(ltet_OutputFileName, 0)
        return
    end
end
ltet_outfid = fopen(fullfile(ltet_OutputPathName, ltet_OutputFileName), 'w');
if ltet_outfid == -1
    error('lfp_TrialEventsTabulation:noOutfid1', ...
        'Could not open sequences output file');
end

% Print headers
fprintf(ltet_outfid, 'TrialID\tstart\tend');
for ltet_nameix = 1:length(ltet_eventnames)
    fprintf(ltet_outfid, '\t%s', ltet_eventnames{ltet_nameix});
end
fprintf(ltet_outfid, '\n');

for ltet_trialnum = ltet_trials
    ltet_startix = lfp_TrialIndex(ltet_trialnum,1);
    ltet_endix = lfp_TrialIndex(ltet_trialnum,2);
    fprintf(ltet_outfid, '%s\t%.6f\t%.6f', ...
        lfp_getTrialID(ltet_trialnum), ...
        lfp_Events(ltet_startix,1), lfp_Events(ltet_endix,1) );
    for ltet_nameix = 1:length(ltet_eventnames)
        ltet_evtix = find(lfp_Events(ltet_startix+1:ltet_endix-1, 2) ...
            == eval(ltet_eventnames{ltet_nameix}) );
        if isempty(ltet_evtix)
            fprintf(ltet_outfid, '\t');
        else
            ltet_evtix = ltet_evtix(1) + ltet_startix;
            fprintf(ltet_outfid, '\t%.6f', lfp_Events(ltet_evtix, 1));
        end
    end
    fprintf(ltet_outfid, '\n');
end
fclose(ltet_outfid);
