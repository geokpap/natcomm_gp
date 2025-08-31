function lfp_TrialParamsTabulation(varargin)
%lfp_TrialParamsTabulation
%lfp_TrialParamsTabulation('preset', filename)

% Creates a tab-delimited spreadsheet with a variable number of columns: "TrialID" (Unique
% Trial ID), and a variable number of columns with names of the form
% "p<n>", where <n> denotes an integer indicating the trial parameter
% number.  Only enough "p<n>" columns are created to hold the longest value
% in lfp_TrialParams.
% OPTIONS
% 'preset' - bypasses GUI for selecting output filename and uses <filename>
%   instead; note that if <filename> is not an absolute pathname, then the
%   output file will be written to the current working directory.

%$Rev: 45 $
%$Date: 2009-02-11 15:56:11 -0500 (Wed, 11 Feb 2009) $
%$Author: dgibson $

lfp_declareGlobals;

presetflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'preset'
            presetflag = true;
            tpOutputFileName = varargin{argnum + 1};
            OutputPathName = '';
            argnum = argnum + 1;
        otherwise
            error('lfp_SeqTabulation:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));

if ~presetflag
    tpOutputFileName = 'TrialParams.xls';
    [tpOutputFileName, OutputPathName] = uiputfile( ...
        fullfile(lfp_DataDir, tpOutputFileName), ...
        'Save fixation table as:' );
    if isequal(tpOutputFileName, 0)
        return
    end
end
tpoutfid = fopen(fullfile(OutputPathName, tpOutputFileName), 'w');
if tpoutfid == -1
    error('lfp_TrialParamsTabulation:noOutfid1', ...
        'Could not open trial parameters output file');
end

% Find longest set of trial params:
maxlen = 0;
for trial = trials
    if length(lfp_TrialParams{trial}) > maxlen
        maxlen = length(lfp_TrialParams{trial});
    end
end

fprintf(tpoutfid, ['TrialID' sprintf('\tp%d', 1:maxlen)]);
fprintf(tpoutfid, '\n');

for trial = trials
    fprintf(tpoutfid, lfp_getTrialID(trial));
    fprintf(tpoutfid, '\t%d', lfp_TrialParams{trial});
    fprintf(tpoutfid, '\n');
end

fclose(tpoutfid);
