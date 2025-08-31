function lfp_seq2File(seqs, varargin)
%LFP_SEQ2FILE writes one or more sequences to a text spreadsheet file.
%lfp_seq2File(seqs)
%lfp_seq2File(..., 'headers', headers)
%lfp_seq2File(..., 'preset', filename)

% If <seqs> is a cell vector, it is a set of fixation sequences as returned
% by lfp_TargSeq; if it is a 2D cell array, then each column is a set of
% fixation sequences. Creates a tab-delimited spreadsheet with Unique Trial
% ID in the first column, and for each set of sequences, a column
% containing one sequence per row formatted as a string.
% OPTIONS
% 'preset' - bypasses GUI for selecting output filename and uses <filename>
%   instead; note that if <filename> is not an absolute pathname, then the
%   output file will be written to the current working directory.
% 'headers' - column header strings supplied in the cell vector <headers>.
%   Default is to write data only, without headers.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

headersflag = false;
presetflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'headers'
            headersflag = true;
            headers = varargin{argnum + 1};
            argnum = argnum + 1;
        case 'preset'
            presetflag = true;
            seqOutputFileName = varargin{argnum + 1};
            OutputPathName = '';
            argnum = argnum + 1;
        otherwise
            error('lfp_seq2File:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% Reconcile numbers of trials, get number of sets
trials = lfp_enabledTrials(1:length(lfp_SelectedTrials));
if numel(seqs) == length(seqs)
    % vector
    numtrials = length(seqs);
    numsets = 1;    % number of sets of sequences
    seqs = reshape(seqs,[],1);  % put in column orientation
else
    % 2D
    numtrials = size(seqs,1);
    numsets = size(seqs,2);
end
if length(trials) ~= numtrials
    error('lfp_seqTabulation:badSelect', ...
        'Trial selection state has changed since generating <fixations>' );
end

% Reconcile headers
if headersflag
    if length(headers) ~= size(seqs, 2)
        error('lfp_seq2File:badargs1', ...
            'The number of headers must match the number of sequences.');
    end
end

% Get output filename from user
if ~presetflag
    seqOutputFileName = 'sequences.xls';
    [seqOutputFileName, OutputPathName] = uiputfile( ...
        fullfile(lfp_DataDir, seqOutputFileName), ...
        'Save fixation table as:' );
    if isequal(seqOutputFileName, 0)
        return
    end
end

seqoutfid = fopen(fullfile(OutputPathName, seqOutputFileName), 'w');
if seqoutfid == -1
    error('lfp_seq2File:noOutfid1', ...
        'Could not open sequences output file');
end

% Print headers
if headersflag
    fprintf(seqoutfid, 'TrialID');
    for k = 1:length(headers)
        fprintf(seqoutfid, '\t%s', headers{k});
    end
    fprintf(seqoutfid, '\n');
end

for trialidx = 1:size(seqs,1)
    fprintf(seqoutfid, '%s', ...
        lfp_getTrialID(trials(trialidx)) );
    for setnum = 1:numsets
        seqstr = '';
        if ~isempty(seqs{trialidx, setnum})
            seqstr = sprintf('%02d', seqs{trialidx, setnum}(1,2));
            if size(seqs{trialidx, setnum},1) > 1
                seqstr = ...
                    [ seqstr sprintf(' %02d', seqs{trialidx, setnum}(2:end,2)) ];
            end
        end
        fprintf(seqoutfid, '\t%s', ...
            seqstr );
    end % for setnum
    fprintf(seqoutfid, '\n');
end % for trialidx

fclose(seqoutfid);
