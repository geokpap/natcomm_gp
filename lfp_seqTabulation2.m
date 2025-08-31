function lfp_SeqTabulation2(seqs, varargin)
%lfp_SeqTabulation(seqs)
%lfp_SeqTabulation(..., 'headers', headers)
%lfp_SeqTabulation(..., 'preset', filename)

% If <seqs> is a cell vector, it is a set of fixation sequences as returned
% by lfp_TargSeq; if it is a 2D cell array, then each column is a set of
% fixation sequences, and the 'headers' option must be supplied.
% For a single set of fixation sequences, creates a tab-delimited
% spreadsheet with 3 columns: "TrialID" (Unique Trial ID), "rawseq" (raw
% sequence string formatted from <seqs> with repeated targets deleted),
% "renumseq" ("rawseq" with colored targets renumbered); may optionally add
% "1stonlySL" and "1stonlyRO").  Produces those same
% output columns for each set of seqs when <seqs> contains multiple sets,
% except that the "TrialID" column is not duplicated.
%
% Note that lfp_SelectedTrials and lfp_BadTrials must be unchanged since
% <seqs> was generated.  Unfortunately, all that we can do to enforce that
% condition is check that there is still the same number of trials
% selected, which is necessary but not sufficient.
%
% OPTIONS
% 'preset' - bypasses GUI for selecting output filename and uses <filename>
%   instead; note that if <filename> is not an absolute pathname, then the
%   output file will be written to the current working directory.
% '1stonly' - add the "1stonlySL" and "1stonlyRO" columns to the output.
% 'headers' - overrides the default column header strings with whatever
%   values are supplied in the cell vector <headers>.  This applies only to
%   the sequence columns; the first column header remains "TrialID".  There
%   must be exactly twice as many headers as sequences if '1stonly' is not
%   given, and exactly three times as many if it is.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

headers = {'rawseq', 'renumseq', '1stonlySL', '1stonlyRO'};
headersflag = false;
presetflag = false;
firstonlyflag = false;
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
        case '1stonly'
            firstonlyflag = true;
        otherwise
            error('lfp_SeqTabulation:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~firstonlyflag && ~headersflag
    headers(end) = [];
end

if headersflag
    if (firstonlyflag && length(headers) ~= 2 * size(seqs, 2)) ...
            || (~firstonlyflag && length(headers) ~= 2 * size(seqs, 2))
        error('lfp_SeqTabulation:badargs1', ...
            'The number of headers must match the number of sequences.');
    end
else
    if (numel(seqs) ~= length(seqs))
        error('lfp_SeqTabulation:badargs2', ...
            'You must supply ''headers'' when tabulating multiple sequences.');
    end
end

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
    error('lfp_SeqTabulation:badSelect', ...
        'Trial selection state has changed since generating <fixations>' );
end

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
    error('lfp_SeqTabulation:noOutfid1', ...
        'Could not open sequences output file');
end

for setnum = 1:numsets
    s2(:,setnum) = lfp_deleteRepeatTargs(seqs(:,setnum));
    s3(:,setnum) = lfp_RenumberTargSeq(s2(:,setnum), 10);
    if firstonlyflag
        s4(:,setnum) = lfp_SeqSubst(s3(:,setnum), 0, []);
        % seq of only 11, 12, 13 (1st, 2nd, 3rd targets)
        s5(:,setnum) = lfp_SeqDelete(s4(:,setnum), [1:1:9]);
        % eliminate consecutive repeats
        s6(:,setnum) = lfp_deleteRepeatTargs(s5(:,setnum));
        % reduce sequence to just the first occurrence of each target
        s7(:,setnum) = lfp_seqFirstOnly(s6(:,setnum));
        % S8 is s7 renumbered in terms of absolute spatial position of each
        % target
        s8(:,setnum) = lfp_RecoverTargSeq(s7(:,setnum));
    end
end

% Print headers
fprintf(seqoutfid, 'TrialID');
for k = 1:length(headers)
    fprintf(seqoutfid, '\t%s', headers{k});
end
fprintf(seqoutfid, '\n');

for trialidx = 1:length(trials)
    fprintf(seqoutfid, '%s', ...
        lfp_getTrialID(trials(trialidx)) );
    for setnum = 1:numsets
        rawseqstr = '';
        renumseqstr = '';
        firstonlySLstr = '';
        firstonlyROstr = '';
%         if ~isempty(s2{trialidx, setnum})
%             rawseqstr = sprintf('%02d', s2{trialidx, setnum}(1,2));
%             if size(s2{trialidx, setnum},1) > 1
%                 rawseqstr = ...
%                     [ rawseqstr sprintf(' %02d', s2{trialidx, setnum}(2:end,2)) ];
%             end
%         end
%         if ~isempty(s3{trialidx, setnum})
%             renumseqstr = sprintf('%02d', s3{trialidx, setnum}(1,2));
%             if size(s3{trialidx, setnum},1) > 1
%                 renumseqstr = ...
%                     [ renumseqstr sprintf(' %02d', s3{trialidx, setnum}(2:end,2)) ];
%             end
%         end
        if firstonlyflag && ~isempty(s8{trialidx, setnum})
            firstonlySLstr = sprintf('%02d', s8{trialidx, setnum}(1,2));
            if size(s8{trialidx, setnum},1) > 1
                firstonlySLstr = ...
                    [ firstonlySLstr sprintf(' %02d', s8{trialidx, setnum}(2:end,2)) ];
            end
            if ~isempty(s7{trialidx, setnum})
                firstonlyROstr = sprintf('%02d', s7{trialidx, setnum}(1,2));
                if size(s7{trialidx, setnum},1) > 1
                    firstonlyROstr = ...
                        [ firstonlyROstr sprintf(' %02d', s7{trialidx, setnum}(2:end,2)) ];
                end
            end
        end
%         fprintf(seqoutfid, '\t%s\t%s', ...
%             rawseqstr, renumseqstr);
        if firstonlyflag
            fprintf(seqoutfid, '\t%s', firstonlySLstr);
            fprintf(seqoutfid, '\t%s', firstonlyROstr);
        end
    end % for setnum
    fprintf(seqoutfid, '\n');
end % for trialidx

fclose(seqoutfid);