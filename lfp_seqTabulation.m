function lfp_SeqTabulation(seqs, varargin)
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
% "renumno0", which is "renumseq" with 0's deleted).  Produces those same
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
% 'renumno0' - add the "renumno0" and "renumno0to9" columns to the output.
% 'headers' - overrides the default column header strings with whatever
%   values are supplied in the cell vector <headers>.  This applies only to
%   the sequence columns; the first column header remains "TrialID".  There
%   must be exactly twice as many headers as sequences if 'renumno0' is not
%   given, and exactly three times as many if it is.
% 'db', <DSN>, <targtable> - send data directly to an ODBC database instead
%   of a text file. <DSN> is a string containing the Windows Data Source
%   Name for the database, <targtable> is a string containing the name of
%   the target table in the database to which the data should be added.
%   The target table must already contain an entry for each TrialID, and
%   must already contain a column for each name in the header row that
%   would normally appear in the output text file.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

headers = {'rawseq', 'renumseq', 'renumno0', 'renumno0to9'};
dbflag = false;
headersflag = false;
presetflag = false;
renumno0flag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'db'
            if length(varargin) < argnum+2
                error('makeOrderedHist:badargs', ...
                    'The ''db'' option requires 2 additional arguments' );
            end
            dbflag = true;
            DSN = varargin{argnum+1};
            targtable = varargin{argnum+2};
            argnum = argnum + 2;
        case 'headers'
            headersflag = true;
            headers = varargin{argnum + 1};
            argnum = argnum + 1;
        case 'preset'
            presetflag = true;
            seqOutputFileName = varargin{argnum + 1};
            OutputPathName = '';
            argnum = argnum + 1;
        case 'renumno0'
            renumno0flag = true;
        otherwise
            error('lfp_SeqTabulation:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~renumno0flag && ~headersflag
    headers(end) = [];
end

if headersflag
    if (renumno0flag && length(headers) ~= 4 * size(seqs, 2)) ...
            || (~renumno0flag && length(headers) ~= 2 * size(seqs, 2))
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
if dbflag
    logintimeout(5);
    conn = database(DSN, '', '');
    setdbprefs('DataReturnFormat', 'cellarray');
    setdbprefs('NullStringRead', '');
else
    seqoutfid = fopen(fullfile(OutputPathName, seqOutputFileName), 'w');
    if seqoutfid == -1
        error('lfp_SeqTabulation:noOutfid1', ...
            'Could not open sequences output file');
    end
end

for setnum = 1:numsets
    s2(:,setnum) = lfp_deleteRepeatTargs(seqs(:,setnum));
    s3(:,setnum) = lfp_RenumberTargSeq(s2(:,setnum), 10);
    if renumno0flag
        s4(:,setnum) = lfp_SeqSubst(s3(:,setnum), 0, []);
        % added here only111213 seq
        s5(:,setnum) = lfp_SeqDelete(s4(:,setnum), [1:1:9]);
        % eliminate consecutive repeats
        s6(:,setnum) = lfp_deleteRepeatTargs(s5(:,setnum));
    end
end

% Print headers
if ~dbflag
    fprintf(seqoutfid, 'TrialID');
    for k = 1:length(headers)
        fprintf(seqoutfid, '\t%s', headers{k});
    end
    fprintf(seqoutfid, '\n');
end

for trialidx = 1:length(trials)
    colidx = 1; % needed for 'db' option
    if ~dbflag
        fprintf(seqoutfid, '%s', ...
            lfp_getTrialID(trials(trialidx)) );
    end
    for setnum = 1:numsets
        % Create sequence strings
        rawseqstr = '';
        renumseqstr = '';
        renumno0str = '';
        renumno0to9str = '';
        if ~isempty(s2{trialidx, setnum})
            rawseqstr = sprintf('%02d', s2{trialidx, setnum}(1,2));
            if size(s2{trialidx, setnum},1) > 1
                rawseqstr = ...
                    [ rawseqstr sprintf(' %02d', s2{trialidx, setnum}(2:end,2)) ];
            end
        end
        if ~isempty(s3{trialidx, setnum})
            renumseqstr = sprintf('%02d', s3{trialidx, setnum}(1,2));
            if size(s3{trialidx, setnum},1) > 1
                renumseqstr = ...
                    [ renumseqstr sprintf(' %02d', s3{trialidx, setnum}(2:end,2)) ];
            end
        end
        if renumno0flag && ~isempty(s4{trialidx, setnum})
            renumno0str = sprintf('%02d', s4{trialidx, setnum}(1,2));
            if size(s4{trialidx, setnum},1) > 1
                renumno0str = ...
                    [ renumno0str sprintf(' %02d', s4{trialidx, setnum}(2:end,2)) ];
            end
            if ~isempty(s6{trialidx, setnum})
                renumno0to9str = sprintf('%02d', s6{trialidx, setnum}(1,2));
                if size(s6{trialidx, setnum},1) > 1
                    renumno0to9str = ...
                        [ renumno0to9str sprintf(' %02d', s6{trialidx, setnum}(2:end,2)) ];
                end
            end
        end
        % Output the sequence strings
        if dbflag
            colvalpairs = sprintf('%s %s=''%s'', %s=''%s''', ...
                colvalpairs, headers{colidx}, rawseqstr, ...
                headers{colidx+1}, renumseqstr );
            if renumno0flag
                colvalpairs = sprintf(', %s %s=''%s'', %s=''%s''', ...
                    colvalpairs, headers{colidx+2}, renumno0str, ...
                    headers{colidx+3}, renumno0to9str );
                colidx = colidx + 4;
            else
                colidx = colidx + 2;
            end
            if colidx <= length(headers)
                colvalpairs = [ colvalpairs ', ' ];
            end
        else
            fprintf(seqoutfid, '\t%s\t%s', ...
                rawseqstr, renumseqstr);
            if renumno0flag
                fprintf(seqoutfid, '\t%s', renumno0str);
                fprintf(seqoutfid, '\t%s', renumno0to9str);
            end
        end
    end % for setnum
    if dbflag
        sqlstr = sprintf('UPDATE %s SET %s WHERE TrialID = ''%s''', ...
            targtable, colvalpairs, lfp_getTrialID(trials(trialidx)));
        curs = exec(conn, sqlstr);
        if ~isempty(curs.Message)
            close(conn);
            error('lfp_SeqTabulation:badSQL', ...
                'The resulting SQL statement\n%s\ncould not be executed by database %s', ...
                sqlstr, DSN );
        end    else
        fprintf(seqoutfid, '\n');
    end
end % for trialidx

fclose(seqoutfid);