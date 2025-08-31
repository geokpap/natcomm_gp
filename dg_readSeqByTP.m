function [seqs, numbers] = dg_readSeqByTP(filename)
%DG_READSEQBYTP reads grouped subtotals exported from FileMaker.
%[seqs, numbers] = dg_readSeqByTP(filename)
% The file pointed to by <filename> should be a tab-delimited ASCII
% spreadsheet containing variable-length string data in the first column,
% used as the innermost sort field in FileMaker, additional grouping
% parameters in columns. 2 - 4 (aka cols. B - D), and summary data such as
% counts in column 5 (E). <seqs> is a cell string column vector containing
% column 1 of <filename>. <numbers> is a numeric array with as many rows as
% the spreadsheet, and one column for each of spreadsheet columns 2 - 5 (B
% - E). Where the spreadsheet contains entries in columns 2 - 4 (B - D),
% the elements in <numbers> are those entries.  Where the spreadsheet
% contains empty cells, the value from the previous row is propagated, so
% the effect is to "Fill Down" each entry to fill all the empty cells until
% the next entry.  It is assumed that every cell in column 5 contains a
% value and that every cell on the first row contains a value.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

% Pre-allocate storage:
BLOCKSIZE = 1024;
NUMCOLS = 4;
seqs = cell(2*BLOCKSIZE,1);
numbers = zeros(2*BLOCKSIZE,NUMCOLS);

% Read and "fill-down" the file:
infid = fopen(filename);
if isequal(infid, -1)
    error('dg_readSeqByTP:badOutput', ...
        'Could not open input file %s', filename );
end
line = fgets(infid);
linenum = 1;
while ~isequal(line, -1)
    % allocate more memory if needed
    if linenum > size(numbers, 1)
        numbers = [ numbers; zeros(BLOCKSIZE,NUMCOLS) ];
        seqs = [ seqs; cell(BLOCKSIZE,1) ];
    end
    % parse input line
    tabs = regexp(line, '\t');
    seqs{linenum} = line(1 : tabs(1) - 1);
    for tabidx = 1 : (length(tabs) - 1)
        if tabs(tabidx+1)-1 >= tabs(tabidx)+1
            numbers(linenum,tabidx) = str2num(line(...
                tabs(tabidx)+1 : tabs(tabidx+1)-1 ));
        else
            numbers(linenum,tabidx) = numbers(linenum-1,tabidx);
        end
    end
    numbers(linenum, end) = str2num(line(tabs(end)+1:end));
    line = fgets(infid);
    linenum = linenum + 1;
end
fclose(infid);

% Trim off unused pre-allocated rows, if any.  <linenum> is now 1 more than
% the number of rows used.
if linenum <= size(numbers, 1)
    numbers(linenum:end,:) = [];
    seqs(linenum:end,:) = [];
end

