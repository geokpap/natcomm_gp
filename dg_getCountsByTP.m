function [seqs, numbers] = dg_getCountsByTP( ...
    DSN, paramnames, targtable, targname, select_cond )
%DG_GETCOUNTSBYTP gets trial counts grouped by Trial Param from an ODBC
%database.
%[seqs, numbers] = [seqs, numbers] = dg_getCountsByTP( ...
%    DSN, paramnames, targtable, targname, select_cond )
%INPUTS:
% DSN: string containing name of ODBC data source
% paramnames: cell string vector containing column names to be used for
%   grouping result tabulation
% targtable: string containing name of target table
% targname: string containing name of target column
% select_cond: the predicate of the WHERE clause in a SQL SELECT statement,
%   e.g. 'p17 = 100 AND p40 IS NULL'; may be an empty string (''). 
%OUTPUTS & DESCRIPTION:
% The database pointed to by <DSN> must contain a table named 'params' with
% a text unique key column named 'TrialID' and numeric columns named
% <paramnames>, and another table named <targtable> with the same 'TrialID'
% column and a text column named <targname>.  Returns <seqs> and <numbers>
% in the same format as dg_readSeqByTP:  <seqs> is a cell string column
% vector containing the values of <targname> that are unique within each
% unique combination of values of the <paramnames>.  <numbers>
% contains one row corresponding to each row of <seqs>, with columns that
% contain the values of the <paramnames> and a final column containing the
% count of all records having the unique combination of values in <seqs>
% and <paramnames>.  (<seqs> could contain values from any text column in
% the database, not necessarily only sequences.)

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

colnames = '';
if isempty(select_cond)
    colconditions = '';
else
    colconditions = sprintf('AND (%s)', select_cond);
end
for k = 1:length(paramnames)
    colnames = [ colnames paramnames{k} ', '];
end
colnames = [colnames targname];
sqlstr = sprintf(...
    ['SELECT %s, COUNT(*) FROM %s t, params p ' ...
    'WHERE t.TrialID = p.TrialID %s GROUP BY %s'], ...
    colnames, targtable, colconditions, colnames );

logintimeout(5);
conn = database(DSN, '', '');
setdbprefs('DataReturnFormat', 'cellarray');
setdbprefs('NullStringRead', '');
curs = exec(conn, sqlstr);
if ~isempty(curs.Message)
    close(conn);
    error('dg_getCountsByTP:badquery', ...
        'The resulting query\n%s\ncould not be answered by database %s.\n"%s"', ...
        sqlstr, DSN, curs.Message );
end
curs = fetch(curs);
seqs = curs.Data(:, end-1);
% skip second to last column, which was <seqs>:
numbers = [ cell2mat(curs.Data(:,1:end-2)) ...
    cell2mat(curs.Data(:,end)) ];
end
close(conn);
