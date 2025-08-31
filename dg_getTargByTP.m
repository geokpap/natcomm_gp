function [seqs, numbers] = dg_getTargByTP( ...
    DSN, paramnames, targtable, targname )
%DG_GETTARGBYTP gets subtotals grouped by Trial Param from an ODBC
%database.
%[seqs, numbers] = dg_getTargByTP(DSN, paramnames, targtable, targname)
%INPUTS:
% DSN: string containing name of ODBC data source
% paramnames: cell string vector containing column names to be used for
%   grouping result tabulation
% targtable: string containing name of target table
% targname: string containing name of target column
%OUTPUTS & DESCRIPTION:
% The database pointed to by <DSN> must contain a table named 'params' with
% a unique key column named 'TrialID' and columns named <paramnames>, and
% another table named <targtable> with the same 'TrialID' column and a
% column named <targname>.  Returns <seqs> and <numbers> in the same format
% as dg_readSeqByTP:  <seqs> is a cell string column vector containing the
% values of <targname> that are unique within each unique combination of
% non-zero values of the <paramnames>.  <numbers> contains one row
% corresponding to each row of <seqs>, with columns that contain the values
% of the <paramnames> and a final column containing the count of all
% records having the unique combination of values in <seqs> and
% <paramnames>.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

colnames = '';
for k = 1:length(paramnames)
    colnames = [ colnames paramnames{k} ', '];
end
colnames = [colnames targname];
colconditions = 

logintimeout(5);
conn = database(DSN, '', '');
setdbprefs('DataReturnFormat', 'cellarray', 'NullStringRead', '');

sqlstr = sprintf(...
    ['SELECT %s, COUNT(*) FROM %s t, params p ' ...
    'WHERE t.TrialID = p.TrialID AND %s GROUP BY %s'], ...
    colnames, targtable, colconditions, colnames );
