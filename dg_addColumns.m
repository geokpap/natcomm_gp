function dg_addColumns(DSN, tablename, filename, datatype)
%dg_addColumns(DSN, filename, datatype)
% Opens an ODBC connection to <DSN> and checks each column header in
% <filename> to see if it already exists in the table <tablename> in <DSN>.
% If the column doesn't exist, a new column is created of type <datatype>.
% If the table doesn't exist, nothing happens.
%INPUTS
% DSN: a string that is an ODBC Data Source Name that requires no user name
%   or password.
% tablename: a string that is the name of a table in <DSN>.
% filename: a string that is the pathname (absolute or relative) of a
%   tab-delimited text file whose first row is a header (i.e. contains
%   column names).
% datatype: a string that denotes a valid data type for <DSN>.
%
%NOTES
% We rely on the fact that in both Microsoft Access and MySQL, when an
% attempt is made to add a column name that already exists in the specified
% table, it silently fails regardless of the data type specified
% (curs.message is "No ResultSet was produced" in case of success for both
% DBs, but the error messages differ).

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

% Get column names
fid = fopen(filename);
if isequal(fid, -1)
    error('dg_addColumns:badinfile1', ...
        'Could not open input file "%s"', filename );
end
header = fgetl(fid);
fclose(fid);
if isequal(header, -1) || isempty(header)
    error('dg_addColumns:badinfile2', ...
        'File "%s" does not contain a complete header line', filename );
end
tabidx = regexp(header, '\t');

% Add column names as needed
logintimeout(5);
conn = database(DSN, '', '');
setdbprefs('DataReturnFormat','cellarray');
for colnum = 1 : (length(tabidx) + 1)
    if colnum > 1
        startidx = tabidx(colnum - 1) + 1;
    else
        startidx = 1;
    end
    if colnum <= length(tabidx)
        endidx = tabidx(colnum) - 1;
    else
        endidx = length(header);
    end
    sqlstring = sprintf('ALTER TABLE %s ADD %s %s', ...
        tablename, header(startidx:endidx), datatype );
    curs = exec(conn, sqlstring);
end
close(conn);
