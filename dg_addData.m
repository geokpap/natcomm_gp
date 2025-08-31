function dg_addData(DSN, tablename, filename, datatype)
%dg_addColumns(DSN, filename, datatype)
% Opens an ODBC connection to <DSN> and checks each record in the file
% <filename> to see if it already exists in the table <tablename> in <DSN>.
% The first column in <filename> is assumed to be the primary key, of a
% text data type. If the record doesn't exist, a new record is inserted
% containing data for all the columns that exist in <filename>.  If the
% record does exist, then the values for all the columns are updated. If
% the table doesn't exist, nothing happens.
%INPUTS
% DSN: a string that is an ODBC Data Source Name that requires no user name
%   or password.
% tablename: a string that is the name of a table in <DSN>.
% filename: a string that is the pathname (absolute or relative) of a
%   tab-delimited text file whose first row is a header (i.e. contains
%   column names).
% datatype: either 'text' (case-insensitive) or anything else; data are
%   interpreted as numeric if <datatype> is not 'text'.
%
%NOTES
% We rely on the fact that in both Microsoft Access and MySQL, when an
% attempt is made to add a column name that already exists in the specified
% table, it silently fails regardless of the data type specified
% (curs.message is "No ResultSet was produced" in case of success for both
% DBs, but the error messages differ).
%
% We use the "fastinsert" function, so this code requires Matlab 7.1 or
% higher.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

textdata = isequal(upper(datatype), 'TEXT');

% First, find the line terminator string
fid = fopen(filename);
if isequal(fid, -1)
    error('dg_addColumns:badinfile1', ...
        'Could not open input file "%s"', filename );
end
lineWithEOL = fgets(fid);
fclose(fid);
fid = fopen(filename);
line = fgetl(fid);
fclose(fid);
EOL = lineWithEOL(length(line)+1:end);
EOLlen = length(EOL);

% Slurp the whole file and find the EOLs (yes, it actually is 4 times
% slower to call fgetl repeatedly):
fid = fopen(filename);
filechars = reshape(fread(fid, Inf, 'uchar'), 1, []);
fclose(fid);
EOLidx = strfind(filechars, EOL);
% Make sure there is a trailing EOLidx, so that there will be exactly one
% EOLidx per line:
if ~isequal(filechars(end-length(EOL)+1:end), EOL)
    EOLidx(end+1) = length(filechars+1);
end

% Get key and data column names
line = char(filechars(1:EOLidx(1)-1));
tabidx{1} = regexp(line, '\t');
keyname = line(1:tabidx{1}(1)-1);
numcols = length(tabidx{1});
for colnum = 1:numcols
    if colnum == length(tabidx{1})
        colnames{colnum} = line(tabidx{1}(colnum) + 1 : end);
    else
        colnames{colnum} = line(tabidx{1}(colnum) + 1 : ...
            tabidx{1}(colnum+1) - 1) ;
    end
end

logintimeout(5);
conn = database(DSN, '', '');
if isequal(conn.Handle, 0)
    error('dg_addData:nodb', ...
        'Could not open database %s', DSN );
end
setdbprefs('DataReturnFormat','numeric');

% Find which key values already exist, to sort lines into inserts and
% updates:
hWaitBar = waitbar(0, 'Creating list of keys', ...
    'Name', 'Searching for pre-existing records');
for linenum = 2:length(EOLidx)
    line = char(filechars(EOLidx(linenum-1)+EOLlen : EOLidx(linenum)-1));
    tabidx{linenum} = regexp(line, '\t');
    if length(tabidx{linenum}) ~= numcols
        % missing or extra columns will be ignored
        warning('dg_addData:badline1', ...
            'Line %d contains %d columns instead of %d.', ...
            linenum, length(tabidx{linenum}), numcols );
    end
    keyval{linenum, 1} = line(1:tabidx{linenum}(1)-1);
    % This is abominably slow:
%     sqlstr = sprintf('SELECT COUNT(*) FROM %s WHERE %s = ''%s''', ...
%         tablename, keyname, keyval{linenum} );
%     curs = exec(conn, sqlstr);
%     curs = fetch(curs);
%     recExists(linenum) = curs.data;
    waitbar(linenum/length(EOLidx), hWaitBar);
end
waitbar(1, hWaitBar, 'Querying database');
% The 'TEMPORARY' keyword does not work in Access 9.0.3821:
success = false;
while ~success
    temptablename = sprintf('dg_addData%015.0f', fix(1e15*rand));
    curs = exec(conn, sprintf( ...
        'CREATE TABLE %s (keyfield TEXT NOT NULL, linenum INT)', temptablename ));
    success = isequal(curs.Message, 'No ResultSet was produced');
end
fastinsert(conn, temptablename, {'keyfield', 'linenum'}, ...
    [ keyval(2:end) ...
    mat2cell((2:length(EOLidx))', ones(size(keyval(2:end))), 1) ] );
sqlstr = sprintf('SELECT linenum FROM %s a, %s b WHERE a.%s = b.keyfield', ...
    tablename, temptablename, keyname );
curs = exec(conn, sqlstr);
curs = fetch(curs);
if ~isequal(curs.data, {'No Data'})
    lines2update = (curs.data)';
else
    lines2update = [];
end
curs = exec(conn, sprintf('DROP TABLE %s', temptablename ));
close(hWaitBar);

% Do the inserts
lines2insert = setdiff(2:length(EOLidx), lines2update);
insertdata = cell(length(lines2insert), numcols + 1);
insertdata(:,2:end) = {'null'};
hWaitBar = waitbar(0, 'Preparing records', ...
    'Name', 'Inserting new records');
for recnum = 1:length(lines2insert)
    linenum = lines2insert(recnum);
    line = char(filechars(EOLidx(linenum-1)+EOLlen : EOLidx(linenum)-1));
    % Put key values in insertdata column 1:
    insertdata{recnum, 1} = keyval{linenum};
    % Put data values in insertdata columns 2 - numcols+1:
    for colnum = 1:length(tabidx{linenum})
        if colnum == length(tabidx{linenum})
            insertdata{recnum, colnum + 1} = ...
                line(tabidx{linenum}(colnum) + 1 : end) ;
        else
            insertdata{recnum, colnum + 1} = ...
                line(tabidx{linenum}(colnum) + 1 : ...
                tabidx{linenum}(colnum+1) - 1) ;
        end
        if isempty(insertdata{recnum, colnum + 1})
            insertdata{recnum, colnum + 1} = 'null';
        else
            if ~textdata
                insertdata{recnum, colnum + 1} = ...
                    str2num(insertdata{recnum, colnum + 1});
            end
        end
    end
    waitbar(recnum/length(lines2insert), hWaitBar);
end
waitbar(1, hWaitBar, 'Performing insertion');
fastinsert(conn, tablename, [{keyname} colnames], insertdata);
close(hWaitBar);

% Select formatting strings for updating existing records:
if textdata
    colfmt = '%s %s=''%s'',';
    lastcolfmt = '%s %s=''%s''';
else
    colfmt = '%s %s=%s,';
    lastcolfmt = '%s %s=%s';
end

% Do the updates
hWaitBar = waitbar(0, '', ...
    'Name', 'Updating existing records');
for linenum = lines2update
    line = char(filechars(EOLidx(linenum-1)+EOLlen : EOLidx(linenum)-1));
    sqlstr = sprintf('UPDATE %s SET ', tablename);
    for colnum = 1:length(tabidx{linenum})
        if colnum == length(tabidx{linenum})
            sqlstr = sprintf(lastcolfmt, ...
                sqlstr, colnames{colnum}, ...
                line(tabidx{linenum}(colnum) + 1 : end) );
        else
            sqlstr = sprintf(colfmt, ...
                sqlstr, colnames{colnum}, ...
                line(tabidx{linenum}(colnum) + 1 : ...
                tabidx{linenum}(colnum+1) - 1) );
        end
    end
    sqlstr = sprintf('%s WHERE %s = ''%s''', sqlstr, keyname, keyval{linenum});
    curs = exec(conn, sqlstr);
    waitbar(linenum/length(EOLidx), hWaitBar);
end
close(hWaitBar);
close(conn);
