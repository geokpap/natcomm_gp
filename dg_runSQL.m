function [data msg1 msg2] = dg_runSQL(datasrc, sqlstr)
%[data msg1 msg2] = dg_runSQL(datasrc, sqlstr)
% Executes the SQL statement <sqlstr> on database <datasrc> and then does a
% 'fetch'.  Empirically, it seems that NULL is represented for numeric
% fields as NaN, and for text fields as 'null' in <data>.  In either case,
% the data cell will NOT be empty.  This unfortunately makes it impossible
% to distinguish the legitimate non-NULL string value 'null' from the empty
% value NULL.  (The empty string is also a legitimate non-NULL value.)
%INPUTS
% datasrc - may be either a connection object or a string containing an ODBC
%   Data Source Name.  In the latter case, a new connection is opened
%   and then closed again.
% sqlstr - a string containing an executable SQL statement.
%OUTPUTS
% data - a cell array containing the data returned by the 'fetch'.
% msg1 - the message returned by the SQL execution; empty if no errors.
% msg2 - the message returned by the 'fetch'; empty if no errors.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

newconnection = false;
msg1 = '';
msg2 = '';
data = {};
switch class(datasrc)
    case 'database'
        conn = datasrc;
    case 'char'
        conn = database(datasrc, '', '');
        newconnection = true;
        if isequal(conn.Handle, 0)
            error('dg_runSQL:baddb', ...
                'Could not open database "%s".\n"%s"', datasrc, conn.Message );
        end
    otherwise
        error('dg_runSQL:badarg', ...
            '<datasrc> must be either a connection object or DSN string.' );
end

% extract the SQL command
sqlstr = strtrim(sqlstr);
delimidx = regexp(sqlstr, '\s');
sqlcmd = sqlstr(1:(delimidx-1));

try
    curs = exec(conn, sqlstr);
    msg1 = curs.Message;
    if isempty(msg1) && ~isempty(regexpi(sqlcmd, '^(select$)'))
        curs = fetch(curs);
        msg2 = curs.Message;
        data = curs.Data;
    end
catch
    if newconnection
        close(conn);
    end
    return
end

if newconnection
    close(conn);
end
