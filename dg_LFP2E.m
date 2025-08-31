function LFP2E = dg_LFP2E(sessiondir, animalID, patchtable, ...
    fnameregexp, csctable, recurse)
%dg_LFP2E(sessiondir, animalID, patchtable)
%dg_LFP2E(..., fnameregexp, csctable, recurse)
% Renames all the LFP files in <sessiondir> with names of the form
% E<n>.<ext>, where <n> is the number of the Single Electrode (SE) channel
% that was the source of the LFP recording. The SE numbers are looked up in
% the union of two tables: the tab-delimited  "CSC channel to ADChannel to
% Sc channel" file saved by 'makeconfig.pl' in <sessiondir> (or in its root
% dir if <sessiondir> is a fragment) and the Excel workbook read from
% <patchtable>.
%INPUTS
% sessiondir: The session root dir of <sessiondir> must be like
%   '2005-1-13_17-53-28' (i.e. a match to the regexp
%   '^\d{4}-\d+-\d+_\d+-\d+-\d+$') or 'G011705' (i.e. a match to the regexp
%   '\d{6}$').
% animalID: any number of letters, may not contain non-letters.
% patchtable: absolute or relative filename for Excel workbook containing
%   the patch panel information for CSC numbers 1 - 24 (this is for Joey
%   data).  The workbook must contain spreadsheets whose names are the
%   concatenation of an alphabetic animal ID (any number of letters) and an
%   optional string of non-letters.  <animalID> must exactly (but not
%   case-sensitively) match the animal ID portion of the relevant
%   spreadsheet(s).  An error is raised if there is no entry corresponding
%   to <sessiondir>, <animalID> in the workbook. Alternatively,
%   <patchtable> may point to a tab-delimited text file containing only the
%   relevant sheet, in which case the filename must end with ".tab" or
%   ".txt" (this is for Theresa data).  In this case, there must not be any
%   extra columns beyond the end of the patch table, and the patch table
%   can contain as many columns as needed.  The CSC numbers must be listed
%   on row 4, and dates and electrode numbers must begin on row 5.
%   Regardless of the file format (text or Excel), dates must be in col. 1
%   and electrode numbers must start in col. 2.
% fnameregexp: (optional) <fnameregexp> is a regular expression that
%   case-insensitively matches the filenames of LFP files; default value
%   'LFP\d+\.(MAT|NCS)'.  It is assumed that any filenames that match
%   <fnameregexp> will also match '\D\d+\.', i.e. will contain a string of
%   nothing but digits terminated by a period, and that the matching string
%   of digits represents the LFP number.
% csctable: (optional) Name (relative to the root session directory) of
%   the "CSC channel to ADChannel to Sc channel" table saved by
%   'makeconfig.pl'; default value 'csctable.xls'.  Must NOT be an absolute
%   pathname.
% recurse: (optional) A flag that controls whether or not subdirectories
%   are also processed; default true.
%OUTPUTS
% LFP2E is a numeric vector that maps from lfpnum to SE number (i.e. each
% element contains the corresponding SE number).

%$Rev: 239 $
%$Date: 2016-03-10 19:30:16 -0500 (Thu, 10 Mar 2016) $
%$Author: dgibson $

if nargin < 6
    recurse = true;
end
if nargin < 5
    csctable = 'csctable.xls';
end
if nargin < 4
    fnameregexp = 'LFP\d+\.(MAT|NCS)';
end
animalID = upper(animalID);
[~, ~, ext] = fileparts(patchtable);
if ismember(lower(ext), {'.txt', '.tab'})
    % Read tab-delimited text file
    txt = dg_tabreadtxt(patchtable);
    % extract <dates> and <trodes>
    col1data = txt(5:end, 1);
    gotdate = ~cellfun(@isempty, regexp(col1data, '^[0-9/]*$'));
    dates = NaN(size(gotdate));
    dates(gotdate) = datenum(col1data(gotdate), 23);
    trodes = str2double(txt(5:end, 2:end));
else
    % Read the Excel file
    [typ, desc] = xlsfinfo(patchtable);
    if ~isequal(typ, 'Microsoft Excel Spreadsheet')
        error('dg_LFP2E:badxls', ...
            '%s is not an Excel workbook', patchtable );
    end
    desc = upper(desc);
    pattern = sprintf('^%s[^A-Za-z]', animalID);
    rightAnimal = find(~cell2mat(dg_mapfunc(@isempty, regexp(desc, pattern))));
    if isempty(rightAnimal)
        error('dg_LFP2E:animalID', ...
            'No such animal ID');
    end
    % extract <dates> and <trodes>
    dates = [];
    trodes = [];
    for k = reshape(rightAnimal, 1, [])
        [num, txt] = xlsread(patchtable, desc{k});
        if size(num,1) ~= size(txt,1)
            % xlsread fails to pad num with trailing rows of NaN when num
            % contains trailing rows of empty cells; fix it
            addon = NaN((size(txt,1)-size(num,1)), size(num,2));
        else
            addon = [];
        end
        dates = [dates; NaN; datenum(txt(2:end,1), 23)]; %#ok<*AGROW>
        trodes = [trodes; num; addon];
    end
end

% Find the right row of patchtable
sessionroot = dg_findSessionRoot(sessiondir);
[~, sessionID] = fileparts(sessionroot);
if regexp(sessionID, '^\d{4}-\d+-\d+_\d+-\d+-\d+$')
    % Cheetah style
    underscore = regexp(sessionID, '_');
    sessiondatenum = datenum(sessionID(1:underscore-1), 29);
elseif regexp(sessionID, '\d{6}$')
    % Theresa style
    sessiondatenum = datenum( ...
        [sessionID(end-5:end-4) '/' sessionID(end-3:end-2) '/' sessionID(end-1:end)], ...
        23);
else
    error('dg_LFP2E:badfmt', ...
        'Unrecognizable format of session ID "%s"', sessionID);
end
rownum = find(dates == sessiondatenum);
if length(rownum) > 1
    if size(txt,2) < 26
        error('dg_LFP2E:patchtable', ...
            'Patchtable contains ambiguous entries');
    end
    rownum = find(ismember(txt(:,26), sessionID));
    if length(rownum) > 1
        error('dg_LFP2E:patchtable2', ...
            'Patchtable contains multiple listings of %s', sessionID);
    end
end
LFP2E = trodes(rownum, :);

% Add the entries from csctable:
if exist(fullfile(sessionroot, csctable), 'file')
    csctabnum = dg_tabread(fullfile(sessionroot, csctable));
    LFP2E(csctabnum(:,1)) = csctabnum(:,3);
    if any(csctabnum(:,1)<25)
        error('dg_LFP2E:badCSCtable', ...
            'CSC table %s contains CSC numbers that belong to patch panel', ...
            fullfile(sessionroot, csctable) );
    end
else
    warning('dg_LFP2E:noCSCtable', ...
        'The file %s does not exist', fullfile(sessionroot, csctable));
end
recursive_LFP2E(LFP2E, sessiondir, recurse, fnameregexp);
end

function recursive_LFP2E(LFP2E, sessiondir, recurse, fnameregexp)
if ispc
    cmdname = 'move';
else
    cmdname = 'mv';
end
files = dir(sessiondir);
for fnum = 1:length(files)
    fname = files(fnum).name;
    if files(fnum).isdir && ~ismember(fname, {'.' '..'})
        if recurse
            recursive_LFP2E(LFP2E, fullfile(sessiondir, fname), ...
                recurse, fnameregexp);
        end
    elseif regexpi(fname, fnameregexp)
        digits = regexp(fname, '\D\d+\.') + 1;
        dot = regexp(fname(digits:end), '\.');
        lfpnum = str2double(fname(digits:(digits+dot-2)));
        newname = sprintf('E%.0f%s', LFP2E(lfpnum), ...
            fname(digits+dot-1:end));
        % Matlab 'movefile' and the following formulation both silently
        % blow away any existing file at <newname>:
        cmdstr = sprintf('!%s "%s" "%s"', cmdname, ...
            fullfile(sessiondir, fname), fullfile(sessiondir, newname));
        eval(cmdstr);
    end
end
end

