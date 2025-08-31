function [table, trialIDs] = dg_readEyeTab(filename)
%[trialIDs, table] = dg_readEyeTab(filename)
% Reads files produced by lfp_EyeTabulation and returns values in same
%   format returned by lfp_EyeTabulation.
% <filename> is the relative or absolute pathname to a tab-delimited text
%   file.  If the file has seven columns, it is a fixations file.   If it
%   has 8 columns, then it is a saccades file.  If it has some other number
%   of columns, a warning is issued and processing continues.
% <table> is a cell array equal to the <fixations> or <saccades> value
%   returned by lfp_EyeTabulation.
% <trialIDs> is a cell string array of the Unique Trial IDs whose data
%   produced the fixations or saccades file.  Note that trials for which no
%   fixations or saccades were recorded do not appear in the files, so it
%   is possible for there to be a trial that was selected that nonetheless
%   does not appear in <trialIDs>.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

% Get the number of columns
fid = fopen(filename);
if fid == -1
    error('dg_readEyeTab:badfile', ...
        'Could not open "%s" for reading', filename );
end
line = fgetl(fid);
fclose(fid);
if line == -1
    error('dg_readEyeTab:badfile2', ...
        '"%s" is empty.', filename );
end
tabidx = regexp(line, '\t');
ncol = length(tabidx) + 1;
if ~ismember(ncol, [7 8])
    warning('dg_readEyeTab:unknown', ...
        '"%s" contains %d columns.', filename, ncol );
end

fid = fopen(filename);
if fid == -1
    error('dg_readEyeTab:badfile3', ...
        'Could not re-open "%s" for reading', filename );
end
C = textscan(fid, '%s', 'delimiter', '\t\n');
fclose(fid);
% These files may (and usually do) contain a trailing empty line, which
% should be axed:
if isempty(C{1}{end}) && (mod(numel(C{1}) - 1, ncol) == 0)
    C{1}(end) = [];
end
if mod(numel(C{1}), ncol)
    error('dg_readEyeTab:badfile4', ...
        'The text file must contain the same number of fields on every row');
end
C = reshape(C{1}, ncol, [])';

% C is now a cell array of the same shape as the table that was read, with
% trialnums in first column, Unique Trial IDs in the second, and fixations
% or saccades values in the remaining columns.

% Get trialIDs in the same order in which they occur in the file, which is
% not necessarily sorted; but the trialnums in col 1 are sorted, and can
% speed it up:
trialnums = zeros(size(C,1),1);
for k = 2:length(trialnums)
    trialnums(k) = sscanf(C{k,1}, '%d');
end
firstrows = [ 2
    find(trialnums(3:end) ~= trialnums(2:end-1)) + 2 ];
trialIDs = C(firstrows,2);

% Construct output table:
table = cell(1, length(trialIDs));
for trialidx = 1:length(trialIDs)
    if trialidx == length(trialIDs)
        lastrow = size(C,1);
    else
        lastrow = firstrows(trialidx+1) - 1;
    end
    rows = firstrows(trialidx):lastrow;
    data = [];
    for ridx = 1:length(rows)
        for col = 3:ncol
            data(ridx,col-2) = sscanf(C{rows(ridx),col}, '%g');
        end
    end
    table{trialidx} = data;
end

