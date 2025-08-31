function unitlist = dg_xls2unitlist(filename, cols, rootdir, ...
    animalprefix, sessionprefix, varargin)
% unitlist = dg_xls2unitlist(filename, cols, rootdir, ...
%   animalprefix, sessionprefix)
% Reading from columns <cols> in the Excel file <filename>, list all units
% in the file as a lookup table that associates a session directory's
% absolute pathname with a list of unit IDs in that session (i.e. list the
% units in the format used for <theunits> by mycohero5 and other funcs in
% bulkProcess_bill).  <cols> must contain exactly four column numbers,
% which represent animal number, session number, tetrode number, and
% cluster number (in that order).  <rootdir> is a directory that contains
% all the animal directories.  <animalprefix> is a string to which the
% animal number is appended to yield an animal ID, e.g. if <animalprefix>
% is 'S', then animal #17 would be 'S17'.  <sessionprefix> works similarly
% for session IDs, so <sessionprefix> = 'ACQ' plus session #3 yields
% 'ACQ03'.
%OPTIONS
% 'tab' - uses dg_tabreadfix instead of xlsread, so in this case 
%   <filename> must point to a tab-delimited text file, not a native
%   Excel file.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if numel(cols) ~= 4 || ~isequal(class(cols), 'double')
    error('<cols> must contain 4 numbers.');
end

if isempty(animalprefix)
    error('<animalprefix> must not be empty.');
end

if isempty(sessionprefix)
    error('<sessionprefix> must not be empty.');
end

argnum = 1;
tabflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'tab'
            tabflag = true;
        otherwise
            error('lfp_disp:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end


animalprefix = upper(animalprefix);
sessionprefix = upper(sessionprefix);
if tabflag
    num = dg_tabreadfix(filename);
    num(1,:) = [];
else
    num = xlsread(filename);
end
animalnums = unique(num(:, cols(1)))';
if isempty(animalnums)
    error('dg_xls2unitlist:e1', ...
        'There are no animal numbers in column %d.', cols(1) );
end

unitlist = {};
for animaln = animalnums
    animalID = sprintf('%s%02d', animalprefix, animaln);
    animalrows = find(num(:, cols(1)) == animaln);
    sessionnums = unique(num(animalrows, cols(2)))';
    for sessionn = sessionnums
        sessionname = sprintf('%s%02d', sessionprefix, sessionn);
        sessionunitIDs = {};
        sessionrows = animalrows(num(animalrows, cols(2)) == sessionn);
        trodenums = unique(num(animalrows, cols(3)))';
        for troden = trodenums
            troderows = sessionrows(num(sessionrows, cols(3)) == troden);
            clustnums = num(troderows, cols(4))';
            for clustn = clustnums
                unitID = sprintf('%s-T%dC%d', sessionname, troden, clustn);
                sessionunitIDs(1, end+1) = {unitID};
            end
        end
        unitlist(end+1, :) = ...
            { fullfile(rootdir, animalID, sessionname), sessionunitIDs };
    end
end


