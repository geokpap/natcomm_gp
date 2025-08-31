function dg_findDupRoots_processline(thisline, prevline, linenum)
% Execute dg_findDupRoots Step 2 for the file listed in <thisline>.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

global R

if isequal(thisline, prevline)
    warning('dg_findDupRoots_processline:dupdata', ...
        'Skipping duplicate data entry at infile line %d', linenum);
    return
end

toks = regexp(thisline, ...
    '^([a-f0-9]+) .* +[0-9]+ +[a-zA-Z0-9_]+ +[a-zA-Z0-9_]+ +[0-9]+ +[A-Za-z]* *[-:0-9]+ +[-:0-9]+ (.*)$', ...
    'tokens');
md5hash = toks{1}{1};
pathname = toks{1}{2};
[dirpath, basename, ext] = fileparts(pathname);
filename = [basename, ext];

Ridx = dg_findDupRoots_processline_finddir(dirpath);
if isempty(Ridx)
    Ridx = dg_findDupRoots_processline_adddir(dirpath);
end
if isempty(prevline)
    % We are at the beginning of the file.  There is nothing to search, and
    % no duplicate to point to.
else
    % Check md5hash and file name from <prevline> to see if it is a dup of
    % <pathname>.
    prevtoks = regexp(prevline, ...
        '^([a-f0-9]+) .* +[0-9]+ +[a-zA-Z0-9_]+ +[a-zA-Z0-9_]+ +[0-9]+ +[A-Za-z]* *[-:0-9]+ +[-:0-9]+ (.*)$', ...
        'tokens');
    prevhash = prevtoks{1}{1};
    prevpath = prevtoks{1}{2};
    [prevdirpath, basename, ext] = fileparts(prevpath);
    prevfilename = [basename, ext];
    if isequal(prevpath, pathname)
        error('dg_findDupRoots_processline:oops', ...
            'This can''t happen.');
    end
    if isequal(md5hash, prevhash) && isequal(filename, prevfilename)
        % <thisline> denotes a dup of the file in <prevline>.  However,
        % these can be false alarms, so we must verify by actually
        % comparing the files.  There are a few pesky pathnames that end in
        % single-quote, so we axe that if necessary here:
        if pathname(end) == ''''
            pathname(end) = [];
        end
        if prevpath(end) == ''''
            prevpath(end) = [];
        end
        % And then there are pathnames that contain single-quotes, which we
        % make unescapable, or already-escaped, or something:
        pathname = strrep(pathname, '''', '''\''''');
        prevpath = strrep(prevpath, '''', '''\''''');
        cmdstr = sprintf('cmp -s ''%s'' ''%s''', pathname, prevpath);
        %disp(cmdstr);
        rc = system(cmdstr);
        switch rc
            case 0
                % The files are the same; enter the index of the copy from
                % the previous line into R(Ridx).dupidx if it is not
                % already there.  R(Ridx).dupidx is made to be a column
                % vector.
                previdx = dg_findDupRoots_processline_finddir(prevdirpath);
                if isempty(R(Ridx).dupidx)
                    R(Ridx).dupidx = previdx;
                else
                    if ~ismember(previdx, R(Ridx).dupidx)
                        R(Ridx).dupidx(end+1, 1) = previdx;
                    end
                end
            case 1
                % the files are different; skip <thisline>
                return
            otherwise
                if exist(pathname) && exist(prevpath)
                    % something bad must have happened
                    error('dg_findDupRoots_processline:barf', ...
                        '''cmp'' returned exit code %d', rc);
                else
                    % just a missing file; skip <thisline>
                    return
                end
        end
    else
        % <thisline> denotes the first instance of a new file; do nothing.
    end
end
R(Ridx).files{end+1, 1} = filename;

end


function Ridx = dg_findDupRoots_processline_adddir(dirpath)
% Add the file <pathname> to the structure <R>, adding directories as
% needed in the manner of 'mkdir -p'.  This may thus create more than one
% new record in <R>.
global R
Ridx = 1;
splitpath = strsplit(dirpath, filesep);
% Get rid of the useless empty string at the beginning:
splitpath(1) = [];
% Iterate down <splitpath>.  Keep <Ridx> pointing at the record that
% represents the current subdirectory.
for pathidx = 1:length(splitpath)
    partialpath = fullfile(filesep, splitpath{1:pathidx});
    if ismember(partialpath, R(Ridx).subdirs)
        % This is a subdir we have encountered before.  Do nothing to
        % <R>, just update Ridx and keep iterating down the
        % <splitpath>.
        Ridx = R(Ridx).subdiridx(ismember( ...
            R(Ridx).subdirs, partialpath ));
    else
        % New subdir, create a new record in <R>:
        prevRidx = Ridx;
        Ridx = length(R) + 1;
        R(Ridx).path = partialpath; % this makes <R> one element longer
        R(Ridx).depth = R(prevRidx).depth + 1;
        R(Ridx).subdirs = {};
        R(Ridx).files = {};
        R(prevRidx).subdirs{end+1, 1} = partialpath;
        R(prevRidx).subdiridx(end+1) = Ridx;
    end
end

end


function Ridx = dg_findDupRoots_processline_finddir(dirpath)
% Find the directory in <R> whose .path is <dirpath>.
% If no such directory is found, [] is returned.
global R
splitpath = strsplit(dirpath, filesep);
% Get rid of the useless empty string at the beginning:
splitpath(1) = [];
% Now we iterate down <splitpath> and navigate the links in <R> in
% parallel, until we either find the file or cannot continue (indicating
% there is no such file).
Ridx = 1;
for pathidx = 1:length(splitpath)
    % keep going down the directory tree
    Ridx = R(Ridx).subdiridx(ismember( R(Ridx).subdirs, ...
        fullfile(filesep, splitpath{1:pathidx}) ));
    if isempty(Ridx)
        return
    end
end

end

