function [origtrialnums, sessionname] = lfp_getSessionInfo(varargin)
%[origtrialnums, sessionname] = lfp_getSessionInfo(varargin)
%  Reads origtrialnums, sessionname from lfp_FragmentFiles.mat if it exists
%  (indicating that lfp_DataDir points at a fragment directory); otherwise
%  creates sessionname and returns [] for origtrialnums. Removes any blanks
%  from sessionname and converts all chars to lowercase before returning.
%  <varargin> is normally empty, but there are some difficult cases that
%  require more than just lfp_DataDir to determine the sessionname:
%   rodent cluster files: <varargin> = {'filename'}, where <filename> is
%       the name of the .TTn or .Tnn file; the file's basename (i.e.
%       without the directory or extension) is the session ID, and the
%       current directory name is the animal ID.  BUT: this does NOT apply
%       to rodent *.MAT files, which we assume are Multiple Cut Cluster
%       format, and reside in the rodent session directory, and should
%       therefore be handled like rodent raw data files.
%   rodent tracker files: same as for rodent clusters, except that it is
%       assumed that the initial "V" in the filename is removed first
%       before calling lfp_getSessionInfo.

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

lfp_declareGlobals;

names = {'lfp_FragmentFiles.mat' 'lfp_fragmentfiles.mat'};
for nameidx = 1:length(names)
    infofilepath = fullfile(lfp_DataDir, names{nameidx});
    if exist(infofilepath, 'file')
        load(infofilepath, '-mat');
        origtrialnums = trialnums;
        sessionname = sessionName;
        break
    end
end
if ~exist(infofilepath, 'file') || ...
        ~isempty(lfp_GetSessionNameFromDir) && lfp_GetSessionNameFromDir
    if ~exist(infofilepath, 'file')
        origtrialnums = [];
    end
    [pathstr,dirname] = fileparts(lfp_DataDir);
    if (length(varargin) == 1)
        [p,d,ext] = fileparts(varargin{1});
    end
    switch lfp_DirStructure
        case 'rodent'
            if length(varargin) == 0 || ( length(varargin) == 1 ...
                    && isequal(upper(ext), '.MAT') )
                % Must get animal ID from parent directory.
                % lfp_DataDir should be pointing at session directory.
                [parentpathstr,parentname] = fileparts(pathstr);
                sessionname = [ parentname dirname ];
            elseif length(varargin) == 1
                % Must get session ID from file name.
                % lfp_DataDir should be pointing at animal directory.
                [pathstr3,basename,ext] = fileparts(varargin{1});
                if regexpi(ext, '\.T+[0-9]+$')
                    % Rodent cluster file
                    sessionname = [ dirname basename ];
                elseif regexpi(ext, '\.DAT$')
                    % Rodent tracker file, may start with a digit
                    % indicating some kind of data repair
                    if ~isempty(regexp(basename, '^\d', 'once' ))
                        basename = basename(2:end);
                    end
                    sessionname = [ dirname basename ];
                elseif regexpi(ext, '\.(MAT|DWD|T|NSE)$')
                    % Multiple Cut Cluster or Single Cut Cluster Matlab
                    % file, Naotaka Clusters, Single Cut Cluster MClust
                    % file, Single Electrode Neuralynx
                    sessionID = basename(1 : strfind(basename, '-') - 1);
                    sessionname = [ dirname sessionID ];
                else
                    error('lfp_getSessionInfo:badFileExt', ...
                        'File extension not recognized: "%s"', ext );
                end
            end
        case 'monkey'
            % This is for use in Joey's late 2010 environment, where
            % session root directory names are like 'hh092807':
            sessionroot = dg_findSessionRoot(lfp_DataDir);
            [p, sessionname] = fileparts(sessionroot);
        otherwise
            % directory name includes both session ID and animal ID, and
            % all files are in one directory
            sessionname = dirname;
    end 
end
% Kill blanks, which are delimiters in Unique Trial ID lists:
sessionname = dg_substitute(lower(sessionname), ' ', '');
