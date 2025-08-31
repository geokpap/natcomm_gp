function result = dg_fileExists(filedir, filename, icflag)
%result = dg_fileExists(filedir, filename, icflag)
%result = dg_fileExists(filedir, filename)
%INPUTS
% filedir: pathname of directory to search.
% filename: a string to match against the contents of <filedir>.
% icflag: if not given, defaults to <true>.  If true, <filename> is
%   interpreted as a regular expression and matched case-insensitively, and
%   <result> is a cell string array containing the actual filename(s) that
%   matched <filename> (or {} if there are no matches).  If false, returns
%   <filename> in <result> if file exists exactly as specified, and '' if
%   not.
%OUTPUTS
% result: empty if file does not exist, non-empty if it does.  See <icflag>
%   for details.

%$Rev: 182 $
%$Date: 2013-10-20 22:01:58 -0400 (Sun, 20 Oct 2013) $
%$Author: dgibson $

if nargin < 3
    icflag = true;
end
if icflag
    files = dir(filedir);
    files(cell2mat({files.isdir})) = [];
    fnames = {files.name};
    ismatch = ~cellfun(@isempty, regexpi(fnames, filename));
    result = fnames(ismatch);
else
    if exist(fullfile(filedir, filename), 'file')
        result = filename;
    else
        result = '';
    end
end
