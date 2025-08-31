function dg_txt2nev(txtfilename, nevfilename, varargin)
%dg_txt2nev(txtfilename, nevfilename, 'createstrings')
%dg_txt2nev(txtfilename, nevfilename)
%   Reads a text file of the type created by dg_readEvents(...,
%   'text') and writes it in Neuralynx Events file format.  All consecutive
%   lines that begin with a space character are considered to be part of
%   the file header, and the leading space is removed from each line before
%   writing the file.  The first line that does NOT begin with a space is
%   considered to be the start of the data records.
% OPTIONS
%   'createstrings' - ignores the event strings column in <txtfilename> and
%       creates new ones based on the numeric values in the TTL IDs column.

%$Rev: 102 $
%$Date: 2011-04-20 17:17:25 -0400 (Wed, 20 Apr 2011) $
%$Author: dgibson $

createstringsflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'createstrings'
            createstringsflag = true;
        otherwise
            error('dg_txt2nev:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

infid = fopen(txtfilename);
if infid == -1
    error('dg_txt2nev:badfile', ...
        'Could not open %s for reading', txtfilename );
end

allocsize = 1024;
ts = NaN(1,allocsize);
ttl = NaN(1,allocsize);
estr = cell(allocsize,1);

% Read Header
hdr = {};
tline = fgetl(infid);
linenum = 1;
while ~isempty(regexp(tline, '^\s.*$'))
    hdr{end+1, 1} = tline(2:end);
    tline = fgetl(infid);
    linenum = linenum + 1;
end
while isempty(tline) || ~isempty(regexp(tline, '^\s*$'))
    tline = fgetl(infid);
    linenum = linenum + 1;
end

% Read Events Data
nextrecnum = 1;
while ~isequal(tline, -1)
    tabidx = regexp(tline, '\t');
    if length(tabidx) ~= 2
        warning('dg_txt2nev:badinput', ...
            'Line %d does not contain two tabs; skipping.', ...
            linenum );
    else
        if nextrecnum > length(ts)
            ts = [ ts NaN(1,allocsize) ];
            ttl = [ ttl NaN(1,allocsize) ];
            estr = [ estr; cell(allocsize,1) ];
        end
        ts(nextrecnum) = str2num(tline(1:(tabidx(1)-1)));
        ttl(nextrecnum) = str2num(tline((tabidx(1)+1):(tabidx(2)-1)));
        estr(nextrecnum) = {tline((tabidx(2)+1):end)};
        nextrecnum = nextrecnum + 1;
    end
    tline = fgetl(infid);
    linenum = linenum + 1;
end
fclose(infid);
unusedlines = isnan(ts);
ts(unusedlines) = [];
ttl(unusedlines) = [];
estr(unusedlines) = [];
        
% Write Events Data
if createstringsflag
    estr = [];
end
dg_writeNlxEvents(nevfilename, ts, ttl, estr, hdr);


