function wrappedstr = dg_wraptext(str, maxlen, mode)
% Searches for whitespace characters in str that can serve as line breaks
% in order to fit the length of each line within maxlen.  Breaks <str> into
% a series of lines that fit.  If
% no whitespace occurs for <maxlen> chars in a row, then the break is made
% at <maxlen>.  The next line starts with the next non-whitespace char.
% Handling of leading and trailing whitespace is undefined.  <mode> is
% optional and defaults to 1.  If <mode> is 1, then the liens are returned
% as a cell string vector; if <mode> is 2, then the lines are returned as a
% single long string with embedded '\n' chars.

%$Rev: 140 $
%$Date: 2012-01-06 23:01:07 -0500 (Fri, 06 Jan 2012) $
%$Author: dgibson $

if nargin < 3
    mode = 1;
end

if isempty(str)
    switch mode
        case 1
            wrappedstr = {''};
        case 2
            wrappedstr = '';
    end
    return
end

wsidx = regexp(str, '\s+');
startidx = 1;
linenum = 1;
while startidx <= length(str)
    if length(str) - startidx + 1 <= maxlen
        endidx = length(str);
    else
        % <brkidx2> is idx into <wsidx>, and is thus an idx2 into str:
        brkidx2 = find(wsidx > startidx & wsidx <= startidx + maxlen);
        if isempty(brkidx2)
            endidx = min(length(str), startidx + maxlen - 1);
        else
            endidx = wsidx(brkidx2(end)) - 1;
        end
    end
    wrappedcellstr{linenum} = str(startidx:endidx);
    linenum = linenum + 1;
    % find next non-whitespace starting with first char after break
    startidx = endidx + 1;
    while startidx <= length(str) && ~isempty(regexp(str(startidx), '\s'))
        startidx = startidx + 1;
    end
    if startidx > length(str)
        break
    end
end
switch mode
    case 1
        wrappedstr = wrappedcellstr;
    case 2
        wrappedstr = wrappedcellstr{1};
        if length(wrappedcellstr) > 1
            for k = 2:length(wrappedcellstr)
                wrappedstr = ...
                    sprintf('%s\n%s', wrappedstr, wrappedcellstr{k});
            end
        end
end
