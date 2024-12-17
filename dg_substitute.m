function result = dg_substitute(str1, str2, str3)
%DG_SUBSTITUTE performs a global search-and-replace operation on a string.
%dg_substitute(str1, str2, str3)
%  Search str1 for any instances of str2, and replace each instance with
%  str3.  The search is case-sensitive.  If the str1 and str2 contain
%  repeats such that str2 matches some of the same character positions
%  multiple times, only the left-most of the overlapping matches is used.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

result = str1;
matches = strfind(result, str2);
matchends = matches + length(str2) - 1;
if length(matches) > 1
    overlapped = find(matches(2:end) <= matchends(1:end-1)) + 1;
    while ~isempty(overlapped)
        matches(overlapped(1)) = [];
        if length(matches) > 1
            matchends = matches + length(str2) - 1;
            overlapped = find(matches(2:end) <= matchends(1:end-1)) + 1;
        else
            overlapped = [];
        end
    end
end
        
for m = length(matches) : -1 : 1
    result = [ result(1 : matches(m)-1) str3 result(matchends(m)+1 : end)];
end