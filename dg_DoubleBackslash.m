function newstring = dg_DoubleBackslash(string)
%newstring = dg_DoubleBackslash(string)
% Returns a copy of string with each backslash replaced by two backslashes.
% Handy for escaping DOS pathnames in TeX strings.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

newstring = [];
bs = find(string == '\');
previous = 0;
for current = bs
    newstring = [ newstring string(previous+1:current-1) '\\' ];
    previous = current;
end
newstring = [ newstring string(previous+1:end) ];