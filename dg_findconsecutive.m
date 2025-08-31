function [start, finish] = dg_findconsecutive(v, index)
% Returns the index of the first and last elements of v starting at index
% that are in consecutive numerical order.  If there is no such subseries in
% v, finish is returned empty.  <start> and <finish> are scalars.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

start = [];
finish = [];
while isempty(start) && (index < length(v))
    if v(index+1) == v(index) + 1
        start = index;
    end
    index = index + 1;
end
accumulating = false;
while isempty(finish) && index <= length(v)
    if v(index) == v(index-1) + 1
        accumulating = true;
    elseif accumulating && (index < length(v))
        % Implicitly, v(index) ~= v(index-1) + 1
        finish = index - 1;
    elseif accumulating
        % Implicitly, v(index) ~= v(index-1) + 1 && index == length(v)
        finish = length(v) - 1;
    end
    index = index + 1;  
end
if isempty(finish) && accumulating
    % Special case: we ran off the end while still accumulating
    finish = length(v);
end