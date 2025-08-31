function [starts, ends] = lfp_theresa_debounce(rawstarts, rawends, events)
%[starts, ends] = lfp_theresa_debounce(rawstarts, rawends, events)
% Debounce target entries and exits, so that repeated entries and exits
% to/from the same target are reduced to just the first entry and last
% exit.  <rawstarts> and <rawends> are lists of indices into <events>,
% which is in the same format as lfp_Events.  Debounced lists of indices
% are constructed  in <starts> and <ends>.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

starts = [];
ends = [];
k = 1;
k2 = 1;
% there is a fixation end (0x20) that is the target entered
% before the grid turned on and want to ignore that
if events(rawends(1),1) < events(rawstarts(1),1)
    rawends = rawends(2:end);
end
if length(rawends) < (length(rawstarts) - 1)
    warning('lfp_theresa_debounce:mismatch', ...
        '>1 extra start, ignoring; %d starts, %d ends', ...
        length(rawstarts), length(rawends) );
    rawstarts(length(rawends)+1 : end) = [];
end
for i = 1 : length(rawstarts)
    if i == 1 %keep the first fix start
        starts(k) = rawstarts(i);
        k = k + 1;
    end
    if i ~= length(rawstarts)
        if events(rawstarts(i),2) == events(rawstarts(i+1),2) %same target entered next
            timedifference = events(rawstarts(i+1),1) - events(rawstarts(i),1);
            if timedifference > 0.150 %150ms - may have to adjust this
                starts(k) = rawstarts(i+1);
                ends(k2) = rawends(i);
                k = k + 1;
                k2 = k2 +1;
            end
        else %different target entered next
            starts(k) = rawstarts(i+1);
            ends(k2) = rawends(i);
            k = k + 1;
            k2 = k2 +1;
        end
    else % i = length(rawstarts) -- last one
        if length(rawstarts) == length(rawends) %every entrance should have exit
            ends(k2) = rawends(i);
%         else
%             warning('lfp_theresa_debounce:unmatched', ...
%                 '%d rawstarts, %d rawends', ...
%                 length(rawstarts), length(rawends) );
        end
    end
end

