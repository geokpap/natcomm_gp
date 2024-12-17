function [pairs, extraL, extraR] = dg_zip(L, R, allocsize)
%DG_ZIP interleaves the elements of two sorted arrays
%[pairs, extraL, extraR] = dg_zip(L, R, allocsize)
% Finds <pairs> of indices into L and R such that L(pair(n,1)) <=
% R(pair(n,2)) < L(pair(n+1,1)).  I.e. index each pair points to two values
% where the first is less than or equal to the second.  It is assumed that
% L and R are in strictly increasing order when addressed as vectors, and
% each contains at least one element. Behavior for repeated values (i.e.
% non-decreasing order) is undefined. If there is more than one value in
% list R that is bracketed by two successive values in L, then the first
% value is used and the extra values are added to the <extraR> list.
% However, if there are extra values in L, then the LAST value is used and
% the others are added to <extraL>.  The end result of this policy is that
% each <pair> is as close in value as possible. Like <pairs>, <extraL> and
% <extraR> are lists of indices into L and R.  The orientations (row or
% column) of <extraL> and <extraR> are undefined. <allocsize> is optional
% with default = 256, and it specifies the number of rows to allocate at a
% time to <pairs>. Unused rows are trimmed off before returning.

%$Rev: 51 $
%$Date: 2010-04-27 17:36:11 -0400 (Tue, 27 Apr 2010) $
%$Author: dgibson $

if nargin < 3
    allocsize = 256;
end

pairs = zeros(allocsize, 2);
extraL = [];
extraR = [];
Lidx = 1; % index into L
Ridx = 1; % index into R
pairsi = 1; % index into pairs
done = false;   % true when no more pairs can be found

% At the start of each iteration, Lidx and Ridx point to a putative pair:
while ~done
    
    % Find match for Lidx:
    while ~done && (R(Ridx) < L(Lidx))
        if Ridx == numel(R)
            % done with R; there is no match for Lidx
            done = true;
        end
        extraR(end+1) = Ridx;
        Ridx = Ridx + 1;
    end
    if done 
        break
    end
    
    % Now R(Ridx) > L(Lidx), so Ridx is the match for Lidx.  Allocate more pairs if
    % necessary and save the pair:
    if pairsi > size(pairs,1)
        pairs = [ pairs; zeros(allocsize, 2) ];
    end
    pairs(pairsi, :) = [Lidx Ridx];
    pairsi = pairsi + 1;
    
    % Attempt next Lidx:
    Lidx = Lidx + 1;
    if Lidx > numel(L)
        done = true;
    end
    
    % Find match for Ridx; if L(Lidx) is not strictly greater, then update
    % the previously recorded match to this Ridx (this allows equal L(Lidx)
    % and R(Ridx) to be a pair):
    while ~done && (L(Lidx) <= R(Ridx))
        if Lidx == numel(L)
            % done with L; there is no match for Ridx
            done = true;
        end
        extraL(end+1) = pairs(pairsi-1, 1);
        pairs(pairsi-1, 1) = Lidx;
        Lidx = Lidx + 1;
    end
    Ridx = Ridx + 1;
    if Ridx > numel(R)
        done = true;
    end
end

% All pairs have been found.  If there are no extras, then Ridx should point
% at R(end+1) and Lidx should point at L(end+1).  Handle the extras if they
% exist.
if Ridx < numel(R) + 1
    extraR(end+1:end+numel(R)-Ridx+1) = Ridx:numel(R);
    if Lidx < numel(L) + 1
        % Find match for Ridx:
        while (L(Lidx) <= R(Ridx))
            extraL(end+1) = pairs(pairsi-1, 1);
            pairs(pairsi-1, 1) = Lidx;
            Lidx = Lidx + 1;
        end
    end
end
if Lidx < numel(L) + 1
    extraL(end+1:end+numel(L)-Lidx+1) = Lidx:numel(L);
end

% Trim unused rows
unused = find(pairs(:,1)==0);
if ~isempty(unused)
    pairs(unused(1):end,:) = [];
end
