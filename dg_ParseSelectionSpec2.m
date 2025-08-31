function evalstring = dg_ParseSelectionSpec2(specstring)
%dg_ParseSelectionSpec2(specstring)
% <specstring> has same syntax as for dg_ParseSelectionSpec. Returns a
% string that can be passed to eval within dg_scoreTrials to determine
% whether a trial satisfies the selection criteria, i.e. to test if an
% event exists we evaluate
%   any(trialdata(trialnum).events(:,2)==<evtID>)
% For comparison, the string returned by dg_ParseSelectionSpec is
%   trialdata(trialnum).events(<evtID>)~=0
% where <evtID> is some event ID specified in <specstring>. Throws an error
% if parsing fails. 
% <specstring> must be of the form e.g.
%   (31&17)|~18
% or, described in syntax form,
%   specstring ::= '(' specstring ')'
%   specstring ::= specstring '&' specstring
%   specstring ::= specstring '|' specstring
%   specstring ::= '~' specstring
%   specstring ::= integer
% Note that event TTL ID 1 is treated as special, because it
% is the "Disc On" ("Start Recording") event, and therefore always is
% implicitly present and has timestamp 0, whereas for other events a 0
% timestamp means that the event did not occur.
% Adapted from dg_ParseSelectionSpec 04-Apr-2005 by DG.  

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

evalstring = '';
for k = find(specstring == '~')
    if k >= length(specstring) || ( ...
            specstring(k+1) ~= '(' && ...
            ~ismember(specstring(k+1), '1234567890') )
        error('dg_ParseSelectionSpec2:BadTilde', ...
            'Tilde ("~") can only precede an event ID or a parenthesis; we have "%s"', ...
            specstring(k:end));
    end
end
doublethese = find(specstring == '&' | specstring == '|');
for k = (length(doublethese) : -1 : 1)
    specstring = [ specstring(1:doublethese(k)) ...
            specstring(doublethese(k):end) ];
end
while specstring
    % append any leading nondigits to evalstring
    digits = find(ismember(specstring, '1234567890'));
    if length(digits) > 0 && digits(1) > 1
        evalstring = [ evalstring specstring(1 : digits(1) - 1) ];
    elseif length(digits) == 0
        % specstring contains no digits, so it must be trailing punctuation
        evalstring = [ evalstring specstring ];
    elseif digits(1) == 1
        % specstring starts with digits, no action required
    else
        % this should never happen: length(digits) < 0 or digits(1) < 1
        error('dg_ParseSelectionSpec2:BadParsing', ...
            'Program Error parsing selection spec' );
    end
    % get next token and complain if it is not an integer
    [id, specstring] = strtok(specstring, '()&|~');
    nondigits = strtok(id, '1234567890');
    if nondigits
        error('dg_ParseSelectionSpec2:BadEventID', ...
            'Expected an integer at "%s"', id);
    end
    % expand the token and append it to evalstring
    if strcmp(id, '1')
        evalstring = [ evalstring 'true' ];
    elseif id
        evalstring = [ evalstring ...
            'any(trialdata(trialnum).events(:,2)==' id ')'];
    end
end