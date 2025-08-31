function trialnums = lfp_parseTrialStr(trialstr, session)
%LFP_PARSETRIALSTR converts a string of trial IDs to an array of trialnums.
%trialnums = lfp_parseTrialStr(trialstr, session)
%   <trialstr> is a string containing a series of Unique Trial IDs in the
%   format returned by lfp_getTrialID, separated by spaces or colons.  Each
%   Trial ID is converted to an internal trial number, and the array is
%   expanded with Matlab syntax, i.e. the colons invoke the Matlab colon
%   operator, and spaces simply separate elements of the final array.
%   <session>, if given, is a default session name to use in case <trialid>
%   contains only the original trial number (lfp_OrigTrialNums).  Matlab
%   2-colon syntax is not supported; the result of specifying 'n1:n2:n3' is
%   undefined.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin < 2
    session = '';
end

% remove leading and trailing blanks:
trialstr = deblank(trialstr(end:-1:1));
trialstr = deblank(trialstr(end:-1:1));

delimiters = ': ';
rem = trialstr;
trialnums = [];
prevtrial = [];
foundcolon = false;
while ~isempty(rem)
    [token, rem] = strtok(rem, delimiters);
    trial = lfp_getTrialNum(token, session);
    if foundcolon
        trialnums = [trialnums prevtrial:trial];
        foundcolon = false;
        prevtrial = [];
    else
        trialnums = [trialnums prevtrial];
        prevtrial = trial;
    end
    if ~isempty(rem)
        % remove and process leading delimiters
        while ismember(rem(1), delimiters)
            if rem(1) == ':'
                foundcolon = true;
            end
            rem(1) = [];
        end
    end
    if isempty(rem)
        if foundcolon
            error('lfp_parseTrialStr:danglingcolon', ...
                'Trial string ended with colon' );
        end
    end
end
trialnums = [trialnums prevtrial];

