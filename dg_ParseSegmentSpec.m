function specstruct = dg_ParseSegmentSpec(specstring)

%DG_PARSESEGMENTSPEC Accepts a string and returns a structure representing
%a time segment specification.  Throws an error if parsing fails.
% specstruct = dg_ParseSegmentSpec(specstring) returns this structure:
%   spectstruct.start - a TimeSpec (see dg_ParseTimeSpec)
%   specstruct.stop - a TimeSpec
% The input specstring must be of the form e.g.
%   '31,38+200to17,18-15'
% or, described in syntax form,
%   specstring ::= timepoint 'to' timepoint
%   timepoint ::= ID-list ['+'|'-' offset]
%   ID-list ::= ID [',' ID-list]
% Blanks may not be tolerated here, and should be removed from user strings
% before being submitted to this function.  Note that the range of legal
% expressions for 'offset' is determined in this implementation by Matlab's
% str2num function, which is pretty liberal; this should not be relied on.
% Likewise any string of the chars '+', '-', ',' before a timepoint will be
% silently ignored, and this should also not be relied on.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

delimiter = findstr('to', specstring);
if ~ all(size(delimiter))
    error('dg_ParseSegmentSpec:NoDelimiter', ...
        'The segment spec "%s" is missing the keyword "to".', ...
            specstring);
else
	start = specstring(1 : delimiter - 1);
	stop = specstring(delimiter + 2 : length(specstring));
    if length(start) == 0 || length(stop) == 0
        error('dg_ParseSegmentSpec:NoTimeSpec', ...
            'There must be a time spec before and after the word "to".');
    end
	specstruct.start = dg_ParseTimeSpec(start);
	specstruct.stop = dg_ParseTimeSpec(stop);
end