function specstruct = dg_ParseTimeSpec(specstring)

%DG_PARSETIMESPEC Accepts a string and returns a structure representing
%a time specification.  Throws an error if parsing fails.
% See dg_ParseSegmentSpec for syntax.  Blanks are NOT tolerated here.
% The structure specstruct is:
% specstruct.idlist - an array of acceptable event TTL ID's
% specstruct.offset - offset in mS

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

idlist = [];
offset = 0;
delimiter = ',';
while delimiter == ','
	[token, specstring] = strtok(specstring, ',+-');
    id = str2num(token);
    if ~ all(size(id))
        error('dg_ParseTimeSpec:BadNumber', 'Expected a number at "%s"', ...
            token);
    end
    if id ~= floor(id)
        error('dg_ParseTimeSpec:NonInteger', ...
            'Event ID %s should be an integer', token);
    end
	idlist = [idlist; id];
    if ~ all(size(specstring))
        delimiter = [];
    else
        delimiter = specstring(1);
    end
end
if all(size(specstring))
    offset = str2num(specstring(2:length(specstring)));
    if delimiter == '-'
        offset = - offset;
    end
end
specstruct.idlist = idlist;
specstruct.offset = offset;