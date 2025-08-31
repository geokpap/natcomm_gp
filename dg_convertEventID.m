function [TS, TTL, ES, Hdr] = dg_convertEventID(filepath, EID, TTLID)
%[TS, TTL, ES, Hdr] = dg_convertEventID(filepath, EID, TTLID)
% Reads a Neurlynx events file and changes the TTL value assigned to events
% that are recorded with a specified "Event ID" value.
%INPUTS
% filepath: Neuralynx events file to read.
% EID: the "Event ID" code that should be converted.
% TTLID: the TTL value that events having <EID> should be given.
% OUTPUTS
%   TS:  Timestamps in clock ticks (units as recorded)
%   TTL:  TTL IDs.  These are returned as signed integers, with 2^15 bit 
%       (strobe/sync) propagated to the left.  To extract the lower 15
%       bits, add 2^15 to negative values to yield positive integers.
%   ES:  Event strings.

%$Rev: 224 $
%$Date: 2015-07-27 14:39:33 -0400 (Mon, 27 Jul 2015) $
%$Author: dgibson $

try
    if ispc
        [TS, EventID, TTL, ES, Hdr] = Nlx2MatEV_411(filepath, ...
            [1, 1, 1, 0, 1], 1, 1, []);
            % Convert to the "bugged" representation of v4.1.2 and earlier:
            TTL(TTL>=2^15) = TTL(TTL>=2^15) - 2^16;
    elseif ismac || isunix
        [TS, EventID, TTL, ES, Hdr] = Nlx2MatEV_v3(filepath, ...
            [1, 1, 1, 0, 1], 1, 1, []);
    else
        error('dg_readEvents:arch', ...
            'Unrecognized computer platform');
    end
catch e
    errmsg = sprintf( ...
        'Nlx2Mat arguments in dg_readEvents:\nfile:%s', ...
        dg_thing2str(filepath) );
    errmsg = sprintf('%s\n%s\n%s', ...
        errmsg, e.identifier, e.message);
    for stackframe = 1:length(e.stack)
        errmsg = sprintf('%s\n%s\nline %d', ...
            errmsg, e.stack(stackframe).file, ...
            e.stack(stackframe).line);
    end
    error('dg_convertEventID:failed', '%s', errmsg);
end
TTL(EventID == EID) = TTLID;

