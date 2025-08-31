function dg_insertEvtIDs(infile, outfile, evtstable, varargin)
%dg_insertEvtIDs(infile, outfile)
% Reads the Neuralynx events file <infile> and uses <evtstable> to fill in
% new TTL event IDs to match the event strings.  Does not replace the TTL
% field for any records that do not match any entries in <evtstable>.
%INPUTS
% infile - pathname to input file
% outfile - pathname to output file
% evtstable - a cell array with two columns.  The first column gives the
%   numeric values desired for the TTL event ID.  The second column gives
%   the string value that must EXACTLY match the event string in the events
%   file.  The matching is case-sensitive, whitespace-sensitive, in general
%   sensitive to anything to which Matlab's 'ismember' function is
%   sensitive.
%OPTIONS
% 'strobe' - set the strobe (2^15) bit on the TTL event IDs.

%$Rev: 101 $
%$Date: 2011-04-08 18:18:12 -0400 (Fri, 08 Apr 2011) $
%$Author: dgibson $

opentime = now();
strobeflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'strobe'
            strobeflag = true;
        otherwise
            error('dg_insertEvtIDs:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end
[TS, TTL, ES, Hdr] = dg_readEvents(infile);
for k = 1:size(evtstable,1)
    evtidx = ismember(ES, evtstable{k,2});
    TTL(evtidx) = evtstable{k,1} - strobeflag * 2^15;
end
Hdr{2} = sprintf('## File Name: (dg_insertEvtIDs): %s ', outfile);
Hdr{3} = sprintf('## Time Opened:  (m/d/y): %s ', datestr(opentime, 0));
Hdr{4} = sprintf('## Time Closed:  (m/d/y): %s ', datestr(now, 0));
dg_writeNlxEvents(outfile, TS, TTL, ES, Hdr);