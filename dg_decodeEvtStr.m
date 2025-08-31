function dg_decodeEvtStr(infile, outfile, convtable)
%INPUTS
% infile: relative or absolute pathname to a Neuralynx events file.
% outfile: relative or absolute pathname at which to create a new Neuralynx
%   events file.  Any pre-existing file at that location will be
%   overwritten.
% convtable: a two-column cell array where the first column contains TTL
%   codes (positive integers) and the second column contains regular
%   expressions which match the event strings in <filename>.
%OUTPUTS
% Output is sent to the specified output file.  It is a copy of the same
% data in the input file, but with TTLs completely replaced by the numeric
% values specified in <convtable>.

%$Rev: 200 $
%$Date: 2014-06-24 19:31:12 -0400 (Tue, 24 Jun 2014) $
%$Author: dgibson $

[TS, TTL, ES, Hdr] = dg_readEvents(infile);
Hdr = [
    Hdr(1:4)
    {'## Converted by dg_decodeEvtStr '}
    Hdr(5:end)
    ];
unmatched = true(size(ES));
for convidx = 1:size(convtable, 1)
    ismatch = ~cellfun(@isempty, regexp(ES, convtable{convidx, 2}));
    TTL(ismatch) = convtable{convidx, 1};
    unmatched(ismatch) = false;
end
if any(unmatched)
    unmatchedstr = unique(ES(unmatched));
    unmatchedmultiline = unmatchedstr{1};
    for k = 2:length(unmatchedstr)
        unmatchedmultiline = sprintf( '%s\n%s', unmatchedmultiline, ...
            unmatchedstr{k} );
    end
    warning('dg_decodeEvtStr:nomatch', ...
        'The following event strings have no match in <convtable> and have not been decoded:\n%s', ...
        unmatchedmultiline );
end
dg_writeNlxEvents(outfile, TS, TTL, ES, Hdr);

