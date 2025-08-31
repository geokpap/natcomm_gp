function [TS, TTL, ES, Hdr] = dg_decodeEvtBits(infile, outfile, ...
    evtstable, varargin)
%dg_decodeEvtBits(infile, outfile, evtstable)
%dg_decodeEvtBits(infile, outfile, evtstable, delay)
% Converts TTL values including hard-wired one-bit-per-device events into
% numeric codes for bit-on and bit-off for each bit.  Reads and writes
% Neuralynx events file format.
%INPUTS
% infile - pathname to input file
% outfile - pathname to output file
% evtstable - determines which bits in the TTL values are monitored for
%   state transitions, and how the transitions are represented in the
%   output.  <evtstable> contains one row for each bit to monitor.  The
%   first column designates the bit to monitor, where bits are numbered
%   sequentially starting at 1 for the lowest-order bit (usually called
%   "bit0" by engineers).  The second column is the event code to use when
%   that bit goes high (turns on), and the third column is the event code
%   to use when that bit goes low (turns off).  If more than one bit
%   changes between two successive event records, then separate events are
%   inserted for each transition at <delay> microsecond intervals.
% delay - the number of microseconds delay between the decoded events for
%   bits that change simultaneously.  This arg is optional and defaults to
%   1.
%OUTPUTS
% TS, TTL, ES, Hdr - modified values as would be returned by dg_readEvents.
%OPTIONS
% 'values', vs - <vs> is a struct containing values for all of the variables
%   normally returned by dg_readEvents; invoking this option bypasses the
%   call to dg_readEvents and also the call to dg_writeNlxEvents, and uses
%   the values in the struct instead.  Each field in the struct must have
%   the name of the corresponding dg_readEvents returned value (i.e., TS,
%   TTL, ES, Hdr).  The values in the TTL field do NOT get converted to be
%   all non-negative (whereas the default is to convert them, rendering the
%   highest order bit useless).

%$Rev: 256 $
%$Date: 2017-06-08 15:14:53 -0400 (Thu, 08 Jun 2017) $
%$Author: dgibson $

opts2delete = [];
vs = [];

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'values'
            opts2delete(end+1) = argnum; %#ok<*AGROW>
            argnum = argnum + 1;
            opts2delete(end+1) = argnum; %#ok<*AGROW>
            vs = varargin{argnum};
        otherwise
            error('funcname:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end
varargin(opts2delete) = [];


if isempty(varargin)
    delay = 1;
else
    delay = varargin{1};
end

if isempty(vs)
    [TS, TTL, ES, Hdr] = dg_readEvents(infile);
    TTL(TTL<0) = TTL(TTL<0) + 2^16;
else
    names = fieldnames(vs);
    for k = 1:length(names);
        eval(sprintf('%s = vs.%s;', names{k}, names{k}));
    end
end

% Find the bits that differ and are in <evtstable>:
bitdiff = NaN(size(TTL));
bitmask = 0;
bitcodes = NaN(16,2);   % rownum = bitnum (starting at 1); col 1 on, 2 off
for k = 1:size(evtstable,1)
    bitmask = bitor(bitmask, 2^(evtstable(k,1) - 1));
    bitcodes(evtstable(k,1), :) = evtstable(k,2:3);
end
bitdiff(1) = bitand(bitmask, TTL(1));
for k = 2:length(TTL)
    bitdiff(1,k) = bitand(bitmask, bitxor(TTL(k), TTL(k-1)));
end

% Incredibly, the following method of counting bits is 3x faster than a
% for-loop bit shifter:
[f,e]=log2(max(bitdiff)); %#ok<*ASGLU>
for k = 1:length(bitdiff)
    numbitsdiff(k) = sum(rem( floor(bitdiff(k) * ...
        pow2(1-max(1,e):0)), 2 ));
end

% Single bit changes only require a new TTL and matching ES:
singlediffidx = find(numbitsdiff == 1);
for k = 1:length(singlediffidx)
    [f,bitnum] = log2(bitdiff(singlediffidx(k)));
    if bitand(bitdiff(singlediffidx(k)), TTL(singlediffidx(k)))
        % the bit turned on
        TTL(singlediffidx(k)) = bitcodes(bitnum, 1);
        ES{singlediffidx(k)} = sprintf('bit(%d) ON', bitnum);
    else
        % the bit turned off
        TTL(singlediffidx(k)) = bitcodes(bitnum, 2);
        ES{singlediffidx(k)} = sprintf('bit(%d) OFF', bitnum);
    end
end

% Multi-bit changes require surgery to the data vectors:
multidiffidx = find(numbitsdiff > 1);
if ~isempty(multidiffidx)
    TS2 = cell(1, 2*length(multidiffidx)+1);
    TTL2 = cell(1, 2*length(multidiffidx)+1);
    ES2 = {};
    for k = 1:length(multidiffidx)
        % copy the preceding series of events
        if k == 1
            idx = 1 : (multidiffidx(k) - 1);
        else
            idx = multidiffidx(k-1)+1 : multidiffidx(k)-1;
        end
        TTL2{2*k - 1} = TTL(idx);
        ES2 = [ES2; ES(idx)]; %#ok<*AGROW>
        TS2{2*k - 1} = TS(idx);
        % create the new series of events to insert
        bitnumlist = find( ...
            rem(floor(bitdiff(multidiffidx(k))*pow2(0:-1:1-max(1,e))),2) );
        for j = 1:length(bitnumlist)
            bitnum = bitnumlist(j);
            if bitand(2^(bitnum-1), TTL(multidiffidx(k)))
                % the bit turned on
                TTL2{2*k}(j) = bitcodes(bitnum, 1);
                ES2 = [ES2; {sprintf('bit(%d) ON', bitnum)}];
            else
                % the bit turned off
                TTL2{2*k}(j) = bitcodes(bitnum, 2);
                ES2 = [ES2; {sprintf('bit(%d) OFF', bitnum)}];
            end
            % timestamps are in microseconds; add <delay> microseconds for
            % each additional event:
            TS2{2*k}(j) = TS(multidiffidx(k)) + delay * (j - 1);
        end
    end
    % copy the events after the last insertion
    idx = multidiffidx(k)+1 : length(TS);
    TTL2{2*k + 1} = TTL(idx);
    ES2 = [ES2; ES(idx)];
    TS2{2*k + 1} = TS(idx);
    % convert cell arrays to vectors (cell2mat does not work on cell arrays
    % containing cell arrays)
    TTL = cell2mat(TTL2);
    TS = cell2mat(TS2);
    ES = ES2;
end

Hdr = [
    Hdr(1)
    {sprintf('## File Name: (dg_decodeEvtBits): %s ', outfile)}
    {sprintf('## Time Written:  (m/d/y): %s ', datestr(now, 0))}
    Hdr(2:end)
    ];
if isempty(vs)
    dg_writeNlxEvents(outfile, TS, TTL, ES, Hdr);
end

